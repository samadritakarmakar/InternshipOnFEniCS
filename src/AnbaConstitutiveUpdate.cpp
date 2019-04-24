// Copyright (C) 2006-2017 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.

#include <vector>
#include <ufc.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Cell.h>
#include "PlasticityModel.h"
#include "utils.h"
#include "AnbaConstitutiveUpdate.h"
#include "UFLQuadratureFunction.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
AnbaConstitutiveUpdate::AnbaConstitutiveUpdate(
  std::shared_ptr<const UFLQuadratureFunction> strain_func,
  std::shared_ptr<const UFLQuadratureFunction> strain_func_z,
  std::shared_ptr<const UFLQuadratureFunction> strain_func_zz,
  std::shared_ptr<const PlasticityModel> plastic_model)
  : _strain_func(strain_func),
    _strain_func_z(strain_func_z),
    _strain_func_zz(strain_func_zz),
    _plastic_model(plastic_model),
    _De(_plastic_model->elastic_tangent())
{
  _eps_p = std::make_shared<HistoryData>(_strain_func->function_space()->mesh(), _strain_func->element(), 6);
  _eps_p_equiv = std::make_shared<HistoryData>(_strain_func->function_space()->mesh(), _strain_func->element(), 1);
  _eps_p_equiv_z = std::make_shared<HistoryData>(_strain_func->function_space()->mesh(), _strain_func->element(), 1);
  _eps_p_equiv_zz = std::make_shared<HistoryData>(_strain_func->function_space()->mesh(), _strain_func->element(), 1);

  // Get stress UFC element
  auto ufc_element_sigma = _strain_func->element()->ufc_element();
  dolfin_assert(ufc_element_sigma);

  // Get stress dof dimension data
  const std::size_t dim = ufc_element_sigma->space_dimension();

  // Compute number of quadrature points per cell
  _num_ip_dofs = _strain_func->element()->value_dimension(0);
  _num_ip_per_cell = _strain_func->element()->space_dimension()/_num_ip_dofs;

  // Resize _w_stress/tangent
  _prev_w_strain.resize(_num_ip_dofs*_num_ip_per_cell);
  _prev_w_strain_z.resize(_num_ip_dofs*_num_ip_per_cell);
  _prev_w_strain_zz.resize(_num_ip_dofs*_num_ip_per_cell);
  _w_strain.resize(_num_ip_dofs*_num_ip_per_cell);
  _w_strain_z.resize(_num_ip_dofs*_num_ip_per_cell);
  _w_strain_zz.resize(_num_ip_dofs*_num_ip_per_cell);
  _w_stress.resize(_num_ip_dofs*_num_ip_per_cell);
  _w_stress_z.resize(_num_ip_dofs*_num_ip_per_cell);
  _w_stress_zz.resize(_num_ip_dofs*_num_ip_per_cell);
  _w_tangent.resize(_num_ip_dofs*_num_ip_dofs*_num_ip_per_cell);
  _w_tangent_sz_ez.resize(_num_ip_dofs*_num_ip_dofs*_num_ip_per_cell);
  _w_tangent_szz_ez.resize(_num_ip_dofs*_num_ip_dofs*_num_ip_per_cell);
  _w_tangent_szz_ezz.resize(_num_ip_dofs*_num_ip_dofs*_num_ip_per_cell);

  const std::size_t num_cells = strain_func->function_space()->mesh()->num_cells();
  _plastic_last.resize(boost::extents[num_cells][_num_ip_per_cell]);
  std::fill(_plastic_last.data(),
            _plastic_last.data() + _plastic_last.num_elements(),
            false);
}
//-----------------------------------------------------------------------------
AnbaConstitutiveUpdate::~AnbaConstitutiveUpdate()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void AnbaConstitutiveUpdate::update(const dolfin::Cell& cell,
                                const double* vertex_coordinates)
{
  const std::size_t cell_index = cell.index();

  ufc::cell ufc_cell;
  cell.get_cell_data(ufc_cell);
  _strain_func->restrict(&_w_strain[0], *_strain_func->element(), cell, vertex_coordinates, ufc_cell);
  _strain_func_z->restrict(&_w_strain_z[0], *_strain_func_z->element(), cell, vertex_coordinates, ufc_cell);
  _strain_func_zz->restrict(&_w_strain_zz[0], *_strain_func_zz->element(), cell, vertex_coordinates, ufc_cell);

  // check for caching
  if (_prev_cell_index != cell_index || 
        _w_strain != _prev_w_strain ||
        _w_strain_z != _prev_w_strain_z ||
        _w_strain_zz != _prev_w_strain_zz
     ) {
    _prev_cell_index = cell_index;
    _prev_w_strain = _w_strain;
    _prev_w_strain_z = _w_strain_z;
    _prev_w_strain_zz = _w_strain_zz;

    Eigen::Matrix<double, 6, 6> cons_tangent;
    Eigen::Matrix<double, 6, 1> strain, strain_z, strain_zz, strain_p, trial_stress;
    Eigen::Matrix<double, 1, 1> strain_p_eq;
    strain.setZero();
    strain_p.setZero();
    trial_stress.setZero();
    
    // Loop over quadrature points
    for (std::size_t ip = 0; ip < _num_ip_per_cell; ip++)
    {
      // Assume elastic tangent
      cons_tangent = _De;

      // Get plastic strain from previous converged time step
      _eps_p->get_old_values(cell_index, ip, strain_p);

      // Compute strain on physical cell (Voigt notation)
      for (std::size_t d = 0; d < _num_ip_dofs; ++d) {
  	  strain(d) = _w_strain[_num_ip_per_cell*d  + ip];
  	  strain_z(d) = _w_strain_z[_num_ip_per_cell*d  + ip];
  	  strain_zz(d) = _w_strain_zz[_num_ip_per_cell*d  + ip];
      }

      // Compute trial stress
      trial_stress = _De*(strain - strain_p);

      // Get equivalent plastic strain from previous converged time
      // step
      _eps_p_equiv->get_old_values(cell_index, ip, strain_p_eq);

      // Testing trial stresses, if yielding occurs the stresses are
      // mapped back onto the yield surface, and the updated parameters
      // are returned.
      const bool active = _plastic_last[cell_index][ip];
      auto ret_map = return_mapping.closest_point_projection(_plastic_model, cons_tangent,
  					      trial_stress, strain_p,
  					      strain_p_eq(0),
  					      active);
      _plastic_last[cell_index][ip] = false;

      // Update plastic strain history for current load step
      _eps_p->set_new_values(cell_index, ip, strain_p);

      // Update equivalent plastic strain for current load step
      _eps_p_equiv->set_new_values(cell_index, ip, strain_p_eq);
      
  	  

      // Copy data into structures
      for (std::size_t d = 0; d < _num_ip_dofs; ++d)
  	  _w_stress[_num_ip_per_cell*d  + ip] = trial_stress[d];
      for (std::size_t d0 = 0; d0 < _num_ip_dofs; ++d0)
      {
  	for (std::size_t d1 = 0; d1 < _num_ip_dofs; ++d1)
  	{
  	  const std::size_t pos = d0*_num_ip_dofs + d1;
  	  _w_tangent[_num_ip_per_cell*pos  + ip] = cons_tangent(d0, d1);
  	}
      }
      
      Eigen::Matrix<double, 6, 6> w_tangent_sz_ez = _De;
      Eigen::Matrix<double, 6, 6> w_tangent_szz_ez; w_tangent_szz_ez.setZero();
      Eigen::Matrix<double, 6, 6> w_tangent_szz_ezz = _De;
      Eigen::Matrix<double, 6, 6> A_tilde; A_tilde.setZero();
      Eigen::Matrix<double, 6, 1> sigma_z;// = _De * strain_z;
      Eigen::Matrix<double, 6, 1> sigma_zz;// = _De * strain_zz;
      if (ret_map.first) {
  	  Eigen::Matrix<double, 6, 1> df_dsigma;
  	  Eigen::Matrix<double, 6, 6> ddf_ddsigma;
  	  _plastic_model->df(df_dsigma, trial_stress);
  	  _plastic_model->ddg(ddf_ddsigma, trial_stress);
  	  Eigen::Matrix<double, 6, 1> E_df_dsigma = _De * df_dsigma;
  	  Eigen::Matrix<double, 6, 1> B = E_df_dsigma * 2.;
  	  double A = df_dsigma.dot(B) + _plastic_model->hardening_parameter(strain_p_eq(0));
  	  double lambda_z = B.dot(strain_z) / A;
  	  Eigen::Matrix<double, 6, 6> tmp = E_df_dsigma * B.transpose() / A;
  	  w_tangent_sz_ez -= tmp;
  	  w_tangent_szz_ezz -= tmp;
  	  A_tilde = w_tangent_sz_ez * ddf_ddsigma * w_tangent_sz_ez * 4.;
  	  A_tilde -= ((ddf_ddsigma * E_df_dsigma).transpose() * w_tangent_sz_ez).transpose() * B.transpose()  / A * 4.;
  	  Eigen::Matrix<double, 6, 1> tmpv = A_tilde * strain_z;
  	  w_tangent_szz_ez = E_df_dsigma * tmpv.transpose();
      }
      sigma_z = w_tangent_sz_ez * strain_z;
      sigma_zz = w_tangent_szz_ezz * strain_zz + w_tangent_szz_ez * strain_z;
      for (std::size_t d = 0; d < _num_ip_dofs; ++d) {
  	  _w_stress_z[_num_ip_per_cell*d  + ip] = sigma_z[d];
  	  _w_stress_zz[_num_ip_per_cell*d  + ip] = sigma_zz[d];
      }
      for (std::size_t d0 = 0; d0 < _num_ip_dofs; ++d0)
      {
  	for (std::size_t d1 = 0; d1 < _num_ip_dofs; ++d1)
  	{
  	  const std::size_t pos = d0*_num_ip_dofs + d1;
  	  _w_tangent_sz_ez[_num_ip_per_cell*pos  + ip] = w_tangent_sz_ez(d0, d1);
  	  _w_tangent_szz_ez[_num_ip_per_cell*pos  + ip] = w_tangent_szz_ez(d0, d1);
  	  _w_tangent_szz_ezz[_num_ip_per_cell*pos  + ip] = w_tangent_szz_ezz(d0, d1);
  	}
      }
    }
  }
}
//-----------------------------------------------------------------------------
void AnbaConstitutiveUpdate::update_history()
{
  // Update plastic elements
  const boost::multi_array<double, 3>& old_eps = _eps_p_equiv->old_data();
  const boost::multi_array<double, 3>& new_eps = _eps_p_equiv->current_data();

  const std::size_t num_cells = _plastic_last.shape()[0];
  const std::size_t ip_per_cell = _plastic_last.shape()[1];
  dolfin_assert(old_eps.shape()[0] == num_cells);

  for (std::size_t c = 0; c < num_cells; ++c)
  {
    for (std::size_t p = 0; p < ip_per_cell; ++p)
    {
      if ((new_eps[c][p][0] - old_eps[c][p][0] > 0.0))
        _plastic_last[c][p] = true;
      else
        _plastic_last[c][p] = false;
    }
  }

  _eps_p->update_history();
  _eps_p_equiv->update_history();
}
//-----------------------------------------------------------------------------
