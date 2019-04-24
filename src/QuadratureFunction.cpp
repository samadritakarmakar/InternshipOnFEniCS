// Copyright (C) 2009-2012 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2009-10-02
// Last changed: 2012-07-17

#include <dolfin/common/Timer.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include "ConstitutiveUpdate.h"
#include "StateUpdate.h"
#include "QuadratureFunction.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
QuadratureFunction::QuadratureFunction(std::shared_ptr<const dolfin::Mesh> mesh,
                       std::shared_ptr<const dolfin::FiniteElement> element,
                       const std::vector<double>& w)
  : _element(element), _state_updater(0), 
  num_ip_dofs(element->value_dimension(0)), num_ip(element->space_dimension()/element->value_dimension(0)),
  _w(w)
{
  dolfin_assert(element);

  // FIXME: Check that we have a quadrature element

  // Num 'IP' dofs per IP
  //const std::size_t num_ip_dofs = element->value_dimension(0);

  // Number of quadrature points per cell
  //const std::size_t num_ip = element->space_dimension()/num_ip_dofs;

  // Allocate space to hold data
  //_data.resize(boost::extents[mesh->num_cells()][num_ip][num_ip_dofs]);
  //std::fill(_data.data(), _data.data() + _data.num_elements(), 0.0);
}
//-----------------------------------------------------------------------------
QuadratureFunction::QuadratureFunction(std::shared_ptr<const dolfin::Mesh> mesh,
                       std::shared_ptr<const dolfin::FiniteElement> element,
                       std::shared_ptr<StateUpdate> state_updater,
                       const std::vector<double>& w)
  : _element(element), _state_updater(state_updater), 
  num_ip_dofs(element->value_dimension(0)), num_ip(element->space_dimension()/element->value_dimension(0)),
  _w(w)
{
  dolfin_assert(element);

  // FIXME: Check that we have a qudarature element

  // Num 'IP' dofs per IP
  //const std::size_t num_ip_dofs = element->value_dimension(0);

  // Number of quadrature points per cell
  //const std::size_t num_ip = element->space_dimension()/num_ip_dofs;

  // Allocate space to hold data
  //_data.resize(boost::extents[mesh->num_cells()][num_ip][num_ip_dofs]);
  //std::fill(_data.data(), _data.data() + _data.num_elements(), 0.0);
}
//-----------------------------------------------------------------------------
QuadratureFunction::~QuadratureFunction()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void QuadratureFunction::restrict(double* w,
                                  const dolfin::FiniteElement& elemenet,
                                  const dolfin::Cell& cell,
                                  const double* vertex_coordinates,
                                  const ufc::cell& ufc_cell) const
{
  if (_state_updater)
    _state_updater->update(cell, vertex_coordinates);
  std::copy(_w.begin(), _w.end(), w);
}
//-----------------------------------------------------------------------------
void QuadratureFunction::compute_vertex_values(std::vector<double>&,
                                               const dolfin::Mesh&) const
{
  dolfin::error("QuadratureFunction::compute_vertex_values not implemented");
}
//-----------------------------------------------------------------------------
// void QuadratureFunction::compute_mean(dolfin::MeshFunction<double>& mf) const
// {
//   if (_data.shape()[0] != mf.size())
//   {
//     dolfin::error("Data/MeshFunction size mis-match in QuadratureFunction::compute_mean");
//   }
// 
//   if (_data.shape()[2] != 1)
//     dolfin::error("QuadratureFunction::compute_mean only supports scalar data");
// 
//   dolfin_assert(mf.mesh());
//   const dolfin::Mesh& mesh = *mf.mesh();
//   if (mf.size() != mesh.num_cells())
//     dolfin::error("Size mis-match in QuadratureFunction::compute_mean");
// 
//   const std::size_t num_ip = _data.shape()[1];
//   for (dolfin::CellIterator cell(mesh); !cell.end(); ++cell)
//   {
//     double x = 0.0;
//     const std::size_t index = cell->index();
//     for (std::size_t p = 0; p < num_ip; ++p)
//       x += _data[index][p][0];
//     mf[index] = x/num_ip;
//   }
// }
//-----------------------------------------------------------------------------
