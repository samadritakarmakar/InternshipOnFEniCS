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
#include <dolfin/fem/assemble_local.h>
#include <dolfin/function/FunctionSpace.h>
//#include "ConstitutiveUpdate.h"
//#include "StateUpdate.h"
#include "UFLQuadratureFunction.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
UFLQuadratureFunction::UFLQuadratureFunction(
                          std::shared_ptr<const dolfin::Form> ufl_form,
                          std::shared_ptr<const dolfin::Form> ufl_reference_form)
  : _ufl_form(ufl_form),
  _element(ufl_form->function_space(0)->element()),
  num_ip_dofs(_element->value_dimension(0)), num_ip(_element->space_dimension()/_element->value_dimension(0))
{
  dolfin_assert(_element);
  auto mesh = ufl_reference_form->mesh();
  dolfin::CellIterator cell(*mesh);
  //for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)
  assemble_local(_A_e_reference, *ufl_reference_form, *cell);

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
UFLQuadratureFunction::~UFLQuadratureFunction()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void UFLQuadratureFunction::restrict(double* w,
                                  const dolfin::FiniteElement& element,
                                  const dolfin::Cell& cell,
                                  const double* vertex_coordinates,
                                  const ufc::cell& ufc_cell) const
{
  assemble_local(_A_e, *_ufl_form, cell);
  _A_e.array() /= _A_e_reference.array();
  Eigen::Map<Eigen::MatrixXd>(w, _A_e.rows(), _A_e.cols()) = _A_e;
}
//-----------------------------------------------------------------------------
void UFLQuadratureFunction::compute_vertex_values(std::vector<double>&,
                                               const dolfin::Mesh&) const
{
  dolfin::error("QuadratureFunction::compute_vertex_values not implemented");
}
