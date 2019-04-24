// Copyright (C) 2009-2012 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2009-10-02
// Last changed: 2012-07-17

#include <dolfin/common/utils.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include "QuadratureData.h"
#include "UFLReturnMapping.h"

using namespace fenicssolid;
//-----------------------------------------------------------------------------
QuadratureData::QuadratureData(std::shared_ptr<const dolfin::Mesh> mesh,
                         std::shared_ptr<const dolfin::FiniteElement> element,
                         std::size_t size,
                         std::shared_ptr<UFLReturnMapping> state_updater
) : _state_updater(state_updater), recursive_guard(false)
{
  // Number of cells
  const std::size_t num_cells = mesh->num_cells();

  // Num 'IP' dofs per IP
  const std::size_t num_ip_dofs = element->value_dimension(0);

  // Number of quadrature points per cell
  const std::size_t num_ip = element->space_dimension()/num_ip_dofs;

  _cur_vals.resize(boost::extents[num_cells][num_ip][size]);
  std::fill(_cur_vals.data(), _cur_vals.data() + _cur_vals.num_elements(),
            0.0);

  // FIXME: Test assumptions.
}

void QuadratureData::restrict(double* w,
	      const dolfin::FiniteElement& element,
	      const dolfin::Cell& cell,
	      const double* vertex_coordinates,
	      const ufc::cell& ufc_cell) const
{
  if (_state_updater && (!recursive_guard)) {
    recursive_guard = true;
    _state_updater->closest_point_projection(this, cell);
    recursive_guard = false;
  }

  const std::size_t cell_index = cell.index();
  const std::size_t num_ip = _cur_vals.shape()[1];
  const std::size_t num_ip_dofs = _cur_vals.shape()[2];
  for (std::size_t ip = 0; ip < num_ip; ip++)
    for (std::size_t d = 0; d < num_ip_dofs; d++)
      w[num_ip*d + ip] = _cur_vals[cell_index][ip][d];
}

//-----------------------------------------------------------------------------
