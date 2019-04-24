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
#include "HistoryData.h"

using namespace fenicssolid;
#include<iostream>
//-----------------------------------------------------------------------------
HistoryData::HistoryData(std::shared_ptr<const dolfin::Mesh> mesh,
                         std::shared_ptr<const dolfin::FiniteElement> element,
                         std::size_t size)
{
  // Number of cells
  const std::size_t num_cells = mesh->num_cells();

  // Num 'IP' dofs per IP
  const std::size_t num_ip_dofs = element->value_dimension(0);

  // Number of quadrature points per cell
  const std::size_t num_ip = element->space_dimension()/num_ip_dofs;

  _old_vals.resize(boost::extents[num_cells][num_ip][size]);
  _cur_vals.resize(boost::extents[num_cells][num_ip][size]);
  std::fill(_old_vals.data(), _old_vals.data() + _old_vals.num_elements(),
            0.0);
  std::fill(_cur_vals.data(), _cur_vals.data() + _cur_vals.num_elements(),
            0.0);

  // FIXME: Test assumptions.
}
//-----------------------------------------------------------------------------
void HistoryData::update_history()
{
  _old_vals = _cur_vals;
}
//-----------------------------------------------------------------------------
void HistoryData::compute_mean(std::shared_ptr<dolfin::MeshFunction<double>> mf) const
{
  if (_old_vals.shape()[0] != mf->size())
  {
   dolfin::error("Data/MeshFunction size mis-match in QuadratureFunction::compute_mean");
  }

  if (_old_vals.shape()[2] != 1)
    dolfin::error("QuadratureFunction::compute_mean only supports scalar data");

  dolfin_assert(mf->mesh());
  const dolfin::Mesh& mesh = *mf->mesh();
  if (mf->size() != mesh.num_cells())
    dolfin::error("Size mis-match in QuadratureFunction::compute_mean");

  const std::size_t num_ip = _old_vals.shape()[1];
  for (dolfin::CellIterator cell(mesh); !cell.end(); ++cell)
  {
    double x = 0.0;
    const std::size_t index = cell->index();
    for (std::size_t p = 0; p < num_ip; ++p)
      x += _old_vals[index][p][0];
    (*mf)[index] = x/num_ip;
  }
}
//-----------------------------------------------------------------------------
std::size_t HistoryData::hash_old() const
{
  return hash(_old_vals);
}
//-----------------------------------------------------------------------------
std::size_t HistoryData::hash_current() const
{
  return hash(_cur_vals);
}
//-----------------------------------------------------------------------------
std::size_t HistoryData::hash(const boost::multi_array<double, 3>& data)
{
  const std::vector<double> tmp(data.origin(), data.origin() + data.size());
  return dolfin::hash_local(tmp);
}
//-----------------------------------------------------------------------------
