// Copyright (C) 2009-2012 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2009-10-02
// Last changed: 2012-07-17

#ifndef __QUADRATURE_DATA_H
#define __QUADRATURE_DATA_H

#include <Eigen/Dense>
#include <boost/multi_array.hpp>
#include <dolfin/function/GenericFunction.h>

// remove after not inlining restric
#include <dolfin/mesh/Cell.h>

namespace dolfin
{
  class FiniteElement;
  class GenericVector;
  class Mesh;
  template<typename T> class MeshFunction;
}

namespace fenicssolid
{

class UFLReturnMapping;

class QuadratureData : public dolfin::GenericFunction
{
public:

  /// Delete copy constructor and assignement
  QuadratureData& operator=(const QuadratureData&) = delete;  // Disallow copying
  QuadratureData(const QuadratureData&) = delete;

  /// Constructor
  QuadratureData(std::shared_ptr<const dolfin::Mesh> mesh,
              std::shared_ptr<const dolfin::FiniteElement> element,
              const std::size_t,
              std::shared_ptr<UFLReturnMapping> state_updater = 0
  );

  /// Link to state updater
  void set_state_updater(std::shared_ptr<UFLReturnMapping> state_updater) {
    _state_updater = state_updater;
  }

  /// Set current values equal to current local ip values
  void set_new_values(std::shared_ptr<QuadratureData> source) 
  {
    // Update values in current vector with local ip values.
    const std::size_t num_cells = _cur_vals.shape()[0];
    const std::size_t num_ip = _cur_vals.shape()[1];
    const std::size_t num_ip_dofs = _cur_vals.shape()[2];
    for (std::size_t cell = 0; cell < num_cells; cell++)
      for (std::size_t ip = 0; ip < num_ip; ip++)
        for (std::size_t d = 0; d < num_ip_dofs; d++)
          _cur_vals[cell][ip][d] = source->_cur_vals[cell][ip][d];
  }

//   template <typename T>
//   void set_new_values(std::size_t cell, unsigned int ip,
//                       const Eigen::MatrixBase<T>& ip_values)
//   {
//     // Update values in current vector with local ip values.
//     const std::size_t value_dim = Eigen::MatrixBase<T>::RowsAtCompileTime;
//     dolfin_assert(ip_values.cols() == 1);
//     dolfin_assert(cell <_cur_vals.shape()[0]);
//     dolfin_assert(ip <_cur_vals.shape()[1]);
//     dolfin_assert(value_dim == _cur_vals.shape()[2]);
// 
//     // Copy data from ip_values
//     for (unsigned int i = 0; i < Eigen::MatrixBase<T>::RowsAtCompileTime; i++)
//       _cur_vals[cell][ip][i] = ip_values(i);
//   }

  std::shared_ptr<const dolfin::FunctionSpace> function_space() const
  {
    return std::shared_ptr<const dolfin::FunctionSpace>(NULL);
  }

    // -- GenericFunction interface

  std::size_t value_rank() const
  {
    if (_cur_vals.shape()[2] == 1) 
      return 0;
    else
      return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return _cur_vals.shape()[2];
  }

  std::vector<std::size_t> value_shape() const
  {
    std::vector<std::size_t> _shape(this->value_rank(), 1);
    if (value_rank())
      _shape[0] = _cur_vals.shape()[2];
    return _shape;
  }

  void restrict(double* w,
                const dolfin::FiniteElement& element,
                const dolfin::Cell& cell,
                const double* vertex_coordinates,
                const ufc::cell& ufc_cell) const;

  void set_cell_values(double* w,
                const dolfin::Cell& cell)
  {
    const std::size_t cell_index = cell.index();
    const std::size_t num_ip = _cur_vals.shape()[1];
    const std::size_t num_ip_dofs = _cur_vals.shape()[2];
    for (std::size_t ip = 0; ip < num_ip; ip++)
      for (std::size_t d = 0; d < num_ip_dofs; d++)
        _cur_vals[cell_index][ip][d] = w[num_ip*d + ip];
  }

  void add_to_cell_values(double* w,
                const dolfin::Cell& cell)
  {
    const std::size_t cell_index = cell.index();
    const std::size_t num_ip = _cur_vals.shape()[1];
    const std::size_t num_ip_dofs = _cur_vals.shape()[2];
    for (std::size_t ip = 0; ip < num_ip; ip++)
      for (std::size_t d = 0; d < num_ip_dofs; d++)
        _cur_vals[cell_index][ip][d] += w[num_ip*d + ip];
  }

  void compute_vertex_values(std::vector<double>&,
                             const dolfin::Mesh&) const
  {
    dolfin::error("QuadratureData::compute_vertex_values not implemented");
  }

  const boost::multi_array<double, 3>& current_data() const
  { return _cur_vals; }

private:

  // State updater
  std::shared_ptr<const UFLReturnMapping> _state_updater;

  // Hash multi_array data
  static std::size_t hash(const boost::multi_array<double, 3>& data);

  // Array of old and current values at integration points
  boost::multi_array<double, 3> _cur_vals;
  
  mutable bool recursive_guard;
 
};

}

#endif
