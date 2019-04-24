// Copyright (C) 2009-2012 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.

#ifndef __QUADRATURE_FUNCTION_H
#define __QUADRATURE_FUNCTION_H

#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include <dolfin/function/GenericFunction.h>

namespace dolfin
{
  class Cell;
  class FiniteElement;
  class Function;
  class Mesh;
  template<typename T> class MeshFunction;
}

namespace ufc
{
  class cell;
}

namespace fenicssolid
{

  class StateUpdate;

  // Class to ease the handling of integration point values

  class QuadratureFunction : public dolfin::GenericFunction
  {
  public:

    /// Delete copy constructor and assignement
    // QuadratureFunction(QuadratureFunction&&) = delete; // shorter than deleting assignment operator
    QuadratureFunction& operator=(const QuadratureFunction&) = delete;  // Disallow copying
    QuadratureFunction(const QuadratureFunction&) = delete;

    /// Constructor
    QuadratureFunction(std::shared_ptr<const dolfin::Mesh> mesh,
                       std::shared_ptr<const dolfin::FiniteElement> element,
                       const std::vector<double>& w);

    /// Constructor
    QuadratureFunction(std::shared_ptr<const dolfin::Mesh> mesh,
                       std::shared_ptr<const dolfin::FiniteElement> element,
                       std::shared_ptr<StateUpdate> state_updater,
                       const std::vector<double>& w);

    /// Destructor
    ~QuadratureFunction();

    std::shared_ptr<const dolfin::FunctionSpace> function_space() const
    {
      return std::shared_ptr<const dolfin::FunctionSpace>(NULL);
    }

    // -- GenericFunction interface

    std::size_t value_rank() const
    { return _element->value_rank(); }

    std::size_t value_dimension(std::size_t i) const
    { return _element->value_dimension(i); }

    std::vector<std::size_t> value_shape() const
    {
      std::vector<std::size_t> _shape(this->value_rank(), 1);
      for (std::size_t i = 0; i < _shape.size(); ++i)
        _shape[i] = this->value_dimension(i);
      return _shape;
    }

    void restrict(double* w,
                  const dolfin::FiniteElement& element,
                  const dolfin::Cell& cell,
                  const double* coordinates,
                  const ufc::cell& ufc_cell) const;

    void compute_vertex_values(std::vector<double>&,
                               const dolfin::Mesh&) const;

    // -- QuadratureFunction extensions

    /// Compute average value per cell for scalar data
//     void compute_mean(dolfin::MeshFunction<double>& mf) const;

    std::shared_ptr<const dolfin::FiniteElement> element() const
    { return _element; }

  private:

    // Finite element
    std::shared_ptr<const dolfin::FiniteElement> _element;

    std::shared_ptr<StateUpdate> _state_updater;

    // Data storage
    //boost::multi_array<double, 3> _data;

    // Num 'IP' dofs per IP
    const std::size_t num_ip_dofs;
    
    // Number of quadrature points per cell
    const std::size_t num_ip;

    // Coefficient
    const std::vector<double>& _w;

  };

}

#endif
