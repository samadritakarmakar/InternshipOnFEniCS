// Copyright (C) 2009-2012 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.

#ifndef __UFL_QUADRATURE_FUNCTION_H
#define __UFL_QUADRATURE_FUNCTION_H

#include <memory>
#include <vector>
#include <Eigen/Dense>
//#include <boost/multi_array.hpp>

#include <dolfin/function/GenericFunction.h>
#include <dolfin/fem/Form.h>

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

  // Class to compute UFL expression at integration point values

  class UFLQuadratureFunction : public dolfin::GenericFunction
  {
  public:

    /// Delete copy constructor and assignement
    // UFLQuadratureFunction(UFLQuadratureFunction&&) = delete; // shorter than deleting assignment operator
    UFLQuadratureFunction& operator=(const UFLQuadratureFunction&) = delete;  // Disallow copying
    UFLQuadratureFunction(const UFLQuadratureFunction&) = delete;

    /// Constructor
    UFLQuadratureFunction(std::shared_ptr<const dolfin::Form> ufl_form,
                          std::shared_ptr<const dolfin::Form> ufl_reference_form);


    /// Destructor
    ~UFLQuadratureFunction();

    std::shared_ptr<const dolfin::FunctionSpace> function_space() const
    {
      return _ufl_form->function_space(0);
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

    std::shared_ptr<const dolfin::FiniteElement> element() const
    { return _element; }

  private:
  
    // Form
    std::shared_ptr<const dolfin::Form> _ufl_form;

    // Finite element
    std::shared_ptr<const dolfin::FiniteElement> _element;

    //std::shared_ptr<StateUpdate> _state_updater;

    // Num 'IP' dofs per IP
    // const std::size_t num_ip_dofs;
    
    // Number of quadrature points per cell
    // const std::size_t num_ip;

    // Coefficient
    mutable Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _A_e;

    // Coefficient on reference cell
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _A_e_reference;

    // Num 'IP' dofs per IP
    const std::size_t num_ip_dofs;
    
    // Number of quadrature points per cell
    const std::size_t num_ip;
  };

}

#endif
