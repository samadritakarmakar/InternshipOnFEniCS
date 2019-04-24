// Copyright (C) 2006-2011 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2011-02-06

#ifndef __PLASTICITY_PROBLEM_H
#define __PLASTICITY_PROBLEM_H

#include <memory>
#include <vector>
#include <dolfin/nls/NonlinearProblem.h>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/fem/SystemAssembler.h>
#include "HistoryData.h"
#include "PlasticityModel.h"
#include "ReturnMapping.h"

namespace dolfin
{
  class DirichletBC;
  class GenericDofMap;
  class GenericMatrix;
  class GenericVector;
  class Form;
  class Function;
  class NewtonSolver;
}

namespace fenicssolid
{

  class QuadratureFunction;

  class PlasticityProblem : public dolfin::NonlinearProblem
  {
  public:

    /// Delete copy constructor and assignement
    PlasticityProblem& operator=(const PlasticityProblem&) = delete;  // Disallow copying
    PlasticityProblem(const PlasticityProblem&) = delete;

    /// Constructor
    PlasticityProblem(std::shared_ptr<const dolfin::Form> a,
                      std::shared_ptr<const dolfin::Form> L,
                      std::shared_ptr<dolfin::Function> U,
                      std::shared_ptr<QuadratureFunction> tangent,
                      std::shared_ptr<QuadratureFunction> sigma,
                      const std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs);

    /// Constructor
    PlasticityProblem(std::shared_ptr<const dolfin::Form> a,
                      std::shared_ptr<const dolfin::Form> L,
                      std::shared_ptr<dolfin::Function> U,
                      std::shared_ptr<QuadratureFunction> tangent,
                      std::shared_ptr<QuadratureFunction> sigma,
                      const std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs,
                      std::shared_ptr<dolfin::NewtonSolver> newton_solver);

    /// Destructor
    ~PlasticityProblem();

    /// Loop quadrature points and compute local tangents and stresses
    void form(dolfin::GenericMatrix& A, dolfin::GenericMatrix& P,
              dolfin::GenericVector& b, const dolfin::GenericVector& x);

    /// User defined assemble of residual vector
    void F(dolfin::GenericVector& b, const dolfin::GenericVector& x);

    /// User defined assemble of Jacobian matrix
    void J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x);

    /// Update history variables
    void update_history();

    /// Public data members that might be of interest after the Newton
    /// solver has converged
    HistoryData eps_p_data;
    HistoryData eps_p_eq_data;

    /// Parameters
    dolfin::Parameters parameters;

  private:

    // For system after constitutive update
    void form_tensors(dolfin::GenericMatrix& A, dolfin::GenericVector& b,
                      const dolfin::GenericVector& x);

    // Assembler
    dolfin::SystemAssembler assembler;

    // Helper data structures to handle variable updates on
    // integration points
    //std::shared_ptr<QuadratureFunction> _tangent;
    //std::shared_ptr<QuadratureFunction> _sigma;

    // Displacement
    std::shared_ptr<const dolfin::Function> _u;

    // The Newton solver
    std::shared_ptr<dolfin::NewtonSolver> _newton_solver;

  };
}

#endif
