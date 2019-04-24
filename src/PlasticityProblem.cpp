// Copyright (C) 2006-2011 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2011-02-06

#include <dolfin/common/Timer.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/SystemAssembler.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericLinearSolver.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/nls/NewtonSolver.h>

#include "QuadratureFunction.h"
#include "utils.h"
#include "PlasticityProblem.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
PlasticityProblem::PlasticityProblem(std::shared_ptr<const dolfin::Form> a,
                                     std::shared_ptr<const dolfin::Form> L,
                                     std::shared_ptr<dolfin::Function> u,
                                     std::shared_ptr<QuadratureFunction> tangent,
                                     std::shared_ptr<QuadratureFunction> sigma,
                           const std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs)
  : assembler(a, L, bcs), eps_p_data(a->mesh(), sigma->element(), 6),
    eps_p_eq_data(a->mesh(), sigma->element(), 1),
    //_tangent(tangent), _sigma(sigma), 
    _u(u)
{
  // FIXME: Test assumptions
  parameters.add("num_preconditioner_rebuilds", -1);
}
//-----------------------------------------------------------------------------
PlasticityProblem::PlasticityProblem(std::shared_ptr<const dolfin::Form> a,
                                     std::shared_ptr<const dolfin::Form> L,
                                     std::shared_ptr<dolfin::Function> u,
                                     std::shared_ptr<QuadratureFunction> tangent,
                                     std::shared_ptr<QuadratureFunction> sigma,
                                     const std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs,
                                     std::shared_ptr<dolfin::NewtonSolver> newton_solver)
  : assembler(a, L, bcs), eps_p_data(a->mesh(), sigma->element(), 6),
    eps_p_eq_data(a->mesh(), sigma->element(), 1),
    //_tangent(tangent), _sigma(sigma),
    _u(u), _newton_solver(newton_solver)
{
  // FIXME: Test assumptions
  parameters.add("num_preconditioner_rebuilds", -1);
}
//-----------------------------------------------------------------------------
PlasticityProblem::~PlasticityProblem()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void PlasticityProblem::F(dolfin::GenericVector& b,
                          const dolfin::GenericVector& x)

{
  // Do nothing
}
//-----------------------------------------------------------------------------
void PlasticityProblem::J(dolfin::GenericMatrix& A,
                          const dolfin::GenericVector& x)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void PlasticityProblem::form(dolfin::GenericMatrix& A,
                             dolfin::GenericMatrix& P,
                             dolfin::GenericVector& b,
                             const dolfin::GenericVector& x)
{
  dolfin::Timer timer("PlasticityProblem form");

  // Update displacement ghost values
  _u->update();

  // Build A and b tensors
  form_tensors(A, b, x);

  // PC same as operator
  //P = &A;
}
//-----------------------------------------------------------------------------
void PlasticityProblem::update_history()
{
  // Update history parameters
  eps_p_data.update_history();
  eps_p_eq_data.update_history();
}
//-----------------------------------------------------------------------------
void PlasticityProblem::form_tensors(dolfin::GenericMatrix& A,
                                     dolfin::GenericVector& b,
                                     const dolfin::GenericVector& x)
{
  int iterate = -1;
  const int num_pc_rebuilds = parameters["num_preconditioner_rebuilds"];
  if (_newton_solver)
  {
    const std::size_t iterate = _newton_solver->iteration();
    if (num_pc_rebuilds >= 0)
    {
      dolfin::GenericLinearSolver& linear_solver
        = _newton_solver->linear_solver();
      if (iterate < num_pc_rebuilds)
      {
        dolfin::log(13, "Rebuild preconditioner for FEniCS Solid Mechanics");
        linear_solver.parameters("preconditioner")["structure"]
          = "same_nonzero_pattern";
      }
      else
      {
        dolfin::log(13, "Reuse preconditioner for FEniCS Solid Mechanics");
        linear_solver.parameters("preconditioner")["structure"] = "same";
      }
    }
  }

  // Assemble
  dolfin::Timer timer("Assemble PlasticityProblem LHS/RHS");
  dolfin::set_log_active(false);
  assembler.assemble(A, b, x);
  dolfin::set_log_active(true);
}
//-----------------------------------------------------------------------------
