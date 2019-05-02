// Copyright (C) 2006-2011 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2011-02-06

#include <dolfin.h>
#include <FenicsSolidMechanics.h>
#include <math.h>

// Switch between linear and quadratic elements
//#include "../forms/p2_forms/Plas3D.h"
//#include "../forms/p1_forms/Plas3D.h"
#include "Plas3D.h"
using namespace dolfin;

// Displacement right end
class DBval : public Expression
{
  public:

    DBval(const double& t) : t(t), Expression(3) {}

    void eval(Array<double>& values, const Array<double>& x) const
    {
      // Stretch
      values[0] = 0;
      values[1] = 0;
      values[2] = -t;

    }
    const double& t;
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundaryX1 : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      //return std::abs(x[2] - 0.05) < DOLFIN_EPS;
      return x[2] >= 0.05 - DOLFIN_EPS;

  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundaryX0 : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      return (x[2] < DOLFIN_EPS );
  }
};

class DirichletBoundaryX2 : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      return on_boundary && (std::sqrt(x[0]*x[0]+x[1]*x[1]) > 0.0495-DOLFIN_EPS);
      //return on_boundary &&(x[2]>DOLFIN_EPS && x[2]<0.05-DOLFIN_EPS);
  }
};
int main()
{
  Timer timer("Total plasicity solver time");

  //dolfin::parameters["reorder_dofs_serial"] = false;
  //dolfin::parameters["dof_ordering_library"] = "SCOTCH";
  //dolfin::parameters["linear_algebra_backend"] = "Epetra";

  // Create mesh
  //auto mesh = std::make_shared<UnitCubeMesh>(8, 8, 8);
  auto mesh=std::make_shared<Mesh>("cylinder.xml");

  // Young's modulus and Poisson's ratio
  double E = 200.0e9;
  double nu = 0.3;

  // Time parameter
  double t = 0.0;

  // Elastic time step, always one step.
  double Edt  = 0.00001;

  // Load region 0, time step and number of steps
  double dt0 = 0.00001;
  unsigned int dt0_steps = 3;

  // Load region 1, time step and number of steps
  double dt1 = -0.002;
  unsigned int dt1_steps =  1;

  // Load region 2, time step and number of steps
  double dt2 = 0.001;
  unsigned int dt2_steps =  4;

  // Source term, RHS
  auto f = std::make_shared<Constant>(0.0, 0.0, 0.0);

  // Function spaces
  auto V = std::make_shared<Plas3D::CoefficientSpace_f>(mesh);
  dolfin::cout << "Number of dofs: " << V->dim() << dolfin::endl;

  // Extract elements for stress and tangent
  std::shared_ptr<const FiniteElement> element_t;
  {
    //Plas3D::BilinearForm::CoefficientSpace_t Vt(mesh);//ADDED by SAM //Original. Not working!!!
    Plas3D::CoefficientSpace_t Vt(mesh);
    element_t = Vt.element();
  }

  auto Vs = std::make_shared<Plas3D::Form_L::CoefficientSpace_s>(mesh);
  auto element_s = Vs->element();

  // Create boundary conditions (use SubSpace to apply simply
  // supported BCs)
  auto zero = std::make_shared<Constant>(0.0, 0.0, 0.0);
  auto val = std::make_shared<DBval>(t);
  auto zero0 = std::make_shared<Constant>(0.0);
  auto Vx= V->sub(0);
  auto Vy= V->sub(1);
  auto Vz= V->sub(2);

  auto dbX0 = std::make_shared<DirichletBoundaryX0>();
  auto dbX1 = std::make_shared<DirichletBoundaryX1>();
  auto dbX2 = std::make_shared<DirichletBoundaryX2>();
  auto bcX0 = std::make_shared<DirichletBC>(V, zero, dbX0);
  auto bcX1 = std::make_shared<DirichletBC>(V, val, dbX1);
  auto bcX2 = std::make_shared<DirichletBC>(Vx, zero0, dbX2);
  auto bcX3 = std::make_shared<DirichletBC>(Vy, zero0, dbX2);


  std::vector<std::shared_ptr<const DirichletBC>> bcs = {bcX0, bcX1, bcX2, bcX3};

  // Slope of hardening (linear) and hardening parameter
  const double E_t = 0.1*E;
  const double hardening_parameter = E_t/(1.0 - E_t/E);

  // Yield stress
  const double yield_stress = 235.0e6;

  // Solution function
  auto u = std::make_shared<Function>(V);

  // Object of class von Mises
  auto J2 = std::make_shared<const fenicssolid::VonMises>(E, nu, yield_stress,
                                                          hardening_parameter);

  auto eps_form = std::make_shared<Plas3D::Form_eps_form>(Vs);
  eps_form->f = u;
  auto mesh_ref = std::make_shared<Mesh>(UnitTetrahedronMesh::create());
  //auto mesh_ref = std::make_shared<Mesh>("cylinder.xml");
  auto Vs_ref = std::make_shared<Plas3D::CoefficientSpace_s>(mesh_ref);
  auto eps_ref_form = std::make_shared<Plas3D::Form_eps_ref_form>(Vs_ref);

  // Strain UFLQuadratureFunction
  const auto Qdef = std::make_shared<fenicssolid::UFLQuadratureFunction>(eps_form, eps_ref_form);

  // Constituive update
  auto constitutive_update
    = std::make_shared<fenicssolid::ConstitutiveUpdate>(Qdef, J2);

  // Create forms and attach functions
  auto tangent
    = std::make_shared<fenicssolid::QuadratureFunction>(mesh, element_t,
                                                        constitutive_update,
                                                        constitutive_update->w_tangent());

  auto a = std::make_shared<Plas3D::Form_a>(V, V);
  a->t = tangent;

  auto L = std::make_shared<Plas3D::Form_L>(V);
  L->f = f;
  auto stress = std::make_shared<fenicssolid::QuadratureFunction>(mesh, element_s,
                                         constitutive_update->w_stress());
  L->s = stress;

  //ADDED by SAM
  //auto V2 = std::make_shared<Plas3D::CoefficientSpace_s>(mesh);
  auto V2 =std::make_shared<Plas3D::CoefficientSpace_fstrss>(mesh); //ADDED BY Q
  //auto V3 =std::make_shared<Plas3D::CoefficientSpace_s2>(mesh);
  auto strss =std::make_shared<Function>(V2);
  auto aStrss =std::make_shared<Plas3D::Form_aStrss>(V2,V2);
  auto LStrss = std::make_shared<Plas3D::Form_L_Strss>(V2);
  LStrss->s2 =stress;
  //auto f2 = std::make_shared<Constant>(0.0, 0.0, 0.0, );
  std::vector<std::size_t> value_shape;
  value_shape.push_back(3);
  value_shape.push_back(3);
  std::vector<double> values;
  for (int i=0; i<9; i++)
      values.push_back(0.0);

  auto f2=std::make_shared<Constant>(value_shape, values); //ADDED BY Q
  LStrss->fstrss =f2 ; //ADDED BY Q
//
//ADDED BY SAM For Eps
  auto V3= std::make_shared<Plas3D::CoefficientSpace_fEps>(mesh);
  auto Eps= std::make_shared<Function>(V3);
  auto aEps = std::make_shared<Plas3D::Form_aEps>(V3,V3);
  auto LEps = std:: make_shared<Plas3D::Form_LEps>(V3);
  LEps->u2=u;
  auto fEps= std::make_shared<Constant>(value_shape, values);
  LEps->fEps =fEps;
//

//ADDED BY SAM For Eps_p
  auto V4= std::make_shared<Plas3D::CoefficientSpace_fEps_p>(mesh);
  auto Eps_p= std::make_shared<Function>(V4);
  auto aEps_p = std::make_shared<Plas3D::Form_aEps_p>(V4,V4);
  auto LEps_p =std::make_shared<Plas3D::Form_LEps_p>(V4);
  LEps_p->eps_p=constitutive_update->eps_p();
  auto fEps_p= std::make_shared<Constant>(value_shape, values);
  LEps_p->fEps_p =fEps_p;
  //
//ADDED BY SAM For stress2
  auto V5= std::make_shared<Plas3D::CoefficientSpace_fstrss>(mesh);
  auto strss2 =std::make_shared<Function>(V5);
  auto aStrss2 =std::make_shared<Plas3D::Form_aStrss2>(V5, V5);
  auto LStress2 = std::make_shared<Plas3D::Form_L_Strss2>(V5);
  LStress2->u2 = u;
  LStress2->eps_p=constitutive_update->eps_p();
  LStress2->fstrss=f2;
  //

  // Create PlasticityProblem
  auto nonlinear_problem
    = std::make_shared<fenicssolid::PlasticityProblem>(a, L, u, tangent,
                                                       stress, bcs);

  // Create nonlinear solver and set parameters
  dolfin::NewtonSolver nonlinear_solver;
  nonlinear_solver.parameters["convergence_criterion"] = "incremental";
  nonlinear_solver.parameters["maximum_iterations"]    = 50;
  nonlinear_solver.parameters["relative_tolerance"]    = 1.0e-6;
  nonlinear_solver.parameters["absolute_tolerance"]    = 1.0e-15;

  // File names for output
  File file1("output/disp.pvd");
  File file2("output/eq_plas_strain.pvd");
  File file3("output/stress.pvd");
  //File file4("output/eps.pvd");
  XDMFFile file4("output/eps.xdmf");
  XDMFFile file5("output/eps_p.xdmf");
  XDMFFile file6("output/stress2.xdmf");

  // Equivalent plastic strain for visualisation
  auto eps_eq = std::make_shared<MeshFunction<double>>(mesh, mesh->topology().dim());
  //auto eps = std::make_shared<MeshFunction<double>>(mesh, mesh->topology().dim());

  // Load-disp info
  unsigned int step = 0;
  //unsigned int steps = dt0_steps + dt1_steps + dt2_steps + 1;
  unsigned int steps = 10;
  while (step < steps)
  {
    // Use elastic tangent for first time step
    if (step == 0)
      t += Edt;
    else
        t += dt0;

    step++;
    std::cout << "step begin: " << step << std::endl;
    std::cout << "time: " << t << std::endl;

    // Solve non-linear problem
    nonlinear_solver.solve(*nonlinear_problem, *u->vector());
    // Update variables
    constitutive_update->update_history();
    //ADDED BY SAM
    dolfin::solve(*aStrss==*LStrss, *strss);
    dolfin::solve(*aEps==*LEps, *Eps);
    dolfin::solve(*aEps_p==*LEps_p, *Eps_p);
    dolfin::solve(*aStrss2==*LStress2, *strss2);
    // Write output to files
    file1 << *u;
    constitutive_update->eps_p_eq()->compute_mean(eps_eq);
    file2 << *eps_eq;
    file3 << *strss;
    file4.write(*Eps,t);
    file5.write(*Eps_p,t);
    file6.write(*strss2,t);
    //constitutive_update->eps_p();
  /*  for (int i=0; i<constitutive_update->w_stress().size(); i++)
    {
        std::cout<<constitutive_update->w_stress()[i]<<"\n";
    }
    std::cout<<"\n";*/
  }
  //file3 << *mesh;
  cout << "Solution norm: " << u->vector()->norm("l2") << endl;

  timer.stop();
  dolfin::list_timings(dolfin::TimingClear::clear, {dolfin::TimingType::wall});

  return 0;
}
