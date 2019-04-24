// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#include "DruckerPrager.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
DruckerPrager::DruckerPrager(double E, double nu, double friction_angle,
                             double dilatancy_angle, double cohesion,
                             double hardening_parameter)
    : PlasticityModel(E, nu), _hardening_parameter(hardening_parameter),
      alpha_friction(6.0*sin(friction_angle)/(3.0 - sin(friction_angle))),
      k_friction(6.0*cohesion*cos(friction_angle)/(3.0 - sin(friction_angle))),
      alpha_dilatancy(6.0*sin(dilatancy_angle)/(3.0 - sin(dilatancy_angle))),
      A(Eigen::Matrix<double, 6, 6>::Constant(0.0))
{
  A(0, 0) =  2; A(1, 1) =  2; A(2, 2) =  2;
  A(3, 3) =  6; A(4, 4) =  6; A(5, 5) =  6;
  A(0, 1) = -1; A(0, 2) = -1; A(1, 0) = -1;
  A(1, 2) = -1; A(2, 0) = -1; A(2, 1) = -1;
}
//-----------------------------------------------------------------------------
double
DruckerPrager::hardening_parameter(double equivalent_plastic_strain) const
{
  return _hardening_parameter;
}
//-----------------------------------------------------------------------------
double DruckerPrager::f(const Eigen::Matrix<double, 6, 1>& stress,
                        double equivalent_plastic_strain) const
{

  return effective_stress(stress) + alpha_friction/3.0*(stress[0] + stress[1] + stress[2])
         - k_friction - _hardening_parameter*equivalent_plastic_strain;
}
//-----------------------------------------------------------------------------
void DruckerPrager::df(Eigen::Matrix<double, 6, 1>& df_dsigma,
                       const Eigen::Matrix<double, 6, 1>& stress) const
{
  const double _effective_stress = effective_stress(stress);
  df_dsigma[0]  = (2.0*stress[0] - stress[1]
                   - stress[2])/(2.0*_effective_stress) + alpha_friction/3.0;
  df_dsigma[1] = (2.0*stress[1] - stress[0]
                  - stress[2])/(2.0*_effective_stress) + alpha_friction/3.0;
  df_dsigma[2] = (2.0*stress[2] - stress[0]
                  - stress[1])/(2.0*_effective_stress) + alpha_friction/3.0;
  df_dsigma[3] = 6.0*stress[3]/(2.0*_effective_stress);
  df_dsigma[4] = 6.0*stress[4]/(2.0*_effective_stress);
  df_dsigma[5] = 6.0*stress[5]/(2.0*_effective_stress);
}
//-----------------------------------------------------------------------------
void DruckerPrager::dg(Eigen::Matrix<double, 6, 1>& dg_dsigma,
                       const Eigen::Matrix<double, 6, 1>& stress) const
{
  const double _effective_stress = effective_stress(stress);
  dg_dsigma[0] = (2.0*stress[0] - stress[1]
                  - stress[2])/(2.0* _effective_stress) + alpha_friction/3.0;
  dg_dsigma[1] = (2.0*stress[1] - stress[0]
                  - stress[2])/(2.0* _effective_stress) + alpha_friction/3.0;
  dg_dsigma[2] = (2.0*stress[2] - stress[0]
                  - stress[1])/(2.0* _effective_stress) + alpha_friction/3.0;
  dg_dsigma[3] = 6.0*stress[3]/(2.0*_effective_stress);
  dg_dsigma[4] = 6.0*stress[4]/(2.0*_effective_stress);
  dg_dsigma[5] = 6.0*stress[5]/(2.0*_effective_stress);
}
//-----------------------------------------------------------------------------
void DruckerPrager::ddg(Eigen::Matrix<double, 6, 6>& ddg_ddsigma,
                        const Eigen::Matrix<double, 6, 1>& stress) const
{
  const double _effective_stress = effective_stress(stress);

  // First derivative of g with respect to sigma, excluding (+
  // alpha_dilatancy/3.0) on diagonal terms
  Eigen::Matrix<double, 6, 1> dg_dsigma_mod;
  dg_dsigma_mod[0] = (2.0*stress[0] - stress[1]
                      - stress[2])/(2.0*_effective_stress);
  dg_dsigma_mod[1] = (2.0*stress[1] - stress[0]
                      - stress[2])/(2.0*_effective_stress);
  dg_dsigma_mod[2] = (2.0*stress[2] - stress[0]
                      - stress[1])/(2.0*_effective_stress);
  dg_dsigma_mod[3] = 6.0*stress[3]/(2.0*_effective_stress);
  dg_dsigma_mod[4] = 6.0*stress[4]/(2.0*_effective_stress);
  dg_dsigma_mod[5] = 6.0*stress[5]/(2.0*_effective_stress);

  ddg_ddsigma = A/(2.0*_effective_stress) -
    dg_dsigma_mod*dg_dsigma_mod.transpose()/_effective_stress;
}
//-----------------------------------------------------------------------------
double
DruckerPrager::effective_stress(const Eigen::Matrix<double, 6, 1>& stress) const
{
  return sqrt(std::pow((stress[0] + stress[1] + stress[2]), 2)
           - 3.0*(stress[0] * stress[1] + stress[0]*stress[2]
           + stress[1]*stress[2] - stress[3]*stress[3]
          - stress[4]*stress[4] - stress[5]*stress[5]));
}
//-----------------------------------------------------------------------------
