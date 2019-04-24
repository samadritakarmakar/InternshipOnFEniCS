// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#include "VonMises.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
VonMises::VonMises(double E, double nu, double yield_stress,
                   double hardening_parameter)
  : PlasticityModel(E, nu), yield_stress(yield_stress),
    _hardening_parameter(hardening_parameter)
{
  A.setZero();
  A(0, 0) =  2; A(1, 1) =  2; A(2, 2) =  2;
  A(3, 3) =  6; A(4, 4) =  6; A(5, 5) =  6;
  A(0, 1) = -1; A(0, 2) = -1; A(1, 0) = -1;
  A(1, 2) = -1; A(2, 0) = -1; A(2, 1) = -1;
}
//-----------------------------------------------------------------------------
double VonMises::hardening_parameter(double eps_p_eq) const
{
  return _hardening_parameter;
}
//-----------------------------------------------------------------------------
double VonMises::f(const Eigen::Matrix<double, 6, 1>& stress,
                   double eps_p_eq) const
{
  return effective_stress(stress) - yield_stress
    - _hardening_parameter*eps_p_eq;
}
//-----------------------------------------------------------------------------
void VonMises::df(Eigen::Matrix<double, 6, 1>& df_dsigma,
                  const Eigen::Matrix<double, 6, 1>& stress) const
{
  const double _effective_stress = effective_stress(stress);

  df_dsigma[0] = (2*stress[0] - stress[1] - stress[2]) / (2*_effective_stress);
  df_dsigma[1] = (2*stress[1] - stress[0] - stress[2]) / (2*_effective_stress);
  df_dsigma[2] = (2*stress[2] - stress[0] - stress[1]) / (2*_effective_stress);
  df_dsigma[3] = 6*stress[3] / (2*_effective_stress);
  df_dsigma[4] = 6*stress[4] / (2*_effective_stress);
  df_dsigma[5] = 6*stress[5] / (2*_effective_stress);
}
//-----------------------------------------------------------------------------
void VonMises::ddg(Eigen::Matrix<double, 6, 6>& ddg_ddsigma,
                   const Eigen::Matrix<double, 6, 1>& stress) const
{
  Eigen::Matrix<double, 6, 1> dg_dsigma;
  df(dg_dsigma, stress);
  ddg_ddsigma = A/(2*effective_stress(stress))
    - dg_dsigma*dg_dsigma.transpose()/effective_stress(stress);
}
//-----------------------------------------------------------------------------
double
VonMises::effective_stress(const Eigen::Matrix<double, 6, 1>& stress) const
{
  const double trace = stress[0] + stress[1] + stress[2];
  return std::sqrt(trace*trace - 3*(stress[0]*stress[1] + stress[0]*stress[2]
                                 + stress[1]*stress[2] - stress[3]*stress[3]
                                 - stress[4]*stress[4] - stress[5]*stress[5]));
}
//-----------------------------------------------------------------------------
