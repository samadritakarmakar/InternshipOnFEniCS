// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#include "PlasticityModel.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
PlasticityModel::PlasticityModel(double E, double nu)
  : _hardening_parameter(0.0)
{
  // Lame coefficients
  const double lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));
  const double mu = E/(2.0*(1.0 + nu));

  // Create elastic tangent
  _elastic_tangent = Eigen::Matrix<double, 6, 6>::Constant(0.0);
  _elastic_tangent(0, 0) = lambda + 2.0*mu;
  _elastic_tangent(1, 1) = lambda + 2.0*mu;
  _elastic_tangent(2, 2) = lambda + 2.0*mu;
  _elastic_tangent(3, 3) = mu;
  _elastic_tangent(4, 4) = mu;
  _elastic_tangent(5, 5) = mu;
  _elastic_tangent(0, 1) = lambda;
  _elastic_tangent(0, 2) = lambda;
  _elastic_tangent(1, 0) = lambda;
  _elastic_tangent(1, 2) = lambda;
  _elastic_tangent(2, 0) = lambda;
  _elastic_tangent(2, 1) = lambda;
}
//-----------------------------------------------------------------------------
PlasticityModel::~PlasticityModel()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
double PlasticityModel::hardening_parameter(double const eps_p_eq) const
{
  return _hardening_parameter;
}
//-----------------------------------------------------------------------------
double PlasticityModel::kappa(double eps_p_eq,
                              const Eigen::Matrix<double, 6, 1>& stress,
                              double lambda_dot) const
{
  return eps_p_eq += lambda_dot;
}
//-----------------------------------------------------------------------------
void PlasticityModel::dg(Eigen::Matrix<double, 6, 1>& dg_dsigma,
                         const Eigen::Matrix<double, 6, 1>& stress) const
{
  // Assume associative flow (dg/dsigma = df/dsigma)
  df(dg_dsigma, stress);
}
//-----------------------------------------------------------------------------
