// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#include <iostream>

#include <dolfin/common/constants.h>
#include "PlasticityModel.h"
#include "ReturnMapping.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
ReturnMapping::ReturnMapping(const unsigned int maxit) : _maxit(maxit)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ReturnMapping::~ReturnMapping()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::pair<bool, unsigned int>
ReturnMapping::closest_point_projection(std::shared_ptr<const PlasticityModel> plastic_model,
                                        Eigen::Matrix<double, 6, 6>& D,
                                        Eigen::Matrix<double, 6, 1>& trial_stresses,
                                        Eigen::Matrix<double, 6, 1>& plastic_strain,
                                        double& equivalent_plastic_strain,
                                        bool use_plastic_tangent)
{
  // Initialise variables
  double delta_lambda = 0.0;
  bool plastic_flag = false;
  unsigned int num_iterations = 0;

  // Work arrays (statically allocated)
  Eigen::Matrix<double, 6, 6> Q;
  Eigen::Matrix<double, 6, 6> Qinv;
  Eigen::Matrix<double, 6, 6>  R, ddg_ddsigma, inverse_Q;
  Eigen::Matrix<double, 6, 1> Rm, Rn, RinvQ;

  // Variables for return mapping
  Eigen::Matrix<double, 6, 1> df_dsigma, dg_dsigma,
    sigma_current,
    sigma_dot, sigma_residual;


  // Trial stress
  sigma_current = trial_stresses;


  // Auxiliary variables to speed up return mapping

  // Elastic tangent
  const Eigen::Matrix<double, 6, 6>& De = plastic_model->elastic_tangent();

  // Evaluate hardening parameter
  double hardening_parameter
    = plastic_model->hardening_parameter(equivalent_plastic_strain);

  // Evaluate yield function (trial stresses)
  double residual_f = plastic_model->f(sigma_current, equivalent_plastic_strain);

  // Check for yielding
  if (residual_f/sigma_current.norm() > 1.0e-12)
  {
    // Compute normal vectors to yield surface and plastic potential
    plastic_model->df(df_dsigma, sigma_current);
    plastic_model->dg(dg_dsigma, sigma_current);
    plastic_flag = true;

    // Perform Newton iterations to project stress onto yield surface
    while (std::abs(residual_f)/sigma_current.norm() > 1e-12)
    {
      num_iterations++;
      if (num_iterations > _maxit)
        dolfin::error("Return mapping iterations > %d.", _maxit);

      // Reset sigma_residual in first step
      if (num_iterations == 1)
        sigma_residual.setZero();

      // Compute second derivative (with respect to stresses) of
      // plastic potential
      plastic_model->ddg(ddg_ddsigma, sigma_current);

      // Compute auxiliary matrix Q
      Q = Eigen::Matrix<double, 6, 6>::Identity() + delta_lambda*De*ddg_ddsigma;

      // Invert Q
      Qinv = Q.inverse();

      // Compute auxiliary matrix R
      R = Qinv*De;

      // lambda_dot, rate of plastic multiplier
      const double residual_tmp = residual_f - sigma_residual.dot(Qinv*df_dsigma);
      double lambda_dot = residual_tmp/(df_dsigma.dot(R*dg_dsigma)
                                        + hardening_parameter);

      // Compute stress increment
      // FIXME: (GNW) is the below correction? In uBLAS it was prod(x, A)
      sigma_dot = Qinv.transpose()*(-lambda_dot*De*dg_dsigma -sigma_residual);

      // Increment plastic multiplier
      delta_lambda += lambda_dot;

      // Update current stress state
      sigma_current += sigma_dot;

      // Update equivalent plastic strain
      equivalent_plastic_strain
        = plastic_model->kappa(equivalent_plastic_strain, sigma_current,
                              lambda_dot);

      // Compute hardening parameter
      hardening_parameter
        = plastic_model->hardening_parameter(equivalent_plastic_strain);

      // Evaluate yield function at new stress state
      residual_f = plastic_model->f(sigma_current, equivalent_plastic_strain);

      // Compute normal to yield surface at new stress state
      plastic_model->df(df_dsigma, sigma_current);

      // Compute normal normal to plastic potential at new stress state
      plastic_model->dg(dg_dsigma, sigma_current);

      // Compute residual vector
      sigma_residual = sigma_current
        - (trial_stresses - delta_lambda*De*dg_dsigma);
    }

    // Update matrices
    plastic_model->ddg(ddg_ddsigma, sigma_current);

    // Compute matrix Q
    Q = Eigen::Matrix<double, 6, 6>::Identity() + delta_lambda*De*ddg_ddsigma;

    // Invert Q
    Qinv = Q.inverse();

    // Compute matrix R and vector Rn
    R = Qinv*De;
    Rn = R*df_dsigma;

    // Compute consistent tangent operator
    D = R - R*(dg_dsigma*Rn.transpose())/(df_dsigma.dot(R*dg_dsigma)
                                          + hardening_parameter);

    // Stresses for next Newton iteration, trial stresses are
    // overwritten by current stresses
    trial_stresses = sigma_current;

    // Update plastic strains
    plastic_strain += delta_lambda*dg_dsigma;
  }
  else if (use_plastic_tangent)
  {
    // Compute continuum tangent operator
    plastic_model->df(df_dsigma, sigma_current);
    plastic_model->dg(dg_dsigma, sigma_current);
    Rn = De*df_dsigma;
    const double denom = df_dsigma.dot(De*dg_dsigma) + hardening_parameter;
    D = De - De*(dg_dsigma*(De*df_dsigma).transpose())/denom;
  }
  else
  {
    // Use elastic tangent
    D = De;
  }

  return std::make_pair(plastic_flag, num_iterations);
}
//-----------------------------------------------------------------------------
