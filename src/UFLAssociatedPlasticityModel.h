// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#ifndef __PLASTICITY_MODEL_H
#define __PLASTICITY_MODEL_H

#include <Eigen/Dense>

namespace fenicssolid
{

  class PlasticityModel
  {
  public:

    /// Delete copy constructor and assignement
    PlasticityModel& operator=(const PlasticityModel&) = delete;  // Disallow copying
    PlasticityModel(const PlasticityModel&) = delete;

    /// Constructor
    PlasticityModel(double E, double nu);

    /// Destructor
    virtual ~PlasticityModel();

    /// Hardening parameter
    virtual double hardening_parameter(double eps_p_eq) const;

    /// Equivalent plastic strain
    virtual double kappa(double eps_p_eq,
                         const Eigen::Matrix<double, 6, 1>& stress,
                         double lambda_dot) const;

    /// Value of yield function f
    virtual double f(const Eigen::Matrix<double, 6, 1>& stress,
                     double eps_eq_p) const = 0;

    /// First derivative of f with respect to sigma
    virtual void df(Eigen::Matrix<double, 6, 1>& df_dsigma,
                    const Eigen::Matrix<double, 6, 1>& stress) const = 0;

    /// First derivative of g with respect to sigma
    virtual void dg(Eigen::Matrix<double, 6, 1>& dg_dsigma,
                    const Eigen::Matrix<double, 6, 1>& stress) const;

    /// Second derivative of g with respect to sigma
    virtual void ddg(Eigen::Matrix<double, 6, 6>& ddg_ddsigma,
                     const Eigen::Matrix<double, 6, 1>& stress) const = 0;

    friend class ConstitutiveUpdate;
    friend class PlasticityProblem;
    friend class ReturnMapping;

  private:

    // Model parameters
    double _hardening_parameter;

    // Elastic tangent
    Eigen::Matrix<double, 6, 6> elastic_tangent;
  };
}

#endif
