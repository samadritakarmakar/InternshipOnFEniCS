// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#ifndef __Drucker_Prager_H
#define __Drucker_Prager_H

#include <Eigen/Dense>
#include "PlasticityModel.h"

namespace fenicssolid
{

  class DruckerPrager : public PlasticityModel
  {

  public:

    /// Delete copy constructor and assignement
    DruckerPrager& operator=(const DruckerPrager&) = delete;  // Disallow copying
    DruckerPrager(const DruckerPrager&) = delete;

    /// Constructor
    DruckerPrager(double E, double nu, double friction_angle,
                  double dilatancy_angle, double cohesion,
                  double hardening_parameter);

    /// Hardening parameter
    double hardening_parameter(double equivalent_plastic_strain) const;

    /// Value of yield function f
    double f(const Eigen::Matrix<double, 6, 1>& stress,
             double equivalent_plastic_strain) const;

    /// First derivative of f with respect to sigma
    void df(Eigen::Matrix<double, 6, 1>& df_dsigma,
            const Eigen::Matrix<double, 6, 1>& stress) const;

    /// First derivative of g with respect to sigma
    void dg(Eigen::Matrix<double, 6, 1>& dg_dsigma,
            const Eigen::Matrix<double, 6, 1>& stress) const ;

    /// Second derivative of g with respect to sigma
    void ddg(Eigen::Matrix<double, 6, 6>& ddg_ddsigma,
             const Eigen::Matrix<double, 6, 1>& stress) const;

  private:

    // Computes effective stresses
    double effective_stress(const Eigen::Matrix<double, 6, 1>& stress) const;

    // Model parameters
    const double _hardening_parameter;
    const double alpha_friction, k_friction, alpha_dilatancy;

    // Auxiliary variables
    Eigen::Matrix<double, 6, 6> A;

  };
}

#endif
