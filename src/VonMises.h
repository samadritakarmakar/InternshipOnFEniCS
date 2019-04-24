// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#ifndef __VON_MISES_H
#define __VON_MISES_H

#include <Eigen/Dense>
#include "PlasticityModel.h"

namespace fenicssolid
{

  class VonMises : public PlasticityModel
  {
  public:

    /// Delete copy constructor and assignement
    VonMises& operator=(const VonMises&) = delete;  // Disallow copying
    VonMises(const VonMises&) = delete;

    /// Constructor
    VonMises(double E, double nu, double yield_stress,
             double hardening_parameter);

    /// Hardening parameter
    double hardening_parameter(double eps_p_eq) const;

    /// Value of yield function f
    inline double f(const Eigen::Matrix<double, 6, 1>& stress,
                    double eps_p_eq) const;

    /// First derivative of f with respect to sigma
    inline void df(Eigen::Matrix<double, 6, 1>& df_dsigma,
                   const Eigen::Matrix<double, 6, 1>& stress) const;

    /// Second derivative of g with respect to sigma
    inline void ddg(Eigen::Matrix<double, 6, 6>& ddg_ddsigma,
                    const Eigen::Matrix<double, 6, 1>& stress) const;

  private:

    // Computes effective stresses
    inline double
      effective_stress(const Eigen::Matrix<double, 6, 1>& stress) const;

    // Model parameters
    const double yield_stress, _hardening_parameter;

    // Auxiliary variables
    Eigen::Matrix<double, 6, 6> A;
  };
}
#endif
