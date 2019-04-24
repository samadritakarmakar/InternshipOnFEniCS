// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-03-03

#ifndef __RETURN_MAPPING_H
#define __RETURN_MAPPING_H

#include <utility>
#include <Eigen/Dense>
#include <dolfin/common/Variable.h>

namespace fenicssolid
{

  class PlasticityModel;

  class ReturnMapping : public dolfin::Variable
  {
  public:

    /// Delete copy constructor and assignement
    ReturnMapping& operator=(const ReturnMapping&) = delete;  // Disallow copying
    ReturnMapping(const ReturnMapping&) = delete;

    /// Constructor
    ReturnMapping(const unsigned int maxit = 50);

    /// Destructor
    ~ReturnMapping();

    /// Closest point projection return mapping
    std::pair<bool, unsigned int>
      closest_point_projection(std::shared_ptr<const PlasticityModel> plastic_model,
                               Eigen::Matrix<double, 6, 6>& consistent_tangent,
                               Eigen::Matrix<double, 6, 1>& trial_stresses,
                               Eigen::Matrix<double, 6, 1>& plastic_strain,
                               double& equivalent_plastic_strain,
                               bool use_plastic_tangent=false);
  private:

    // Maximum number of iterations
    const unsigned int _maxit;
  };
}

#endif
