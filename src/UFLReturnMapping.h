// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-03-03

#ifndef __UFL_RETURN_MAPPING_H
#define __UFL_RETURN_MAPPING_H

#include <utility>
#include <Eigen/Dense>
#include <dolfin/common/Variable.h>
//#include <dolfin/fem/Form.h>

namespace dolfin {
  class Form;
  class Cell;
}

namespace fenicssolid
{

  class QuadratureData;
  
  class UFLReturnMapping
  {
  public:

    /// Delete copy constructor and assignement
    UFLReturnMapping& operator=(UFLReturnMapping&) = delete;  // Disallow copying
    UFLReturnMapping(const UFLReturnMapping&) = delete;

    /// Constructor
    UFLReturnMapping(std::shared_ptr<const dolfin::Form> return_equation,
	    	std::shared_ptr<const dolfin::Form> return_equation_J,
		const unsigned int maxit = 50,
		const double tol = 1.E-12);

    /// Destructor
    ~UFLReturnMapping();

    /// Closest point projection return mapping
    std::pair<bool, unsigned int>
      closest_point_projection(const QuadratureData* x, const dolfin::Cell& cell) const;
  private:
    // Return Mapping Forms
    std::shared_ptr<const dolfin::Form> _ufl_return_mapping;
    std::shared_ptr<const dolfin::Form> _ufl_return_mapping_Jacobian;

    // Maximum number of iterations
    const unsigned int _maxit;
    
    // Tolerance
    const double _tol;

    // Equations
    mutable Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _A_e;

    // Jacobian matrix
    mutable Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _A_J;

    // Solution
    mutable Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _sol;

    // Solver
    mutable Eigen::PartialPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                      Eigen::RowMajor>> _lu;

  };
}

#endif
