// Copyright (C) 2009-2017 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.

#ifndef __FENICS_SOLID_UTILS_H
#define __FENICS_SOLID_UTILS_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>

namespace fenicssolid
{

  /// Return git commit hash for library
  std::string git_commit_hash();

};

#endif
