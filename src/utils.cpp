// Copyright (C) 2009-2012 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2009-10-02
// Last changed: 2012-07-17

#include <iostream>

#include <boost/multi_array.hpp>
#include <ufc.h>
#include <ufc_geometry.h>
#include <dolfin/fem/FiniteElement.h>
#include "utils.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
std::string fenicssolid::git_commit_hash()
{
  return FENICSSOLID_GIT_COMMIT_HASH;
}
//-----------------------------------------------------------------------------
