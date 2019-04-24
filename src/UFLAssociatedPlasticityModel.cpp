// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#include "PlasticityModel.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
PlasticityModel::PlasticityModel(std::shared_ptr<const UFLQuadratureFunction> f,
std::shared_ptr<const UFLQuadratureFunction> de_f)
  : _f(f), _de_f(de_f)
{}

//-----------------------------------------------------------------------------
PlasticityModel::~PlasticityModel()
{
  // Do nothing
}
