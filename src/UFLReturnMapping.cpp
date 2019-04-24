// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#include <iostream>

#include <dolfin/common/constants.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/fem/assemble_local.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/Form.h>

#include "UFLReturnMapping.h"
#include "QuadratureData.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
UFLReturnMapping::UFLReturnMapping(std::shared_ptr<const dolfin::Form> return_equation,
	    	std::shared_ptr<const dolfin::Form> return_equation_J,
		const unsigned int maxit,
		const double tol) : 
		_ufl_return_mapping(return_equation),
		_ufl_return_mapping_Jacobian(return_equation_J),
		_maxit(maxit),
		_tol(tol)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
UFLReturnMapping::~UFLReturnMapping()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::pair<bool, unsigned int>
UFLReturnMapping::closest_point_projection(const QuadratureData* x, const dolfin::Cell& cell) const
{
  bool plastic_flag = false;
  unsigned int num_iterations = 0;
  
  assemble_local(_A_e, *_ufl_return_mapping, cell);

  _sol.resize(_A_e.rows(), 1);
  //if (_A_e.norm() > _tol) std::cerr << _A_e  << "\n------------------\n" ;
  while (_A_e.norm() > _tol) {
    plastic_flag = true;
    num_iterations++;
    if (num_iterations > _maxit)
      dolfin::error("Return mapping iterations > %d.", _maxit);
    // Get bilinear form dofmaps
    std::array<std::shared_ptr<const dolfin::GenericDofMap>, 2> dofmaps_J
      = {{_ufl_return_mapping_Jacobian->function_space(0)->dofmap(), _ufl_return_mapping_Jacobian->function_space(1)->dofmap()}};
    dolfin_assert(dofmaps_J[0] and dofmaps_J[1]);

    //auto dofs_J1 = dofmaps_J[1]->cell_dofs(cell.index());

    assemble_local(_A_J, *_ufl_return_mapping_Jacobian, cell);
    _lu.compute(_A_J);
    _sol = _lu.solve(_A_e);
    const_cast<QuadratureData*>(x)->add_to_cell_values(_sol.data(), cell);
    assemble_local(_A_e, *_ufl_return_mapping, cell);
    //if (_A_e.norm() > _tol) std::cerr << _A_e  << "\n................\n" ;
}
  //if (num_iterations > 0) std::cerr << _A_e << "\nConverged in " << num_iterations << " " << _A_e.norm() << std::endl;

  return std::make_pair(plastic_flag, num_iterations);
}

//-----------------------------------------------------------------------------
