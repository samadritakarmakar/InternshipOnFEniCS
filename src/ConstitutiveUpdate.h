// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#ifndef __CONSTITUTIVE_UPDATE_H
#define __CONSTITUTIVE_UPDATE_H

#include <vector>
#include <Eigen/Dense>

#include <boost/multi_array.hpp>
#include "HistoryData.h"
#include "ReturnMapping.h"
#include "StateUpdate.h"

namespace dolfin
{
  class Cell;
  class FiniteElement;
  class Function;
  class GenericDofMap;
}

namespace fenicssolid
{

  class PlasticityModel;
  class UFLQuadratureFunction;

  class ConstitutiveUpdate : public StateUpdate
  {
  public:

    /// Delete copy constructor and assignement
    ConstitutiveUpdate& operator=(const ConstitutiveUpdate&) = delete;  // Disallow copying
    ConstitutiveUpdate(const ConstitutiveUpdate&) = delete;

    /// Constructor
    ConstitutiveUpdate(std::shared_ptr<const UFLQuadratureFunction> strain_func,
                       std::shared_ptr<const PlasticityModel> plastic_model);

    /// Destructor
    ~ConstitutiveUpdate();

    /// Update stress for cell
    void update(const dolfin::Cell& cell,
                const double* vertex_coordinates);

    /// Update history variables
    void update_history();

    const std::vector<double>& w_stress() const
    { return _w_stress; }

    const std::vector<double>& w_tangent() const
    { return _w_tangent; }

    std::shared_ptr<const HistoryData> eps_p_eq() const
    { return _eps_p_equiv; }

    std::shared_ptr<const HistoryData> eps_p() const
    { return _eps_p; }

  private:
    // Index of previously computed cell
    std::size_t _prev_cell_index;
    // Previously computed strain
    std::vector<double> _prev_w_strain;

    
    // Restriction to UFC cell
    std::vector<double> _w_strain;
    std::vector<double> _w_stress;
    std::vector<double> _w_tangent;

    // Strain tensor UFLQuadratureFunction
    std::shared_ptr<const UFLQuadratureFunction> _strain_func;

    // Plasticity model
    std::shared_ptr<const PlasticityModel> _plastic_model;

    // Return mapping object
    ReturnMapping return_mapping;

    // Internal variables
    std::shared_ptr<HistoryData> _eps_p;
    std::shared_ptr<HistoryData> _eps_p_equiv;

    // Track points that deformed plastically at end of last increment
    boost::multi_array<bool, 2> _plastic_last;

    std::size_t _num_ip_dofs;
    std::size_t _num_ip_per_cell;

    // Stress quadrature points in reference coordinates
    boost::multi_array<double, 2> _X;

    // Elastic tangent
    const Eigen::Matrix<double, 6, 6>& _De;

  };
}

#endif
