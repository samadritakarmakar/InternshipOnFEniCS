// Copyright (C) 2009-2012 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2009-10-02
// Last changed: 2012-07-17

#ifndef __STATE_UPDATE_SOLID_H
#define __STATE_UPDATE_SOLID_H

namespace dolfin
{
  class Cell;
}

namespace fenicssolid
{

  class StateUpdate
  {
  public:

    virtual void update(const dolfin::Cell& cell,
                        const double* vertex_coordinates) = 0;

  };

}

#endif
