/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

 This file is part of ASPECT.

 ASPECT is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 ASPECT is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/base_particle.h>
#include <list>

namespace aspect
{
  namespace Particle
  {
    template <int dim>
    inline
    BaseParticle<dim>::BaseParticle (const Point<dim> &new_loc,
                             const double &new_id)
      :
      location (new_loc),
      id (new_id),
      is_local (true),
      check_vel (true),
      val (0)
    {
    }

    template <int dim>
    inline
    BaseParticle<dim>::BaseParticle ()
      :
      location (),
      velocity (),
      id (0),
      is_local (true),
      check_vel (true),
      val(0)
    {
    }


    template <int dim>
    inline
    BaseParticle<dim>::~BaseParticle ()
    {
    }


    template <int dim>
    void
    BaseParticle<dim>::set_location (const Point<dim> &new_loc)
    {
      location = new_loc;
    }

    template <int dim>
    void
    BaseParticle<dim>::set_data_len (const unsigned int data_len)
    {
      val.resize(data_len - dim - dim - 1);
    }

    template <int dim>
    void
    BaseParticle<dim>::write_data (std::vector<double> &data) const
    {
      // Write location data
      for (unsigned int i = 0; i < dim; ++i)
        {
          data.push_back(location(i));
        }
      // Write velocity data
      for (unsigned int i = 0; i < dim; ++i)
        {
          data.push_back(velocity(i));
        }

      data.push_back(id);

      for (unsigned int i = 0; i < val.size(); i++)
        {
          data.push_back(val [i]);
        }
    }

    template <int dim>
    unsigned int
    BaseParticle<dim>::read_data(const std::vector<double> &data, const unsigned int &pos)
    {
      unsigned int p = pos;
      // Read location data
      for (unsigned int i = 0; i < dim; ++i)
        {
          location (i) = data[p++];
        }
      // Write velocity data
      for (unsigned int i = 0; i < dim; ++i)
        {
          velocity (i) = data[p++];
        }
      id = data[p++];

      for (unsigned int i = 0; i < val.size(); i++)
        {
          val [i] = data[p++];
        }

      return p;
    }

    template <int dim>
    const unsigned int
    BaseParticle<dim>::data_len () const
    {
      const unsigned int length = dim + dim + 1 + val.size();

      return length;
    }


    template <int dim>
    Point<dim>
    BaseParticle<dim>::get_location () const
    {
      return location;
    }

    template <int dim>
    void
    BaseParticle<dim>::set_velocity (Point<dim> new_vel)
    {
      velocity = new_vel;
    }

    template <int dim>
    Point<dim>
    BaseParticle<dim>::get_velocity () const
    {
      return velocity;
    }

    template <int dim>
    double
    BaseParticle<dim>::get_id () const
    {
      return id;
    }

    template <int dim>
    void
    BaseParticle<dim>::set_properties (const std::vector<double> new_properties)
    {
      val = new_properties;
    }

    template <int dim>
    const
    std::vector<double>
    BaseParticle<dim>::get_properties () const
    {
      return val;
    }

    template <int dim>
    std::vector<double>&
    BaseParticle<dim>::get_properties ()
    {
      return val;
    }

    template <int dim>
    bool
    BaseParticle<dim>::local () const
    {
      return is_local;
    }

    template <int dim>
    void
    BaseParticle<dim>::set_local (bool new_local)
    {
      is_local = new_local;
    }

    template <int dim>
    bool
    BaseParticle<dim>::vel_check () const
    {
      return check_vel;
    }

    template <int dim>
    void
    BaseParticle<dim>::set_vel_check (bool new_vel_check)
    {
      check_vel = new_vel_check;
    }

    template class BaseParticle<2>;
    template class BaseParticle<3>;
  }
}
