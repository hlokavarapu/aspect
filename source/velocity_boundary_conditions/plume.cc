/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/velocity_boundary_conditions/plume.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/parameter_handler.h>
#include <fstream>
#include <iostream>
#include <utility>
#include <limits>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {

// ------------------------------ Box -------------------

    template <int dim>
    void
    Plume<dim>::
    initialize ()
    {
      // verify that the geometry is in fact a box since only
      // for this geometry do we know for sure what boundary indicators it
      // uses and what they mean
      Assert (dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model())
              != 0,
              ExcMessage ("This boundary model is only implemented if the geometry is "
                          "in fact a box."));

      lookup.reset(new BoundaryTemperature::internal::PlumeLookup<dim>(data_directory+plume_file_name,
                                                                       this->get_pcout()));
    }

    template <int dim>
    void
    Plume<dim>::
    update ()
    {
      Point<dim> extents = (dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model()))->get_extents();
      extents(dim-1) = 0;
      extents /= 2.0;
      plume_position = lookup->plume_position(this->get_time()) + extents;
    }

    template <int dim>
    Tensor<1,dim>
    Plume<dim>::
    boundary_velocity (const Point<dim> &position) const
    {
      // T=T_0*exp-(r/r_0)**2
      const double vz = V0 * std::exp(-std::pow((position-plume_position).norm()/R0,2));

      Tensor<1,dim> velocity;
      velocity[dim-1] = vz;

      return velocity;
    }

    template <int dim>
    void
    Plume<dim>::declare_parameters (ParameterHandler &prm)
    {
//      prm.enter_subsection("Boundary temperature model");
//      {
        prm.enter_subsection("Plume");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/boundary-temperature/plume/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the model data. This path "
                             "may either be absolute (if starting with a '/') or relative to "
                             "the current directory. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Plume position file name", "Tristan.sur",
                             Patterns::Anything (),
                             "The file name of the plume position data.");
          prm.declare_entry ("Inflow velocity", "0",
                             Patterns::Double (),
                             "Magnitude of the velocity inflow. Units: K.");
          prm.declare_entry ("Radius", "0",
                             Patterns::Double (),
                             "Radius of the anomaly. Units: m.");
//        }
//        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Plume<dim>::parse_parameters (ParameterHandler &prm)
    {
      // Query the unit system for time since we may have to convert below.
      // Note that we can't use this->convert_output_to_years() since this
      // requires the SimulatorAccess base object to have been initialized,
      // but this hasn't happened yet when we get into this function.
      const bool
      use_years_instead_of_seconds
        = prm.get_bool ("Use years in output instead of seconds");


  //    prm.enter_subsection("Boundary temperature model");
  //    {
        prm.enter_subsection("Plume");
        {
          // Get the path to the data files. If it contains a reference
           // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
           // as a #define
           data_directory        = prm.get ("Data directory");
           {
             const std::string      subst_text = "$ASPECT_SOURCE_DIR";
             std::string::size_type position;
             while (position = data_directory.find (subst_text),  position!=std::string::npos)
               data_directory.replace (data_directory.begin()+position,
                                       data_directory.begin()+position+subst_text.size(),
                                       ASPECT_SOURCE_DIR);
           }

           plume_file_name    = prm.get ("Plume position file name");
           V0 = prm.get_double ("Inflow velocity");
           R0 = prm.get_double ("Radius");

           if (use_years_instead_of_seconds == true)
             {
               V0 /= year_in_seconds;
             }
        }
        prm.leave_subsection ();
 //     }
 //     prm.leave_subsection ();
    }


  }
}

// explicit instantiations
namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(Plume,
                                                 "plume",
                                                 "A model in which the velocity is chosen to be zero "
                                                 "everywhere plus a gaussian inflow at the position "
                                                 "of a plume. The plume position is read in from a "
                                                 "file and updated over time.")
  }
}
