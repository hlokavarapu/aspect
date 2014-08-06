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

#include <aspect/global.h>
#include <aspect/boundary_temperature/plume.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/parameter_handler.h>
#include <fstream>
#include <iostream>
#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {

    namespace internal
     {
      template <int dim>
       PlumeLookup<dim>::PlumeLookup(const std::string &filename,
                                const ConditionalOStream &pcout)
       {
         pcout << std::endl << "   Loading Plume position file for boundary temperature: "
               << filename << "." << std::endl << std::endl;

         // Check whether file exists, we do not want to throw
         // an exception in case it does not, because it could be by purpose
         // (i.e. the end of the boundary condition is reached)
         AssertThrow (fexists(filename),
                      ExcMessage (std::string("Plume position file <")
                                  +
                                  filename
                                  +
                                  "> not found!"));

         const double Myr_in_seconds = 1e6 * year_in_seconds;
         const double km_in_m = 1e3;

         std::string temp;
         std::ifstream in(filename.c_str(), std::ios::in);
         AssertThrow (in,
                      ExcMessage (std::string("Couldn't open file <") + filename));

         unsigned int i = 0;
         double start_time;
         while (!in.eof())
           {
             double time,x,y;
             in >> time >> x >> y;

             if (i == 0)
               start_time = time;

             getline(in, temp);
             if (in.eof())
               break;

             Point<dim> position;
             switch(dim)
             {
             case 2:
               position = {x,0};
               break;
             case 3:
               position = {x,y,0};
               break;
             default:
               AssertThrow(false,ExcNotImplemented());
               break;
             }


             plume_positions.push_back(position*km_in_m);

             // Correct times to years instead of Ma
             times.push_back((start_time-time)*Myr_in_seconds);

             i++;
           }
         Assert(i==plume_positions.size(), ExcMessage("Plume positions table size not consistent with number of lines in file."));
         Assert(i==times.size(), ExcMessage("Times table size not consistent with number of lines in file."));
       }

      template <int dim>
       bool
       PlumeLookup<dim>::fexists(const std::string &filename) const
       {
         std::ifstream ifile(filename.c_str());
         return !(!ifile); // only in c++11 you can convert to bool directly
       }

       template <int dim>
       Point<dim>
       PlumeLookup<dim>::plume_position(const double time) const
       {
         Point<dim> position;
         if (time <= times.front())
           return plume_positions.front();
         else if (time >= times.back())
           return plume_positions.back();
         else
           {
             for (unsigned int i = 0; i < times.size() - 1; i++)
               {
                 if ((time > times[i])
                     && time < times[i+1])
                   {
                     const double timestep = times[i+1]-times[i];

                     position = ((time - times[i]) / timestep) * plume_positions[i]
                         + ((times[i+1] - time) / timestep) * plume_positions[i+1];
                   }
               }
           }
         return position;
       }
     }

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

      lookup.reset(new internal::PlumeLookup<dim>(data_directory+plume_file_name,
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
    double
    Plume<dim>::
    temperature (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &location) const
    {

      Assert (boundary_indicator<2*dim, ExcMessage ("Unknown boundary indicator."));

      double perturbation = 0.0;
      unsigned int bottom = 2;
      if (dim == 3)
        bottom = 4;

      if (boundary_indicator == bottom)
        // T=T_0*exp-(r/r_0)**2
        perturbation = T0 * std::exp(-std::pow((location-plume_position).norm()/R0,2));

        return temperature_[boundary_indicator] + perturbation;
    }


    template <int dim>
    double
    Plume<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
        return *std::min_element(temperature_, temperature_+2*dim);
      else
        {
          double min = maximal_temperature(fixed_boundary_ids);
          for (typename std::set<types::boundary_id>::const_iterator
               p = fixed_boundary_ids.begin();
               p != fixed_boundary_ids.end(); ++p)
            if (p != fixed_boundary_ids.end())
              min = std::min(min,temperature_[*p]);
          return min;
        }
    }



    template <int dim>
    double
    Plume<dim>::
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
        return *std::max_element(temperature_, temperature_+2*dim);
      else
        {
          double max = -std::numeric_limits<double>::max();
          for (typename std::set<types::boundary_id>::const_iterator
               p = fixed_boundary_ids.begin();
               p != fixed_boundary_ids.end(); ++p)
            if (p != fixed_boundary_ids.end())
              max = std::max(max,temperature_[*p]);
          return max;
        }
    }

    template <int dim>
    void
    Plume<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
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
          prm.declare_entry ("Magnitude", "0",
                             Patterns::Double (),
                             "Magnitude of the temperature anomaly. Units: K.");
          prm.declare_entry ("Radius", "0",
                             Patterns::Double (),
                             "Radius of the temperature anomaly. Units: m.");
          prm.declare_entry ("Left temperature", "1",
                             Patterns::Double (),
                             "Temperature at the left boundary (at minimal x-value). Units: K.");
          prm.declare_entry ("Right temperature", "0",
                             Patterns::Double (),
                             "Temperature at the right boundary (at maximal x-value). Units: K.");
          prm.declare_entry ("Bottom temperature", "0",
                             Patterns::Double (),
                             "Temperature at the bottom boundary (at minimal z-value). Units: K.");
          prm.declare_entry ("Top temperature", "0",
                             Patterns::Double (),
                             "Temperature at the top boundary (at maximal x-value). Units: K.");
          if (dim==3)
            {
              prm.declare_entry ("Front temperature", "0",
                                 Patterns::Double (),
                                 "Temperature at the front boundary (at minimal y-value). Units: K.");
              prm.declare_entry ("Back temperature", "0",
                                 Patterns::Double (),
                                 "Temperature at the back boundary (at maximal y-value). Units: K.");
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Plume<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
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
           T0 = prm.get_double ("Magnitude");
           R0 = prm.get_double ("Radius");
          switch (dim)
            {
              case 2:
                temperature_[0] = prm.get_double ("Left temperature");
                temperature_[1] = prm.get_double ("Right temperature");
                temperature_[2] = prm.get_double ("Bottom temperature");
                temperature_[3] = prm.get_double ("Top temperature");
                break;

              case 3:
                temperature_[0] = prm.get_double ("Left temperature");
                temperature_[1] = prm.get_double ("Right temperature");
                temperature_[2] = prm.get_double ("Front temperature");
                temperature_[3] = prm.get_double ("Back temperature");
                temperature_[4] = prm.get_double ("Bottom temperature");
                temperature_[5] = prm.get_double ("Top temperature");
                break;

              default:
                Assert (false, ExcNotImplemented());
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(Plume,
                                               "plume",
                                               "A model in which the temperature is chosen constant on "
                                               "all the sides of a box.")
  }
}
