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
#include <aspect/velocity_boundary_conditions/box_plates.h>

#include <aspect/geometry_model/box.h>
#include <aspect/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {

    namespace internal
    {
      template <int dim>
      BoxPlatesLookup<dim>::BoxPlatesLookup(const std::string &filename,
                                   const Point<dim> &grid_extent_,
                                   const double dx,
                                   const double dy,
                                   const double interpolation_width_,
                                   const ConditionalOStream &pcout)
        :
        ids(0,0),
        old_ids(0,0),
        id_values (&ids),
        old_id_values (&old_ids),
        current_time_index(0),
        delta_x(dx),
        delta_y(dy),
        grid_extent(grid_extent_),
        nx((grid_extent_(0) / dx)),
        ny((dim == 3) ?
            (grid_extent_(1) / dy)
            : nx),
        interpolation_width(interpolation_width_)
      {

        pcout << std::endl << "   Loading Plate velocity file: "
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
        const double cmyr_si = year_in_seconds * 100;

        std::string temp;
        std::ifstream in(filename.c_str(), std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't open file <") + filename));

        getline(in, temp); // eat first line


        double old_time = -1;
        std::map<unsigned char, plate_velocity > velocity_slice;

        while (!in.eof())
          {
            double start_time,end_time,vx,vy,omega;
            unsigned char plate;
            in >> start_time >> end_time >> plate >> vx >> vy >> omega;

            getline(in, temp);
            if (in.eof())
              {
                velocity_values.push_back(velocity_slice);
                break;
              }

            plate_velocity velocity;
            switch(dim)
            {
            case 2:
              velocity.velocity[0] = vx / cmyr_si;
              break;
            case 3:
              velocity.velocity[0] = vx / cmyr_si;
              velocity.velocity[1] = vy / cmyr_si;
              break;
            default:
              AssertThrow(false,ExcNotImplemented());
              break;
            }

            // omega is in [cm/yr/km] scale it to [1/s]
            velocity.rotation = omega / cmyr_si / 1000;

            if (start_time > old_time + 1e-16)
              {
                if (times.size() > 0)
                  {
                  velocity_values.push_back(velocity_slice);
                  velocity_slice = std::map<unsigned char, plate_velocity > ();
                  }

                // Correct times to years instead of Ma
                // (start_time-time)*Myr_in_seconds
                times.push_back(start_time*Myr_in_seconds);

                old_time = start_time;
              }

            velocity_slice.insert(std::pair<unsigned char,plate_velocity >
            (plate,velocity));

          }

        // We read in the velocity values from today backwards, therefore if
        // we start the model run in the past we start at the last velocity_values
        // slice and move to the first element over time
        current_time_index = velocity_values.size()-2;
      }

      template <int dim>
      void BoxPlatesLookup<dim>::screen_output(const ConditionalOStream &pcout) const
      {
        std::ostringstream output;

        output << std::setprecision (3) << std::setw(3) << std::fixed << std::endl
               << "   Set up Box plates boundary velocity module."  << std::endl
               << "   The grid extents to: " << grid_extent
               << std::endl;

        pcout << output.str();
      }

      template <int dim>
      void BoxPlatesLookup<dim>::update(const double time, const double first_velocity_file_time,
          const ConditionalOStream &pcout)
      {
        if (current_time_index > 0)
          {
          if (first_velocity_file_time - time < times[current_time_index])
            {
              --current_time_index;
              pcout << std::endl << "Set current time index to: " << current_time_index << std::endl;
            }
          velocity_time_weight = (time - (times.back() - times[current_time_index+1])) / (times[current_time_index+1] - times[current_time_index]);
          pcout << std::endl << "New velocity time weight is: " << velocity_time_weight;
          }

      }

      template <int dim>
      bool
      BoxPlatesLookup<dim>::fexists(const std::string &filename) const
      {
        std::ifstream ifile(filename.c_str());
        return !(!ifile); // only in c++11 you can convert to bool directly
      }

      template <int dim>
      void
      BoxPlatesLookup<dim>::load_file(const std::string &filename, const ConditionalOStream &pcout)
      {
        pcout << std::endl << "   Loading Box plates boundary velocity file "
              << filename << "." << std::endl << std::endl;

        // Check whether file exists, we do not want to throw
        // an exception in case it does not, because it could be by purpose
        // (i.e. the end of the boundary condition is reached)
        AssertThrow (fexists(filename),
                     ExcMessage (std::string("Box plates file <")
                                 +
                                 filename
                                 +
                                 "> not found!"));

        std::string temp;
        std::ifstream in(filename.c_str(), std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't find velocities. Is file native gpml format for velocities?")));

        // swap pointers to old and new values, we overwrite the old ones
        // and the new ones become the old ones
        std::swap (id_values, old_id_values);

        (*id_values).reinit(nx,ny);


        unsigned int i = 0;
        char sep;
        while (!in.eof())
          {
            unsigned char plate_id;
            const double cmyr_si = 3.1557e9;

            in >> plate_id;

            if (in.eof())
              break;

            const unsigned int idx_x = i % nx;
            const unsigned int idx_y = i / nx;

            (*id_values)[idx_x][idx_y] = plate_id;

            i++;
          }

        AssertThrow(i == nx*ny,
                    ExcMessage (std::string("Number of read in points does not match number of assumed points. File corrupted?")));
      }

      template <int dim>
      Tensor<1,dim>
      BoxPlatesLookup<dim>::surface_velocity(const Point<dim> &position,
                                      const double time_weight) const
      {

        Tensor<1,dim> surf_vel;

        unsigned int spatial_index[2] = {0,0};
        calculate_spatial_index(spatial_index,position);

        // always use closest point, ensured by the check_termination = false setting
        double n_interpolation_weight = add_interpolation_point(surf_vel,position,spatial_index,time_weight,false);

        // If interpolation is requested, loop over points in theta, phi direction until we exit a circle of
        // interpolation_width radius around the evaluation point position
        if (interpolation_width > 0)
          {
            unsigned int i = 0;
            unsigned int j = 1;
            do
              {
                do
                  {
                    const unsigned int idx[2][2] = { {(unsigned int)(spatial_index[0]+i), (unsigned int)(spatial_index[1]+j)},
                      {(unsigned int)(spatial_index[0]-i), (unsigned int)(spatial_index[1]-j)}
                    };

                    const double weight_1 = add_interpolation_point(surf_vel,position,idx[0],time_weight,true);
                    const double weight_2 = add_interpolation_point(surf_vel,position,idx[1],time_weight,true);

                    // termination criterion, our interpolation points are outside of the defined circle
                    if ((std::abs(weight_1 + 1) < 1e-14) || (std::abs(weight_2 + 1) < 1e-14))
                      {
                        if (j == 0)
                          {
                            // If we are outside of the interpolation circle even without extending phi (j), we have
                            // visited all possible interpolation points. Return current surf_vel.
                            return surf_vel / n_interpolation_weight;
                          }
                        else
                          {
                            // We are outside of the interpolation circle, but resetting j to 0 (phi to original index)
                            // and moving forward in theta direction (i) might show us more interpolation points.
                            break;
                          }
                      }

                    n_interpolation_weight += weight_1;
                    n_interpolation_weight += weight_2;

                    j++;
                  }
                while (j < id_values->n_cols() / 2);
                j = 0;
                i++;
              }
            while (i < id_values->n_rows() / 2);
          }

        // We have interpolated over the whole dataset of the provided velocities, which is ok.
        // Velocity will be constant for all evaluation points in this case.

        return surf_vel / n_interpolation_weight;
      }

      template <int dim>
      double
      BoxPlatesLookup<dim>::add_interpolation_point(Tensor<1,dim>       &surf_vel,
                                             const Point<dim> &position,
                                             const unsigned int          grid_index[2],
                                             const double       time_weight,
                                             const bool         check_termination) const
      {
        const double distance = (position - get_grid_point_position(grid_index[0],
                                                                   grid_index[1])).norm();

        // termination criterion, our interpolation points are outside of the defined circle
        if (check_termination && (distance > interpolation_width))
          {
            return -1;
          }

        const double point_weight = 1.0;

        const unsigned char plate_id = (*id_values)[grid_index[0]][grid_index[1]];
        const unsigned char old_plate_id = (*old_id_values)[grid_index[0]][grid_index[1]];

        Tensor<1,dim> velocity = velocity_values[current_time_index].find(plate_id)->second.velocity;
        Tensor<1,dim> old_velocity = velocity_values[current_time_index+1].find(plate_id)->second.velocity;

        const double omega = velocity_values[current_time_index].find(plate_id)->second.rotation;
        const double old_omega = velocity_values[current_time_index+1].find(plate_id)->second.rotation;


        Tensor<1,dim> rotation_velocity;
        Tensor<1,dim> old_rotation_velocity;

        if (dim == 3)
          {
            rotation_velocity[0] = -omega*position[1];
            rotation_velocity[1] = omega*position[0];
            old_rotation_velocity[0] = -old_omega*position[1];
            old_rotation_velocity[1] = old_omega*position[0];
          }
        else if (dim == 2)
          {
            rotation_velocity[0] = -omega*position[1];
            old_rotation_velocity[0] = -old_omega*position[1];
          }

        velocity += rotation_velocity;
        old_velocity += old_rotation_velocity;
        //            phase_func = 0.5*(1.0 + std::tanh(pressure_deviation / pressure_width));

        const double depth_weight = 0.5*(1.0 + std::tanh((position(dim-1)-460000)/50000));
        surf_vel += point_weight * depth_weight * velocity_time_weight * velocity;

        if (time_weight < 1.0 - 1e-7)
          surf_vel += point_weight * depth_weight * (1-velocity_time_weight) * old_velocity;

        return point_weight;
      }


      template <int dim>
      Point<dim>
      BoxPlatesLookup<dim>::get_grid_point_position(const unsigned int x_index,
                                                    const unsigned int y_index) const
      {
        Point<dim> spherical_position;
        switch (dim)
        {
        case 2:
          return Point<dim> (x_index * delta_x, 0);
          break;
        case 3:
          return Point<dim> (x_index * delta_x, y_index * delta_y, 0);
          break;
        default:
          AssertThrow(false,
              ExcNotImplemented());
          break;
        }
        return Point<dim>();
      }

      template <int dim>
      void
      BoxPlatesLookup<dim>::calculate_spatial_index(unsigned int *index,
                                             const Point<dim> &position) const
      {
        switch (dim)
        {
        case 2:
          index[0] = std::min(static_cast<unsigned int> (position[0]/delta_x),nx-1);
          index[1] = 0;
          break;
        case 3:
          index[0] = std::min(static_cast<unsigned int> (position[0]/delta_x),nx-1);
          index[1] = std::min(static_cast<unsigned int> (position[1]/delta_y),ny-1);
          break;
        default:
          AssertThrow(false,ExcNotImplemented());
        }
      }
    }


    template <int dim>
    BoxPlates<dim>::BoxPlates ()
      :
      time_relative_to_vel_file_start_time(0.0),
      current_time_step(-2),
      velocity_file_start_time(0.0),
      time_step(0.0),
      time_weight(0.0),
      time_dependent(true),
      interpolation_width(0.0),
      lookup()
    {}


    template <int dim>
    void
    BoxPlates<dim>::initialize ()
    {
      Assert (dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model()) != 0,
              ExcMessage ("This boundary condition can only be used if the geometry "
                          "is a spherical shell."));

      const GeometryModel::Box<dim> *geometry_model = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model());
      const Point<dim> grid_extent = geometry_model->get_extents();

      lookup.reset(new internal::BoxPlatesLookup<dim>(data_directory+velocity_file_name,
                                                 grid_extent,
                                                 x_step,
                                                 y_step,
                                                 interpolation_width,
                                                 this->get_pcout()));


    }



    template <int dim>
    std::string
    BoxPlates<dim>::create_filename (const int timestep) const
    {
      std::string templ = data_directory+id_file_names;
      const int size = templ.length();
      char *filename = (char *) (malloc ((size + 10) * sizeof(char)));
      snprintf (filename, size + 10, templ.c_str (),
                (const unsigned int) ((first_velocity_file_time - timestep * time_step)/(1e6*year_in_seconds)));
      std::string str_filename (filename);
      free (filename);
      return str_filename;
    }


    template <int dim>
    void
    BoxPlates<dim>::update ()
    {
      Interface<dim>::update ();

      time_relative_to_vel_file_start_time = this->get_time() - velocity_file_start_time;

      lookup->update(time_relative_to_vel_file_start_time,first_velocity_file_time,this->get_pcout());

      // If the boundary condition is constant, switch off time_dependence end leave function.
      // This also sets time_weight to 1.0
      if ((create_filename (current_time_step) == create_filename (current_time_step+1)) && time_dependent)
        {
          lookup->screen_output(this->get_pcout());

          lookup->load_file(create_filename (current_time_step),
                             this->get_pcout());
          end_time_dependence (current_time_step);
          return;
        }

      if (time_dependent && (time_relative_to_vel_file_start_time >= 0.0))
        {
          if (current_time_step < 0)
            lookup->screen_output (this->get_pcout());

          if (static_cast<int> (time_relative_to_vel_file_start_time / time_step) > current_time_step)
            {
              update_velocity_data();
            }

          time_weight = time_relative_to_vel_file_start_time / time_step - current_time_step;

          Assert ((0 <= time_weight) && (time_weight <= 1),
                  ExcMessage (
                    "Error in set_current_time. Time_weight has to be in [0,1]"));
        }
    }

    template <int dim>
    void
    BoxPlates<dim>::update_velocity_data ()
    {
      const int old_time_step = current_time_step;
      current_time_step =
        static_cast<unsigned int> (time_relative_to_vel_file_start_time / time_step);
      // Load next velocity file for interpolation
      // If the time step was large enough to move forward more
      // then one velocity file, we need to load both current files
      // to stay accurate in interpolation
      if (current_time_step > old_time_step + 1)
        try
          {
            lookup->load_file (create_filename (current_time_step),
                               this->get_pcout());
          }
        catch (...)
          // If loading current_time_step failed, end time dependent part with old_time_step.
          {
            try
              {
                end_time_dependence (old_time_step);
              }
            catch (...)
              {
                // If loading the old file fails (e.g. there was no old file), cancel the model run.
                // We might get here, if the time step is so large that step t is before the whole boundary condition
                // while step t+1 is already behind all files in time.
                AssertThrow (false,
                             ExcMessage (
                               "Loading new and old velocity file did not succeed. "
                               "Maybe the time step was so large we jumped over all files "
                               "or the files were removed during the model run. "
                               "Another possible way here is to restart a model with "
                               "previously time-dependent boundary condition after the "
                               "last file was already read. Aspect has no way to find the "
                               "last readable file from the current model time. Please "
                               "prescribe the last velocity file manually in such a case. "
                               "Cancelling calculation."));
              }
          }

      // Now load the velocity file. This part is the main purpose of this function.
      try
        {
          lookup->load_file (create_filename (current_time_step+1),
                             this->get_pcout());
        }

      // If loading current_time_step + 1 failed, end time dependent part with current_time_step.
      // We do not need to check for success here, because current_time_step was guaranteed to be
      // at least tried to be loaded before, and if it fails, it should have done before (except from
      // hard drive errors, in which case the exception is the right thing to be thrown).

      catch (...)
        {
          end_time_dependence (current_time_step);
        }
    }

    template <int dim>
    void
    BoxPlates<dim>::end_time_dependence (const int timestep)
    {
      // Next velocity file not found --> Constant velocities
      // by simply loading the old file twice
      lookup->load_file (create_filename (timestep), this->get_pcout());
      // no longer consider the problem time dependent from here on out
      // this cancels all attempts to read files at the next time steps
      time_dependent = false;
      // this cancels the time interpolation in lookup
      time_weight = 1.0;
      // Give warning if first processor
      this->get_pcout() << std::endl
                        << "   Loading new velocity file did not succeed." << std::endl
                        << "   Assuming constant boundary conditions for rest of model run."
                        << std::endl << std::endl;
    }

    template <int dim>
    Tensor<1,dim>
    BoxPlates<dim>::
    boundary_velocity (const Point<dim> &position) const
    {
      if (time_relative_to_vel_file_start_time >= 0.0)
        return scale_factor * lookup->surface_velocity(position,time_weight);
      else
        return Tensor<1,dim> ();
    }

    template <int dim>
    void
    BoxPlates<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Boundary velocity model");
      {
        prm.enter_subsection ("Box plates model");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/velocity-boundary-conditions/tristan_plume/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the model data. This path "
                             "may either be absolute (if starting with a '/') or relative to "
                             "the current directory. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Id file names", "plates_%d",
                             Patterns::Anything (),
                             "The file name of the id data. Provide file in format: "
                             "(Velocity file name).\\%d.gpml where \\%d is any sprintf integer "
                             "qualifier, specifying the format of the current file number.");
          prm.declare_entry ("Velocity file name", "velocities",
                             Patterns::Anything (),
                             "The file name of the material data. Provide file in format: "
                             "(Velocity file name).\\%d.gpml where \\%d is any sprintf integer "
                             "qualifier, specifying the format of the current file number.");
          prm.declare_entry ("Time step", "1e6",
                             Patterns::Double (0),
                             "Time step between following velocity files. "
                             "Depending on the setting of the global 'Use years in output instead of seconds' flag "
                             "in the input file, this number is either interpreted as seconds or as years. "
                             "The default is one million, i.e., either one million seconds or one million years.");
          prm.declare_entry ("X grid spacing", "20e3",
                             Patterns::Double (0),
                             "Distance between two grid points in x direction.");
          prm.declare_entry ("Y grid spacing", "20e3",
                             Patterns::Double (0),
                             "Distance between two grid points in y direction.");
          prm.declare_entry ("Velocity file start time", "0",
                             Patterns::Double (0),
                             "Time at which the velocity file with number 0 shall be loaded. Previous to this "
                             "time, a no-slip boundary condition is assumed. "
                             "Depending on the setting of the global 'Use years in output instead of seconds' flag "
                             "in the input file, this number is either interpreted as seconds or as years.");
          prm.declare_entry ("First velocity file time", "140e6",
                             Patterns::Double (0),
                             "Time at which the velocity file with number 0 shall be loaded. Previous to this "
                             "time, a no-slip boundary condition is assumed. "
                             "Depending on the setting of the global 'Use years in output instead of seconds' flag "
                             "in the input file, this number is either interpreted as seconds or as years.");
          prm.declare_entry ("Scale factor", "1",
                             Patterns::Double (0),
                             "Scalar factor, which is applied to the boundary velocity. "
                             "You might want to use this to scale the velocities to a "
                             "reference model (e.g. with free-slip boundary) or another "
                             "plate reconstruction.");
          prm.declare_entry ("Interpolation width", "0",
                             Patterns::Double (0),
                             "Determines the width of the velocity interpolation zone around the current point. "
                             "Currently equals the arc distance between evaluation point and velocity data point that "
                             "is still included in the interpolation. The weighting of the points currently only accounts "
                             "for the surface area a single data point is covering ('moving window' interpolation without "
                             "distance weighting).");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    BoxPlates<dim>::parse_parameters (ParameterHandler &prm)
    {
      // Query the unit system for time since we may have to convert below.
      // Note that we can't use this->convert_output_to_years() since this
      // requires the SimulatorAccess base object to have been initialized,
      // but this hasn't happened yet when we get into this function.
      const bool
      use_years_instead_of_seconds
        = prm.get_bool ("Use years in output instead of seconds");

      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("Box plates model");
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

          velocity_file_name    = prm.get ("Velocity file name");
          id_file_names    = prm.get ("Id file names");

          interpolation_width   = prm.get_double ("Interpolation width");
          scale_factor          = prm.get_double ("Scale factor");
          x_step                = prm.get_double ("X grid spacing");
          y_step                = prm.get_double ("Y grid spacing");

          time_step             = prm.get_double ("Time step");
          velocity_file_start_time = prm.get_double ("Velocity file start time");
          first_velocity_file_time = prm.get_double ("First velocity file time");

          if (use_years_instead_of_seconds == true)
            {
              time_step                *= year_in_seconds;
              velocity_file_start_time *= year_in_seconds;
              first_velocity_file_time *= year_in_seconds;
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(BoxPlates,
                                                 "box plates",
                                                 "Implementation of a model in which the boundary "
                                                 "velocity is derived from files that are generated "
                                                 "by the GPlates program.")
  }
}
