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
                                            const unsigned int components,
                                            const double time_scale_factor,
                                            const double velocity_scale_factor)
        :
          components(components),
          data(components),
          time_scale_factor(time_scale_factor),
          velocity_scale_factor(velocity_scale_factor)
      {
        // Check whether file exists, we do not want to throw
        // an exception in case it does not, because it could be by purpose
        // (i.e. the end of the boundary condition is reached)
        AssertThrow (Utilities::fexists(filename),
                     ExcMessage (std::string("Plate velocity file <")
                                 +
                                 filename
                                 +
                                 "> not found!"));

        std::string temp;
        std::ifstream in(filename.c_str(), std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't open file <") + filename));
        getline(in, temp); // eat first line


        double old_time (std::numeric_limits<double>::quiet_NaN());
        velocity_map velocity_slice;
        double start_time,end_time,vx,vy,omega;
        unsigned char plate_id;
        while (in >> start_time >> end_time >> plate_id >> vx >> vy >> omega)
          {
            // scale all properties
            vx *= velocity_scale_factor;
            vy *= velocity_scale_factor;
            // omega is in [velocity/km] scale it to [1/s]
            omega *= velocity_scale_factor / 1000;
            start_time *= time_scale_factor;
            end_time *= time_scale_factor;

            // If we have just read in the next time slice
            // save the last one and start a new one
            if (start_time > old_time + 1e-16)
              {
                velocity_values.push_back(std::make_pair(old_time,velocity_slice));
                velocity_slice = velocity_map();
              }

            // each time will contain a map from plate character to velocity
            plate_velocity velocity;
            switch(dim)
            {
            case 1:
              velocity.first[0] = vx;
              break;
            case 2:
              velocity.first[0] = vx;
              velocity.first[1] = vy;
              break;
            default:
              AssertThrow(false,ExcNotImplemented());
              break;
            }

            velocity.second = omega;

            velocity_slice.insert(std::pair<unsigned char,plate_velocity >
            (plate_id,velocity));

            old_time = start_time;
          }

        //save the last slice
        velocity_values.push_back(std::make_pair(start_time,velocity_slice));
      }

      template <int dim>
      void
      BoxPlatesLookup<dim>::load_file(const std::string &filename,
                                      const double time)
      {
        const double time_until_end = velocity_values.back().first - time;
        unsigned int old_index,next_index;
        double velocity_time_weight;

        if (time_until_end > velocity_values.back().first)
          {
            old_index = velocity_values.size();
            next_index = velocity_values.size();
            velocity_time_weight = 0.0;
          }
        else if (time_until_end < 0.0)
          {
            old_index = 0;
            next_index = 0;
            velocity_time_weight = 1.0;
          }
        else
          {
            for (unsigned int i = velocity_values.size() - 2; i >= 0; i--)
                if (time_until_end >= velocity_values[i].first)
                  {
                    old_index = i + 1;
                    next_index = i;
                    velocity_time_weight = (velocity_values[old_index].first - time_until_end) / (velocity_values[old_index].first - velocity_values[next_index].first);
                    break;
                  }
          }

        const velocity_map old_map = velocity_values[old_index].second;
        const velocity_map next_map = velocity_values[next_index].second;

        Assert((0.0 <= velocity_time_weight) && (1.0 >= velocity_time_weight),
               ExcMessage ("Velocity time weight is wrong"));

        // Check whether file exists, we do not want to throw
        // an exception in case it does not, because it could be by purpose
        // (i.e. the end of the boundary condition is reached)
        AssertThrow (Utilities::fexists(filename),
                     ExcMessage (std::string("Ascii data file <")
                                 +
                                 filename
                                 +
                                 "> not found!"));

        std::ifstream in(filename.c_str(), std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't open data file <"
                                             +
                                             filename
                                             +
                                             ">.")));

        // Read header lines and if necessary reinit tables
        while (in.peek() == '#')
          {
            std::string line;
            getline(in,line);
            std::stringstream linestream(line);
            std::string word;
            while (linestream >> word)
              if (word == "POINTS:")
                for (unsigned int i = 0; i < dim; i++)
                  {
                    unsigned int temp_index;
                    linestream >> temp_index;

                    if (table_points[i] == 0)
                      table_points[i] = temp_index;
                    else
                      AssertThrow (table_points[i] == temp_index,
                                   ExcMessage("The file grid must not change over model runtime. "
                                              "Either you prescribed a conflicting number of points in "
                                              "the input file, or the POINTS comment in your data files "
                                              "is changing between following files."));
                  }
          }

        /**
         * Table for the new data. This peculiar reinit is necessary, because
         * there is no constructor for Table, which takes TableIndices as
         * argument.
         */
        Table<dim,double> data_table;
        data_table.TableBase<dim,double>::reinit(table_points);
        std::vector<Table<dim,double> > data_tables(components+dim,data_table);

        // Read data lines
        unsigned int line = 0;
        double temp_data;

        while (!in.eof())
          {
            Point<dim> position;
             for (unsigned int i = 0; i < dim; i++)
               {
                 if (!(in >> position[i]))
                   break;
                 data_tables[i](compute_table_indices(line)) = position[i];
               }

             char plate_id;
             if (!(in >> plate_id))
               break;

            Tensor<1,dim+1> old_velocity = old_map.find(plate_id)->second.first;
            const double old_omega = old_map.find(plate_id)->second.second;

            Tensor<1,dim+1> velocity;
            double omega;

            // It might happen that the current plate disappears in the next time step
            // In case the plate is not longer there, use the old velocity
            if (old_map.find(plate_id) != old_map.end())
              {
                velocity = next_map.find(plate_id)->second.first;
                omega    = next_map.find(plate_id)->second.second;
              }
            else
              {
                velocity = old_velocity;
                omega    = old_omega;
              }


             Tensor<1,dim+1> rotation_velocity;
             Tensor<1,dim+1> old_rotation_velocity;

             if (dim == 2)
               {
                 rotation_velocity[0] = -omega*position[1];
                 rotation_velocity[1] = omega*position[0];
                 old_rotation_velocity[0] = -old_omega*position[1];
                 old_rotation_velocity[1] = old_omega*position[0];
               }
             else if (dim == 1)
               {
                 rotation_velocity[0] = -omega*position[0];
                 old_rotation_velocity[0] = -old_omega*position[0];
               }

             velocity += rotation_velocity;
             old_velocity += old_rotation_velocity;

             const Tensor<1,dim+1> surface_velocity = velocity_time_weight * velocity
                                                      + (1-velocity_time_weight) * old_velocity;

             for (unsigned int i = 0; i < dim+1; i++)
               data_tables[dim+i](compute_table_indices(line)) = surface_velocity[i];

            line++;

            // TODO: add checks for coordinate ordering in data files
          }

        AssertThrow(line == data_table.n_elements(),
                    ExcMessage (std::string("Number of read in points does not match number of expected points. File corrupted?")));

        std_cxx11::array<unsigned int,dim> table_intervals;

        for (unsigned int i = 0; i < dim; i++)
          {
            table_intervals[i] = table_points[i] - 1;

            TableIndices<dim> idx;
            grid_extent[i].first = data_tables[i](idx);
            idx[i] = table_points[i] - 1;
            grid_extent[i].second = data_tables[i](idx);
          }

        for (unsigned int i = 0; i < components; i++)
          {
            if (data[i])
              delete data[i];
            data[i] = new Functions::InterpolatedUniformGridData<dim> (grid_extent,
                                                                       table_intervals,
                                                                       data_tables[dim+i]);
          }
      }


      template <int dim>
      TableIndices<dim>
      BoxPlatesLookup<dim>::compute_table_indices(const unsigned int line) const
      {
        TableIndices<dim> idx;
        idx[0] = line % table_points[0];
        if (dim >= 2)
          idx[1] = (line / table_points[0]) % table_points[1];
        if (dim == 3)
          idx[2] = line / (table_points[0] * table_points[1]);

        return idx;
      }

      template <int dim>
      double
      BoxPlatesLookup<dim>::get_data(const Point<dim> &position,
                                     const unsigned int component) const
      {
        return data[component]->value(position);
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
    {
      AssertThrow(false,
          ExcMessage("This plugin is currently not functional!!! Debug it before use."));
    }


    template <int dim>
    void
    BoxPlates<dim>::initialize ()
    {
      Assert (dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model()) != 0,
              ExcMessage ("This boundary condition can only be used if the geometry "
                          "is a box."));

      const GeometryModel::Box<dim> *geometry_model = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model());
      const Point<dim> grid_extent = geometry_model->get_extents();

//      lookup.reset(new internal::BoxPlatesLookup<dim>(data_directory+velocity_file_name,
//                                                      grid_extent,
//                                                      x_step,
//                                                      y_step,
//                                                      interpolation_width,
//                                                      this->get_pcout(),
//                                                      year_in_seconds * 1e6,
//                                                      scale_factor));
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

//      lookup->update(time_relative_to_vel_file_start_time);

      // If the boundary condition is constant, switch off time_dependence end leave function.
      // This also sets time_weight to 1.0
      if ((create_filename (current_time_step) == create_filename (current_time_step+1)) && time_dependent)
        {
//          lookup->screen_output(this->get_pcout());
//
//          lookup->load_file(create_filename (current_time_step),
//                             this->get_pcout());
          end_time_dependence (current_time_step);
          return;
        }

      if (time_dependent && (time_relative_to_vel_file_start_time >= 0.0))
        {
//          if (current_time_step < 0)
//            lookup->screen_output (this->get_pcout());

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
                               0.0);
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
                             0.0);
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
      lookup->load_file (create_filename (timestep), 0.0);
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
        return Tensor<1,dim> ();
        //return lookup->surface_velocity(position,time_weight);
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
          prm.declare_entry ("Time scale factor", "1e6",
                             Patterns::Double (0),
                             "Determines the factor applied to the times in the velocity file. The unit is assumed to be"
                             "s or yr depending on the 'Use years in output instead of "
                             "seconds' flag. If you provide time in Ma set this "
                             "factor to 1e6.");
          prm.declare_entry ("Scale factor", "1",
                             Patterns::Double (0),
                             "Scalar factor, which is applied to the boundary velocity. "
                             "You might want to use this to scale the velocities to a "
                             "reference model (e.g. with free-slip boundary) or another "
                             "plate reconstruction.");
          prm.declare_entry ("Interpolation width", "0",
                             Patterns::Double (0),
                             "Determines the width of the velocity interpolation zone around the current point. "
                             "Currently equals the distance between evaluation point and velocity data point that "
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
          id_file_names         = prm.get ("Id file names");

          interpolation_width   = prm.get_double ("Interpolation width");
          scale_factor          = prm.get_double ("Scale factor");
          x_step                = prm.get_double ("X grid spacing");
          y_step                = prm.get_double ("Y grid spacing");

          time_step             = prm.get_double ("Time step");
          velocity_file_start_time = prm.get_double ("Velocity file start time");
          first_velocity_file_time = prm.get_double ("First velocity file time");

          time_scale_factor     = prm.get_double ("Time scale factor");
          scale_factor          = prm.get_double ("Scale factor");


          // If the input and output are used in years, convert it to seconds
          // for internal usage.
          if (this->convert_output_to_years())
            {
              time_step                *= year_in_seconds;
              velocity_file_start_time *= year_in_seconds;
              first_velocity_file_time *= year_in_seconds;
              time_scale_factor        *= year_in_seconds;
              scale_factor             /= year_in_seconds;
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
    namespace internal
    {
      template class BoxPlatesLookup<1>;
      template class BoxPlatesLookup<2>;
      template class BoxPlatesLookup<3>;
    }

    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(BoxPlates,
                                                 "box plates",
                                                 "Implementation of a model in which the boundary "
                                                 "velocity is derived from files that are formatted "
                                                 "as a single velocity file and one plate id file per "
                                                 "time step. Example data files are available in "
                                                 "data/velocity-boundary-conditions/box_plates.")
  }
}
