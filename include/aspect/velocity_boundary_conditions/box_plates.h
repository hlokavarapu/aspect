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


#ifndef __aspect__velocity_boundary_conditions_box_plates_h
#define __aspect__velocity_boundary_conditions_box_plates_h

#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>

namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    using namespace dealii;

    namespace internal
    {
      /**
       * AsciiDataLookup reads in files containing input data
         in ascii format. Note the required format of the
         input data: The first lines may contain any number of comments
         if they begin with '#', but one of these lines needs to
         contain the number of grid points in each dimension as
         for example '# POINTS: 3 3'. The order of the columns
         has to be 'coordinates data' with @p dim coordinate columns
         and @p components data columns. Note that the data in the input
         files need to be sorted in a specific order:
         the first coordinate needs to ascend first,
         followed by the second and so on in order to
         assign the correct data to the prescribed coordinates.
       */
      template <int dim>
      class BoxPlatesLookup
      {
        public:
          BoxPlatesLookup(const std::string &filename,
                          const unsigned int components,
                          const double time_scale_factor,
                          const double velocity_scale_factor);

          /**
           * Loads a data text file. Throws an exception if the file does not exist,
           * if the data file format is incorrect or if the file grid changes over model runtime.
           */
          void
          load_file(const std::string &filename,
                    const double time);

          /**
           * Returns the computed velocity
           * in cartesian coordinates.
           *
           * @param position The current position to compute the data (velocity, temperature, etc.)
           * @param component Number of the component that is requested
           */
          double
          get_data(const Point<dim> &position,
                   const unsigned int component) const;

        private:
          /**
           * The number of data components read in (=columns in the data file).
           */
          const unsigned int components;

          /**
           * Interpolation functions to access the data.
           */
          std::vector<Functions::InterpolatedUniformGridData<dim> *> data;

          /**
           * Model size
           */
          std_cxx11::array<std::pair<double,double>,dim> grid_extent;

          /**
           * Number of points in the data grid.
           */
          TableIndices<dim> table_points;

          /**
           * Scales the data times by a scalar factor. Can be
           * used to transform the unit of the data.
           */
          const double time_scale_factor;

          /**
           * Scales the data boundary velocity by a scalar factor. Can be
           * used to transform the unit of the data.
           */
          const double velocity_scale_factor;

          /**
           * Any plate velocity consists of a velocity and the rotation
           * of the associated plate in the used map projection
           */
          typedef std::pair<Tensor<1,dim+1>,double> plate_velocity;

          typedef std::map<unsigned char,plate_velocity > velocity_map;

          /**
           * Table for the plate velocities of all times.
           */
          std::vector<std::pair<double, velocity_map> > velocity_values;

          /**
           * Computes the table indices of each entry in the input data file.
           * The index depends on dim, grid_dim and the number of components.
           */
          TableIndices<dim>
          compute_table_indices(const unsigned int line) const;
      };
    }

    /**
     * A class that implements prescribed velocity boundary conditions
     * determined from input files in the form of a single velocity file
     * and several plate id files (one per time step).
     *
     * @ingroup VelocityBoundaryConditionsModels
     */
    template <int dim>
    class BoxPlates : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        BoxPlates ();

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function returns the value from the lookup class.
         */
        Tensor<1,dim>
        boundary_velocity (const Point<dim> &position) const;

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        virtual
        void
        initialize ();

        /**
         * A function that is called at the beginning of each time step. For
         * the current plugin, this function loads the next velocity files if
         * necessary and outputs a warning if the end of the set of velocity
         * files is reached.
         */
        virtual
        void
        update ();

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * A variable that stores the current time of the simulation, but
         * relative to the velocity_file_start_time.
         */
        double time_relative_to_vel_file_start_time;

        /**
         * A variable that stores the currently used velocity file of a
         * series. It gets updated if necessary by set_current_time.
         */
        int  current_time_step;

        /**
         * Time at which the velocity file with number 0 is loaded.
         * Previous to this time, a no-slip boundary condition is assumed.
         */
        double velocity_file_start_time;

        /**
         * Time at which the velocity file with number 0 shall be loaded.
         * Previous to this time, a no-slip boundary condition is assumed.
         */
        double first_velocity_file_time;

        /**
         * Directory in which the velocity files are present.
         */
        std::string data_directory;

        /**
         * Filename of velocity file.
         */
        std::string velocity_file_name;

        /**
         * First part of filename of id files. The files have to have
         * the pattern velocity_file_name.n.gpml where n is the number of the
         * current timestep (starts from 0).
         */
        std::string id_file_names;

        /**
         * Time in model units (depends on other model inputs) between two
         * velocity files.
         */
        double time_step;

        /**
         * Distance between two grid points in x direction.
         */
        double x_step;

        /**
         * Distance between two grid points in y direction.
         */
        double y_step;

        /**
         * Weight between velocity file n and n+1 while the current time is
         * between the two values t(n) and t(n+1).
         */
        double time_weight;

        /**
         * State whether we have time_dependent boundary conditions. Switched
         * off after finding no more velocity files to suppress attempts to
         * read in new files.
         */
        bool time_dependent;

        /**
         * Scale the velocity boundary condition by a scalar factor.
         */
        double scale_factor;

        /**
         * Scale the time steps of the velocity file by a scalar factor. Can be
         * used to transform the unit of the time (if they are not
         * specified in the default unit (s or yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double time_scale_factor;


        /**
         * Determines the width of the velocity interpolation zone around the
         * current point. Currently equals the distance between evaluation
         * point and velocity data point that is still included in the
         * interpolation. The weighting of the points currently only accounts
         * for the surface area a single data point is covering ('moving
         * window' interpolation without distance weighting).
         */
        double interpolation_width;

        /**
         * Pointer to an object that reads and processes data we get from
         * the data files.
         */
        std_cxx1x::shared_ptr<internal::BoxPlatesLookup<dim-1> > lookup;

        /**
         * Handles the update of the velocity data in lookup.
         */
        void
        update_velocity_data ();

        /**
         * Handles settings and user notification in case the time-dependent
         * part of the boundary condition is over.
         */
        void
        end_time_dependence (const int timestep);

        /**
         * Create a filename out of the name template.
         */
        std::string
        create_filename (const int timestep) const;
    };
  }
}


#endif
