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


namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    using namespace dealii;


    namespace internal
    {
      /**
       * GPlatesLookup handles all kinds of tasks around looking up a certain
       * velocity boundary condition from a gplates .gpml file. This class
       * keeps around the contents of two sets of files, corresponding to two
       * instances in time where GPlates provides us with data; the boundary
       * values at one particular time are interpolated between the two
       * currently loaded data sets.
       */
      template <int dim>
      class BoxPlatesLookup
      {
        public:

          /**
           * Initialize all members and the two pointers referring to the
           * actual velocities. Also calculates any necessary rotation
           * parameters for a 2D model.
           */
          BoxPlatesLookup(const std::string &filename,
                          const Point<dim> &grid_extent_,
                          const double dx,
                          const double dy,
                          const double interpolation_width_,
                          const ConditionalOStream &pcout);

          /**
           * Outputs the GPlates module information at model start. Need to be
           * separated from Constructor because at construction time the
           * SimulatorAccess is not initialized and only Rank 0 should give
           * the screen output.
           */
          void screen_output(const ConditionalOStream &pcout) const;

          /**
           * Check whether a file named @p filename exists.
           */
          bool fexists(const std::string &filename) const;

          /**
           * Loads a gplates .gpml velocity file. Throws an exception if the
           * file does not exist.
           */
          void load_file(const std::string &filename,
                         const ConditionalOStream &pcout);

          /**
           * Called once per timestep. Updates time_index if necessary.
           */
          void update(const double time,
                      const double first_velocity_file_time,
                      const ConditionalOStream &pcout);

          /**
           * Returns the computed surface velocity in cartesian coordinates.
           * Takes as input the position and current time weight.
           *
           * @param position The current position to compute velocity
           * @param time_weight A weighting between the two current timesteps
           * n and n+1
           */
          Tensor<1,dim> surface_velocity(const Point<dim> &position,
                                         const double time_weight) const;

        private:

          /**
           * Tables which contain the velocities
           */
          dealii::Table<2,unsigned char > ids;
          dealii::Table<2,unsigned char > old_ids;

          /**
           * Pointers to the actual tables. Used to avoid unnecessary copying
           * of values. These pointers point to either ids or
           * old_ids.
           */
          dealii::Table<2,unsigned char > *id_values;
          dealii::Table<2,unsigned char > *old_id_values;

          /**
           *
           */
          std::vector<double> times;
          unsigned int current_time_index;

          struct plate_velocity
          {
              Tensor<1,dim> velocity;
              double rotation;
          };

          /**
           * Table for the data point positions.
           */
          std::vector<std::map<unsigned char,plate_velocity > > velocity_values;

          double velocity_time_weight;

          /**
           * Distances between adjacent point in the x/y grid
           */
          const double delta_x;
          const double delta_y;

          const Point<dim> grid_extent;

          const unsigned int nx;
          const unsigned int ny;


          /**
           * Determines the width of the velocity interpolation zone around
           * the current point. Currently equals the arc distance between
           * evaluation point and velocity data point that is still included
           * in the interpolation. The weighting of the points currently only
           * accounts for the surface area a single data point is covering
           * ('moving window' interpolation without distance weighting).
           */
          const double interpolation_width;

          /**
           * calculates the index given a certain position
           *
           * @param index Reference to the index field, which is modified.
           * @param position Input position, which is converted to spatial
           * index
           */
          void
          calculate_spatial_index(unsigned int *index,
                                  const Point<dim> &position) const;

          /**
           * This function adds a certain data point to the interpolated
           * surface velocity at this evaluation point. This includes
           * calculating the interpolation weight and the rotation of the
           * velocity to the evaluation point position (the velocity need to
           * be tangential to the surface).
           */
          double
          add_interpolation_point(Tensor<1,dim>       &surf_vel,
                                  const Point<dim> &position,
                                  const unsigned int          spatial_index[2],
                                  const double       time_weight,
                                  const bool         check_termination) const;


          /**
           * Returns the position (cartesian or spherical depending on last
           * argument) of a data point with a given theta,phi index.
           */
          Point<dim>
          get_grid_point_position(const unsigned int x_index,
                                  const unsigned int y_index) const;

      };
    }


    /**
     * A class that implements prescribed velocity boundary conditions
     * determined from a GPlates input files.
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
         * current class, this function returns value from gplates.
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
         * Time at which the velocity file with number 0 shall be loaded.
         * Previous to this time, a no-slip boundary condition is assumed.
         */
        double velocity_file_start_time;

        /**
         * Time at which the velocity file with number 0 shall be loaded.
         * Previous to this time, a no-slip boundary condition is assumed.
         */
        double first_velocity_file_time;

        /**
         * Directory in which the gplates velocity are present.
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
         * Determines the width of the velocity interpolation zone around the
         * current point. Currently equals the arc distance between evaluation
         * point and velocity data point that is still included in the
         * interpolation. The weighting of the points currently only accounts
         * for the surface area a single data point is covering ('moving
         * window' interpolation without distance weighting).
         */
        double interpolation_width;

        /**
         * Pointer to an object that reads and processes data we get from
         * gplates files.
         */
        std_cxx1x::shared_ptr<internal::BoxPlatesLookup<dim> > lookup;

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
