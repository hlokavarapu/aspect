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


#ifndef __aspect__velocity_boundary_conditions_plume_box_h
#define __aspect__velocity_boundary_conditions_plume_box_h

#include <aspect/velocity_boundary_conditions/interface.h>

// Additional lookup classes are within these
#include <aspect/velocity_boundary_conditions/ascii_data.h>
#include <aspect/velocity_boundary_conditions/box_plates.h>
#include <aspect/boundary_temperature/plume.h>


#include <aspect/simulator_access.h>

namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    using namespace dealii;

    /**
     * A class that implements prescribed velocity boundary conditions
     * determined from a AsciiData input file.
     *
     * @ingroup VelocityBoundaryConditionsModels
     */
    template <int dim>
    class PlumeBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        PlumeBox ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        virtual
        void
        initialize ();

        /**
         * A function that is called at the beginning of each time step. For
         * the current plugin, this function loads the next data files if
         * necessary and outputs a warning if the end of the set of data
         * files is reached.
         */
        virtual
        void
        update ();

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function returns value from the text files.
         */
        Tensor<1,dim>
        boundary_velocity (const Point<dim> &position) const;


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
         * A variable that stores the currently used data file of a
         * series. It gets updated if necessary by set_current_time.
         */
        int  current_file_number;

        /**
         * Time from which on the data file with number 'First data
         * file number' is used as boundary condition. Previous to this
         * time, a no-slip boundary condition is assumed. Depending on the setting
         * of the global 'Use years in output instead of seconds' flag
         * in the input file, this number is either interpreted as seconds or as years."
         */
        double first_data_file_model_time;

        /**
         * Number of the first data file to be loaded when the model time
         * is larger than 'First data file model time'.
         */
        int first_data_file_number;

        /**
         * In some cases the boundary files are not numbered in increasing
         * but in decreasing order (e.g. 'Ma BP'). If this flag is set to
         * 'True' the plugin will first load the file with the number
         * 'First data file number' and decrease the file number during
         * the model run.
         */
        bool decreasing_file_order;

        /**
         * Directory in which the data files are present.
         */
        std::string data_directory;

        /**
         * Filename of data file. The file names can contain
         * the specifiers %s and/or %c (in this order), meaning the name of the
         * boundary and the number of the data file time step.
         */
        std::string data_file_name;

        /**
         * Directory in which the surface velocity files are present.
         */
        std::string surface_data_directory;

        /**
         * Filename of velocity file.
         */
        std::string surface_velocity_file_name;

        /**
         * First part of filename of id files. The files have to have
         * the pattern velocity_file_name.n.gpml where n is the number of the
         * current timestep (starts from 0).
         */
        std::string surface_id_file_names;

        /**
         * Distance between two surface grid points in x direction.
         */
        double x_step;

        /**
         * Distance between two surface grid points in y direction.
         */
        double y_step;

        /**
         * Boundary id for which this plugin is created. Is initialized in
         * initialize() and used to determine the appropriate file name.
         */
        types::boundary_id boundary_id;

        /**
         * Time in model units (depends on other model inputs) between two
         * data files.
         */
        double data_file_time_step;

        /**
         * Number of grid points in data file
         */
        std_cxx1x::array<unsigned int,3> data_points;

        /**
         * Weight between data file n and n+1 while the current time is
         * between the two values t(n) and t(n+1).
         */
        double time_weight;

        /**
         * State whether we have time_dependent boundary conditions. Switched
         * off after finding no more data files to suppress attempts to
         * read in new files.
         */
        bool time_dependent;

        /**
         * Scale the boundary condition by a scalar factor. Can be
         * used to transform the unit of the velocities (if they are not
         * specified in the default unit (m/s or m/yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double scale_factor;

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
         * Scale the time steps of the surface by a scalar factor. Can be
         * used to transform the unit of the velocities (if they are not
         * specified in the default unit (m/s or m/yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double surface_time_scale_factor;

        /**
         * Scale the surface velocity by a scalar factor. Can be
         * used to transform the unit of the velocities (if they are not
         * specified in the default unit (m/s or m/yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double surface_scale_factor;

        /**
         * Pointer to an object that reads and processes data we get from
         * text files.
         */
        std_cxx11::shared_ptr<internal::AsciiDataLookup<dim,dim-1> > lookup;

        /**
         * Pointer to an object that reads and processes data we get from
         * text files.
         */
        std_cxx11::shared_ptr<internal::BoxPlatesLookup<dim> > surface_lookup;

        /**
         * Pointer to an object that reads and processes data we get from
         * gplates files.
         */
        std_cxx1x::shared_ptr<BoundaryTemperature::internal::PlumeLookup<dim> > plume_lookup;

        /**
         * Current plume position. Used to avoid looking up the plume position
         * for every quadrature point, since this is only time-dependent.
         */
        Point<dim> plume_position;

        /**
         * Filename of plume file.
         */
        std::string plume_file_name;

        /**
         * Directory in which the plume file is present.
         */
        std::string plume_data_directory;

        /**
         * Magnitude of the plume velocity anomaly
         */
        double V0;

        /**
         * Radius of the plume velocity anomaly
         */
        double R0;

        /**
         * Handles the update of the data in lookup.
         */
        void
        update_data ();

        /**
         * Handles settings and user notification in case the time-dependent
         * part of the boundary condition is over.
         */
        void
        end_time_dependence (const int timestep);

        /**
         * Create a filename out of the name template for the side and bottom
         * ascii data model.
         */
        std::string
        create_filename (const int filenumber) const;

        /**
         * Create a filename out of the name template for the surface
         * box plates model.
         */
        std::string
        create_surface_filename (const int filenumber) const;

    };
  }
}


#endif
