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


#ifndef __aspect__velocity_boundary_conditions_plume_h
#define __aspect__velocity_boundary_conditions_plume_h

#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/boundary_temperature/plume.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    using namespace dealii;

    /**
     * A class that implements a velocity boundary condition for a box
     * geometry with plume inflow.
     *
     * @ingroup VelocityBoundaryConditionsModels
     */
    template <int dim>
    class Plume : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Loads file.
         */
        virtual
        void
        initialize ();

        /**
         * A function that is called at the beginning of each time step. For
         * the current plugin, this function updates the current plume position.
         */
        virtual
        void
        update ();

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function returns a zero tensor except from the
         * region around the current plume position, where it prescribes a
         * vertical inflow of a certain magnitude with gaussian distribution.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files. This
         * class declares parameters concerning the location of the input file
         * and the strength and radius of the inflow.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Current plume position. Used to avoid looking up the plume position
         * for every quadrature point, since this is only time-dependent.
         */
        Point<dim> plume_position;

        /**
         * Pointer to an object that reads and processes data we get from
         * gplates files.
         */
        std_cxx1x::shared_ptr<BoundaryTemperature::internal::PlumeLookup<dim> > lookup;

        /**
         * Directory in which the plume file is present.
         */
        std::string data_directory;

        /**
         * Filename of plume file.
         */
        std::string plume_file_name;

        /**
         * Magnitude of the temperature anomaly
         */
        double V0;

        /**
         * Radius of the temperature anomaly
         */
        double R0;

    };
  }
}


#endif
