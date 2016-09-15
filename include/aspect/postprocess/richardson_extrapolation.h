/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#ifndef __aspect__postprocess_richardson_extrapolation_h
#define __aspect__postprocess_richardson_extrapolation_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that evaluates the solution vector at individual
     * points.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class RichardsonExtrapolation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        void initialize ();

        /**
         * Evaluate the solution and determine the values at the
         * selected points.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Declare the parameters this class takes through input files.
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

        bool
        read_in_data_for_h_over_2_resolution();

        bool
        save_data_for_h_over_2_resolution();

        /**
         *  Function that computes the error between interpolated (read in) values
         *  with the current solution at current nodal points.
         */
        void
        compute_error ();

        std::string
        get_formatted_file_name();

      private:

        BlockVector::

        /**
         * Run time parameters.
         */
        std::string file_name;
        double end_time;
    };
  }
}


#endif
