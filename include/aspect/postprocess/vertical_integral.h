/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#ifndef __aspect__postprocess_vertical_integral_h
#define __aspect__postprocess_vertical_integral_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A class derived from CellDataVectorCreator that takes an output
     * vector and computes a variable that represents the vertical integral
     * of a given compositional field. This quantity only makes sense at the
     * surface of the domain. Thus, the value is set to zero in all the
     * cells inside of the domain.
     *
     * The member functions are all implementations of those declared in the
     * base class. See there for their meaning.
     */
    template <int dim>
    class VerticalIntegral
      : public Interface<dim>,
        public SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for the vertical integral of a given
         * compositional field.
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

      private:
        /**
         * A parameter that we read from the input file that denotes, which
         * compositional field to integrate.
         */
        std::string name_of_compositional_field;

        /**
         * A parameter that we read from the input file that denotes, at what
         * time this postprocessor will be executed. The default value of 0.0
         * is interpreted to produce output every timestep.
         */
        double time_of_output;

        /**
         * The format in which to produce graphical output. This also
         * determines the extension of the file name to which to write.
         */
        DataOutBase::OutputFormat output_format;

        /**
         * A parameter that can be used to exclude the upper part
         * of the model from integration. All cells with a smaller
         * depth are ignored.
         */
        double minimum_depth;

        /**
         * A parameter that can be used to exclude the lower part
         * of the model from integration. All cells with a larger
         * depth are ignored. The default value of 0 will be
         * replaced by the maximum depth of the model.
         */
        double maximum_depth;
    };
  }
}

#endif
