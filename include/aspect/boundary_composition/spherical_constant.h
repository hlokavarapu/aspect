/*
  Copyright (C) 2013 - 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_boundary_composition_spherical_constant_h
#define _aspect_boundary_composition_spherical_constant_h

#include <aspect/boundary_composition/interface.h>

#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    /**
     * A class that implements a composition boundary condition for a
     * spherical shell geometry in which the composition at the inner and
     * outer surfaces (i.e. at the core-mantle and the mantle-
     * lithosphere/atmosphere boundaries) are constant.
     * This class works for the sphere, spherical shell,
     * chunk and ellipsoidal chunk geometries.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class SphericalConstant : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * This function returns the constant compositions read from the
         * parameter file for the inner and outer boundaries.
         *
         * @copydoc aspect::BoundaryComposition::Interface::boundary_composition()
         */
        virtual
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const;

        /**
         * Return the minimal composition on that part of the boundary on
         * which Dirichlet conditions are posed.
         */
        virtual
        double minimal_composition (const std::set<types::boundary_id> &fixed_boundary_ids) const;

        /**
         * Return the maximal composition on that part of the boundary on
         * which Dirichlet conditions are posed.
         */
        virtual
        double maximal_composition (const std::set<types::boundary_id> &fixed_boundary_ids) const;

        /**
         * Declare the parameters this class takes through input files. This
         * class declares the inner and outer boundary compositions.
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
         * Compositions at the inner and outer boundaries.
         */
        double inner_composition;
        double outer_composition;
    };
  }
}


#endif
