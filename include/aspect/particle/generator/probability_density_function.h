/*
 Copyright (C) 2015 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_generator_probability_density_function_h
#define __aspect__particle_generator_probability_density_function_h

#include <aspect/particle/generator/interface.h>

#include <deal.II/base/parsed_function.h>

#include <boost/random.hpp>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generates a random distribution of particles over the simulation
       * domain. The particle density is determined by a user-defined
       * probability density function in the parameter file. This is done using
       * a "roulette wheel" style selection. Every cell is weighted by the value
       * of the provided function at its center multiplied with the cell
       * volume. Then a map between the accumulated cell weight
       * and the cell index of the current cell is constructed. Consequently,
       * a random number between zero and the global integral of the
       * probability density function uniquely defines one particular cell.
       * Afterwards, every process generates n_global_particles random numbers,
       * but only generates a tracer if it is the owner of the active cell
       * that is associated with this random number.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class ProbabilityDensityFunction : public Interface<dim>
      {
        public:
          /**
           * Constructor.
           *
           */
          ProbabilityDensityFunction();

          /**
           * Generate a set of particles in the current
           * particle world. The particle density is set by an analytically
           * prescribed density function that is set as an input parameter.
           * This function builds a list of probabilities for all local cells
           * and then calls generate_particles_in_subdomain() to generate
           * the local particles.
           *
           * @return A multimap containing cells and their contained particles.
           */
          virtual
          std::multimap<types::LevelInd, Particle<dim> >
          generate_particles();

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
           * Random number generator and an object that describes a
           * uniform distribution on the interval [0,1]. These
           * will be used to generate particle locations at random.
           */
          boost::mt19937            random_number_generator;
          boost::uniform_01<double> uniform_distribution_01;

          /**
           * Number of particles to create
           */
          types::particle_index n_tracers;

          /**
           * A function object representing the particle location probability
           * density.
           */
          Functions::ParsedFunction<dim> function;

          /**
           * Generate a set of particles distributed within the local domain
           * according to the probability density function.
           *
           * @param [in] cells Map between accumulated cell weight and cell index
           * @param [in] global_weight The integrated probability function
           * @param [in] start_weight The starting weight of the first cell of the local process.
           * @param [in] num_particles The total number of particles to generate.
           * @param [in] start_id The starting ID to assign to generated particles of the local process.
           * @return The particle world the particles will exist in.
           *
           */
          std::multimap<types::LevelInd, Particle<dim> >
          generate_particles_in_subdomain (const std::map<double,types::LevelInd> &cells,
                                           const double global_weight,
                                           const double start_weight,
                                           const types::particle_index num_particles,
                                           const types::particle_index start_id);

          /**
           * This function loops over all active cells in the local subdomain
           * and returns a vector of accumulated cell weights. This is done
           * by adding the return value of function multiplied with the cell
           * volume for every cell it loops over.
           */
          std::vector<double>
          local_accumulated_cell_weights () const;
      };

    }
  }
}

#endif
