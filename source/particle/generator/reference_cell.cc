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

#include <aspect/particle/generator/reference_cell.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <deal.II/base/std_cxx11/array.h>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      ReferenceCell<dim>::generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles)
      {
        std_cxx11::array<unsigned int,dim> n_particles_per_direction;
        std_cxx11::array<double,dim> spacing;

        // Calculate separation of particles
        for (unsigned int i = 0; i < dim; ++i)
          {
            n_particles_per_direction[i] = n_particles_per_direction_per_cell;
            spacing[i] = 1.0/n_particles_per_direction[i];
          }


        types::particle_index particle_index = 0;

        typename Triangulation<dim>::active_cell_iterator
        cell = this->get_triangulation().begin_active(),
        endc = this->get_triangulation().end();
        for (; cell != endc; cell++)
          {
            for (unsigned int i = 0; i < n_particles_per_direction[0]; ++i)
              {
                for (unsigned int j = 0; j < n_particles_per_direction[1]; ++j)
                  {
                    if (dim == 2)
                      {
                        Point<dim> position_unit = Point<dim> (i*spacing[0] + spacing[0]/2,j*spacing[1] + spacing[1]/2);
                        Point<dim> position_real = this->get_mapping().transform_unit_to_real_cell(cell, position_unit);
                        // Try to add the particle. If it is not in this domain, do not
                        // worry about it and move on to next point.
                        try
                          {
                            particles.insert(this->generate_particle(cell, position_unit, position_real, particle_index));
                          }
                        catch (ExcParticlePointNotInDomain &)
                          {}
                        particle_index++;
                      }
                    else if (dim == 3)
                      for (unsigned int k = 0; k < n_particles_per_direction[2]; ++k)
                        {
                          Point<dim> position_unit = Point<dim> (i*spacing[0] + spacing[0]/2,j*spacing[1] + spacing[1]/2, k*spacing[2] + spacing[2]/2);
                          Point<dim> position_real = this->get_mapping().transform_unit_to_real_cell(cell, position_unit);
                          // Try to add the particle. If it is not in this domain, do not
                          // worry about it and move on to next point.
                          try
                            {
                              particles.insert(this->generate_particle(cell, position_unit, position_real, particle_index));
                            }
                          catch (ExcParticlePointNotInDomain &)
                            {}
                          particle_index++;
                        }
                    else
                      ExcNotImplemented();
                  }
              }
          }
      }


      template <int dim>
      void
      ReferenceCell<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            prm.declare_entry ("Number of tracers", "1000",
                               Patterns::Double (0),
                               "Total number of tracers to create (not per processor or per element). "
                               "The number is parsed as a floating point number (so that one can "
                               "specify, for example, '1e4' particles) but it is interpreted as "
                               "an integer, of course.");

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Reference cell");
              {
                prm.declare_entry ("Number of tracers per cell per direction", "10",
                                   Patterns::Double (0),
                                   "Total number of tracers to create per element per direction. "
                                   "The number is parsed as a floating point number (so that one can "
                                   "specify, for example, '1e4' particles) but it is interpreted as "
                                   "an integer, of course.");
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      ReferenceCell<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Reference cell");
              {
                n_particles_per_direction_per_cell = prm.get_double ("Number of tracers per cell per direction");
                n_tracers = static_cast<types::particle_index>(std::pow(n_particles_per_direction_per_cell,dim));
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      ASPECT_REGISTER_PARTICLE_GENERATOR(ReferenceCell,
                                         "reference cell",
                                         "Generate a uniform distribution of particles "
                                         "over a rectangular domain in 2D or 3D. Uniform here means "
                                         "the particles will be generated with an equal spacing in "
                                         "each spatial dimension. Note that in order "
                                         "to produce a regular distribution the number of generated "
                                         "tracers might not exactly match the one specified in the "
                                         "input file.")
    }
  }
}
