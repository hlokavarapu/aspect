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

#include <aspect/particle/generator/quadrature_points.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      QuadraturePoints<dim>::generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles)
      {
        /**
        * TODO: Need to setup an independent method that deals with globally unique particle ids independently from the
        * different particle generators.
        */
        types::particle_index particle_index = 0;

        const QGauss<dim> quadrature_formula(this->get_parameters().stokes_velocity_degree + 1);

        FEValues<dim> fe_values(this->get_mapping(),
                                this->get_fe(),
                                quadrature_formula,
                                update_values |
                                update_quadrature_points |
                                update_JxW_values);

        typename Triangulation<dim>::active_cell_iterator
        cell = this->get_triangulation().begin_active(),
        endc = this->get_triangulation().end();
        for (; cell != endc; cell++)
          {
            fe_values.reinit(cell);
            std::vector<Point<dim>> quadrature_points = fe_values.get_quadrature_points();
            for (typename std::vector<Point<dim>>::const_iterator q_points_itr = quadrature_points.begin();
                 q_points_itr != quadrature_points.end(); q_points_itr++)
              {
                Point<dim> particle_position_real = (*q_points_itr);
                Point<dim> particle_position_unit = this->get_mapping().transform_real_to_unit_cell(cell, (*q_points_itr));

                try
                  {
                    particles.insert(this->generate_particle(cell, particle_position_unit, particle_position_real, particle_index));
                  }
                catch (ExcParticlePointNotInDomain &)
                  {}
                particle_index++;
              }
          }
        n_tracers = static_cast<types::particle_index> (fe_values.n_quadrature_points*this->get_triangulation().n_active_cells());
      }


      template <int dim>
      void
      QuadraturePoints<dim>::declare_parameters (ParameterHandler &prm)
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
      QuadraturePoints<dim>::parse_parameters (ParameterHandler &prm)
      {
        /**
         * TODO:
         */
//          n_tracers = static_cast<types::particle_index>(this->get_triangulation().n_active_cells() * this->get)
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
      ASPECT_REGISTER_PARTICLE_GENERATOR(QuadraturePoints,
                                         "quadrature points",
                                         "Generate particles at the quadrature points of each active"
                                         "cell of the triangulation mesh. The number of tracers is set at parsing.")
    }
  }
}
