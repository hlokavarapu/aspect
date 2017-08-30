/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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
#include <aspect/particle/integrator/euler.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Euler scheme integrator, where $y_{n+1} = y_n + dt * v(y_n)$.
       * This requires only one step per integration, and doesn't involve any extra data.
       */
      template <int dim>
      void
      Euler<dim>::local_integrate_step(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                                       const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle,
                                       const std::vector<Tensor<1,dim> > &old_velocities,
                                       const std::vector<Tensor<1,dim> > &,
                                       const double dt)
      {
        Assert(static_cast<unsigned int> (std::distance(begin_particle, end_particle)) == old_velocities.size(),
               ExcMessage("The particle integrator expects the velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        typename std::vector<Tensor<1,dim> >::const_iterator old_velocity = old_velocities.begin();

        const QGauss<dim> quadrature_formula(this->get_parameters().stokes_velocity_degree + 1);
        FEValues<dim> fe_values(this->get_mapping(),
                                this->get_fe(),
                                quadrature_formula,
                                update_inverse_jacobians);

        fe_values.reinit (cell);

        std::vector<DerivativeForm<1,dim,dim> > inverse_jacobians = fe_values.get_inverse_jacobians();

        for (typename std::multimap<types::LevelInd, Particle<dim> >::iterator it = begin_particle;
             it != end_particle; ++it, ++old_velocity)
          {
            const Point<dim> loc = it->second.get_location();
            it->second.set_location(loc + dt * (*old_velocity));
          }
      }

//      template <int dim>
//      void
//      Euler<dim>::local_integrate_ref_loc(const typename DoFHandler<dim>::active_cell_iterator &cell,
//                                          const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
//                                          const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle,
//                                          const std::vector<Tensor<1,dim> > &old_velocities,
//                                          const std::vector<Tensor<1,dim> > &velocities,
//                                          const double dt)
//      {
//        Assert(static_cast<unsigned int> (std::distance(begin_particle, end_particle)) == old_velocities.size(),
//               ExcMessage("The particle integrator expects the velocity vector to be of equal size "
//                          "to the number of particles to advect. For some unknown reason they are different, "
//                          "most likely something went wrong in the calling function."));
//
//        typename std::vector<Tensor<1,dim> >::const_iterator old_velocity = old_velocities.begin();
//
////          if (update_flags & update_values)
////              fe_value.get_function_values (this->get_solution(),
////                                            values);
//
//        //it->second.set_reference_location(it->second.get_reference_location() + dt * (*old_velocity));
//      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(Euler,
                                          "euler",
                                          "Explicit Euler scheme integrator, where $y_{n+1} = y_n + dt * v(y_n)$. "
                                          "This requires only one integration substep per timestep.")
    }
  }
}
