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

#include <aspect/particle/generator/uniform_box.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <deal.II/base/std_cxx11/array.h>
#include <aspect/geometry_model/box.h>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      UniformBox<dim>::generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles)
      {
        const Tensor<1,dim> P_diff = P_max - P_min;

        double volume(1.0);
        for (unsigned int i = 0; i < dim; ++i)
          volume *= P_diff[i];

        std_cxx11::array<unsigned int,dim> n_particles_per_direction;
        std_cxx11::array<double,dim> spacing;

        // Calculate separation of particles
        for (unsigned int i = 0; i < dim; ++i)
          {
            n_particles_per_direction[i] = round(std::pow(n_tracers * std::pow(P_diff[i],dim) / volume, 1./dim));
            spacing[i] = P_diff[i] / fmax(n_particles_per_direction[i] - 1,1);
          }

        typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active(),
                endc = this->get_dof_handler().end();

        const GeometryModel::Box<dim> *domain = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
        Point<dim> extent = domain->get_extents();

        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
          {
            if(dim == 2) {
              // Find the dimensions of the current active cell owned by the current processor.

              double local_x_min, local_y_min;
              double local_x_max, local_y_max;

              local_x_min = cell->vertex(0)[0], local_y_min = cell->vertex(0)[1];
              local_x_max = cell->vertex(0)[0], local_y_max = cell->vertex(0)[1];

              for (int i = 1; i < GeometryInfo<dim>::vertices_per_cell; i++) {
                if (cell->vertex(i)[0] < local_x_min)
                  local_x_min = cell->vertex(i)[0];
                if (cell->vertex(i)[1] < local_y_min)
                  local_y_min = cell->vertex(i)[1];
                if (cell->vertex(i)[0] > local_x_max)
                  local_x_max = cell->vertex(i)[0];
                if (cell->vertex(i)[1] > local_y_max)
                  local_y_max = cell->vertex(i)[1];
              }

              /**
               * If the current active cell is not within the dimensions of the uniform box, then
               * we continue to the next active cell in the local domain.
               */

              if ((local_x_max <= P_min[0]) || (local_x_min >= P_max[0]) ||
                  (local_y_max <= P_min[1]) || (local_y_min >= P_max[1]))
                continue;

              /**
               * For each dimension, compute the min and max indicies that bound all possible particles that can be
               * generated within the current cell's domain. These indices are collected in vectors of integers.
               */

              std::vector<int> x_indices;

              int max_index = (int) std::ceil((local_x_max - P_min[0]) / spacing[0]);
              if ((max_index > n_particles_per_direction[0]) ||
                  (P_min[0] + (n_particles_per_direction[0] - 1) * spacing[0] == local_x_max))
                max_index = n_particles_per_direction[0];

              for (int i = (int) std::ceil((local_x_min - P_min[0]) / spacing[0]); i < max_index; i++) {
                double p_x = P_min[0] + spacing[0] * i;
                if (local_x_min <= p_x && p_x <= local_x_max) {
                  x_indices.push_back(i);
                }
              }

              std::vector<int> y_indices;

              max_index = (int) std::ceil((local_y_max - P_min[1]) / spacing[1]);
              if ((max_index > n_particles_per_direction[1]) ||
                  (P_min[1] + (n_particles_per_direction[1] - 1) * spacing[1] == local_y_max))
                max_index = n_particles_per_direction[1];

              for (int i = (int) std::ceil((local_y_min - P_min[1]) / spacing[1]); i < max_index; i++) {
                double p_y = P_min[1] + spacing[1] * i;
                if (local_y_min <= p_y && p_y <= local_y_max) {
                  y_indices.push_back(i);
                }
              }

              /**
               * If the current active cell contains no particles that need to be generated, then
               * we continue to the next active cell in the local domain.
               */

              if (x_indices.size() == 0 || y_indices.size() == 0)
                continue;

              /**
               * Finally, we generate the particles contained within the current cell.
               * We also hard code the edge cases where particles must be generated on the boundaries of the local cell.
               */

              for (int i = 0; i < x_indices.size(); i++) {
                for (int j = 0; j < y_indices.size(); j++) {
                  Point<dim> position(P_min[0] + x_indices[i] * spacing[0], P_min[1] + y_indices[j] * spacing[1]);
                  const Particle<dim> particle(position,
                                               x_indices[i] * n_particles_per_direction[1] + y_indices[j]);
                  const types::LevelInd p_cell(cell->level(), cell->index());
                  std::cout << "position" << position << "id:" << particle.get_id() << std::endl;

                  if (local_x_min <= position[0] && position[0] < local_x_max) {
                    particles.insert(std::make_pair(p_cell, particle));
                  }
                  else if (local_y_min <= position[1] && position[1] < local_y_max) {
                    particles.insert(std::make_pair(p_cell, particle));
                  }
                  else if (position[0] == extent[0]) {
                    particles.insert(std::make_pair(p_cell, particle));
                  }
                  else if (position[1] == extent[1]) {
                    particles.insert(std::make_pair(p_cell, particle));
                  }
                }
              }
            }
            else if (dim == 3)
            {
              // Find the dimensions of the current active cell owned by the current processor.

                double local_x_min, local_y_min, local_z_min;
                double local_x_max, local_y_max, local_z_max;

                local_x_min = cell->vertex(0)[0], local_y_min = cell->vertex(0)[1], local_z_min = cell->vertex(0)[2];
                local_x_max = cell->vertex(0)[0], local_y_max = cell->vertex(0)[1], local_z_max = cell->vertex(0)[2];

                for (int i = 1; i < GeometryInfo<dim>::vertices_per_cell; i++) {
                  if (cell->vertex(i)[0] < local_x_min)
                    local_x_min = cell->vertex(i)[0];
                  if (cell->vertex(i)[1] < local_y_min)
                    local_y_min = cell->vertex(i)[1];
                  if (cell->vertex(i)[2] < local_z_min)
                    local_z_min = cell->vertex(i)[2];
                  if (cell->vertex(i)[0] > local_x_max)
                    local_x_max = cell->vertex(i)[0];
                  if (cell->vertex(i)[1] > local_y_max)
                    local_y_max = cell->vertex(i)[1];
                  if (cell->vertex(i)[2] > local_z_max)
                    local_z_max = cell->vertex(i)[2];
                }

                /**
                * If the current active cell is not within the dimensions of the uniform box, then
                * we continue to the next active cell in the local domain.
                */

                if ((local_x_max <= P_min[0]) || (local_x_min >= P_max[0]) ||
                    (local_y_max <= P_min[1]) || (local_y_min >= P_max[1]) ||
                    (local_z_max <= P_min[2]) || (local_z_min >= P_max[2]))
                  continue;

                /**
                 * For each dimension, compute the min and max indicies that bound all possible particles that can be
                 * generated within the current cell's domain. These indices are collected in vectors of integers.
                 */

                std::vector<int> x_indices;

                int max_index = (int) std::ceil((local_x_max - P_min[0]) / spacing[0]);
                if ((max_index > n_particles_per_direction[0]) ||
                    (P_min[0] + (n_particles_per_direction[0] - 1) * spacing[0] == local_x_max))
                  max_index = n_particles_per_direction[0];

                for (int i = (int) std::ceil((local_x_min - P_min[0]) / spacing[0]); i < max_index; i++) {
                  double p_x = P_min[0] + spacing[0] * i;
                  if (local_x_min <= p_x && p_x <= local_x_max) {
                    x_indices.push_back(i);
                  }
                }

                std::vector<int> y_indices;

                max_index = (int) std::ceil((local_y_max - P_min[1]) / spacing[1]);
                if ((max_index > n_particles_per_direction[1]) ||
                    (P_min[1] + (n_particles_per_direction[1] - 1) * spacing[1] == local_y_max))
                  max_index = n_particles_per_direction[1];

                for (int i = (int) std::ceil((local_y_min - P_min[1]) / spacing[1]); i < max_index; i++) {
                  double p_y = P_min[1] + spacing[1] * i;
                  if (local_y_min <= p_y && p_y <= local_y_max) {
                    y_indices.push_back(i);
                  }
                }

                std::vector<int> z_indices;

                max_index = (int) std::ceil((local_z_max - P_min[2]) / spacing[2]);
                if ((max_index > n_particles_per_direction[2]) ||
                    (P_min[2] + (n_particles_per_direction[2] - 1) * spacing[2] == local_z_max))
                  max_index = n_particles_per_direction[2];

                for (int i = (int) std::ceil((local_z_min - P_min[2]) / spacing[2]); i < max_index; i++) {
                  double p_z = P_min[2] + spacing[2] * i;
                  if (local_z_min <= p_z && p_z <= local_z_max) {
                    z_indices.push_back(i);
                  }
                }

                /**
                 * If the current active cell contains no particles that need to be generated, then
                 * we continue to the next active cell in the local domain.
                 */

                if (x_indices.size() == 0 || y_indices.size() == 0 || z_indices.size() == 0)
                  continue;

                /**
                 * Finally, we generate the particles contained within the current cell.
                 * We also hard code the edge cases where particles must be generated on the boundaries of the local cell.
                 */

                for (int i = 0; i < x_indices.size(); i++) {
                  for (int j = 0; j < y_indices.size(); j++) {
                    for (int k = 0; k < z_indices.size(); k++) {
                      Point<dim> position(P_min[0] + x_indices[i] * spacing[0], P_min[1] + y_indices[j] * spacing[1], P_min[2] + z_indices[k] * spacing[2]);
                      int particle_id = x_indices[i] * n_particles_per_direction[2] * n_particles_per_direction[1] + y_indices[j] * n_particles_per_direction[2] + z_indices[k];
                      const Particle<dim> particle(position, particle_id);
                      const types::LevelInd p_cell(cell->level(), cell->index());
                      if (local_x_min <= position[0] && position[0] < local_x_max &&
                              local_z_min <= position[2] && position[2] < local_z_max) {
                        particles.insert(std::make_pair(p_cell, particle));
                      }
                      else if (local_x_min <= position[0] && position[0] < local_x_max &&
                              local_y_min <= position[1] && position[1] < local_y_max) {
                        particles.insert(std::make_pair(p_cell, particle));
                      }
                      else if (local_y_min <= position[1] && position[1] < local_y_max &&
                              local_z_min <= position[2] && position[2] < local_z_max) {
                        particles.insert(std::make_pair(p_cell, particle));
                      }
                      else if (position[0] == extent[0]) {
                        particles.insert(std::make_pair(p_cell, particle));
                      }
                      else if (position[1] == extent[1]) {
                        particles.insert(std::make_pair(p_cell, particle));
                      }
                      else if (position[2] == extent[2]) {
                        particles.insert(std::make_pair(p_cell, particle));
                      }
                    }
                  }
                }

//                const Point<dim> particle_position = Point<dim> (P_min[0]+i*spacing[0],P_min[1]+j*spacing[1],P_min[2]+k*spacing[2]);
//
//                // Try to add the particle. If it is not in this domain, do not
//                // worry about it and move on to next point.
//                try
//                  {
//                    particles.insert(this->generate_particle(particle_position,particle_index));
//                  }
//                catch (ExcParticlePointNotInDomain &)
//                  {}
//                particle_index++;
              }
            else
              ExcNotImplemented();
          }


//        types::particle_index particle_index = 0;
//
//        for (unsigned int i = 0; i < n_particles_per_direction[0]; ++i)
//          {
//            for (unsigned int j = 0; j < n_particles_per_direction[1]; ++j)
//              {
//                if (dim == 2)
//                  {
//                    const Point<dim> particle_position = Point<dim> (P_min[0]+i*spacing[0],P_min[1]+j*spacing[1]);
//
//                    // Try to add the particle. If it is not in this domain, do not
//                    // worry about it and move on to next point.
//                    try
//                      {
//                        particles.insert(this->generate_particle(particle_position,particle_index));
//                      }
//                    catch (ExcParticlePointNotInDomain &)
//                      {}
//                    particle_index++;
//                  }
//                else if (dim == 3)
//                  for (unsigned int k = 0; k < n_particles_per_direction[2]; ++k)
//                    {
//                      const Point<dim> particle_position = Point<dim> (P_min[0]+i*spacing[0],P_min[1]+j*spacing[1],P_min[2]+k*spacing[2]);
//
//                      // Try to add the particle. If it is not in this domain, do not
//                      // worry about it and move on to next point.
//                      try
//                        {
//                          particles.insert(this->generate_particle(particle_position,particle_index));
//                        }
//                      catch (ExcParticlePointNotInDomain &)
//                        {}
//                      particle_index++;
//                    }
//                else
//                  ExcNotImplemented();
//              }
//          }
      }


      template <int dim>
      void
      UniformBox<dim>::declare_parameters (ParameterHandler &prm)
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
              prm.enter_subsection("Uniform box");
              {
                prm.declare_entry ("Minimum x", "0",
                                   Patterns::Double (),
                                   "Minimum x coordinate for the region of tracers.");
                prm.declare_entry ("Maximum x", "1",
                                   Patterns::Double (),
                                   "Maximum x coordinate for the region of tracers.");
                prm.declare_entry ("Minimum y", "0",
                                   Patterns::Double (),
                                   "Minimum y coordinate for the region of tracers.");
                prm.declare_entry ("Maximum y", "1",
                                   Patterns::Double (),
                                   "Maximum y coordinate for the region of tracers.");
                prm.declare_entry ("Minimum z", "0",
                                   Patterns::Double (),
                                   "Minimum z coordinate for the region of tracers.");
                prm.declare_entry ("Maximum z", "1",
                                   Patterns::Double (),
                                   "Maximum z coordinate for the region of tracers.");
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
      UniformBox<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            n_tracers    = static_cast<types::particle_index>(prm.get_double ("Number of tracers"));

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Uniform box");
              {
                P_min(0) = prm.get_double ("Minimum x");
                P_max(0) = prm.get_double ("Maximum x");
                P_min(1) = prm.get_double ("Minimum y");
                P_max(1) = prm.get_double ("Maximum y");

                if (dim == 3)
                  {
                    P_min(2) = prm.get_double ("Minimum z");
                    P_max(2) = prm.get_double ("Maximum z");
                  }
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
      ASPECT_REGISTER_PARTICLE_GENERATOR(UniformBox,
                                         "uniform box",
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
