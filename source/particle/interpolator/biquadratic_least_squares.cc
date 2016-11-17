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

#include <aspect/particle/interpolator/biquadratic_least_squares.h>

#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the cell-wise averaged properties of all tracers of the cell containing the
       * given positions.
       */
      template <int dim>
      std::vector<std::vector<double> >
      BiquadraticLeastSquares<dim>::properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                                         const std::vector<Point<dim> > &positions,
                                                         const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        const types::LevelInd cell_index = std::make_pair<unsigned int, unsigned int> (cell->level(),cell->index());
        const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator,
              typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator> particle_range = particles.equal_range(cell_index);

        const unsigned int n_particles = std::distance(particle_range.first,particle_range.second);
        AssertThrow(n_particles != 0,
                    ExcMessage("At least one cell contained no particles. The 'constant "
                               "average' interpolation scheme does not support this case. "));
        const unsigned int n_properties = particles.begin()->second.get_properties().size();

        const unsigned int n_coefficients = 6;

        std::vector<std::vector<double> > properties(positions.size());

        dealii::FullMatrix<double> A(n_particles,n_coefficients);
        A = 0;

        unsigned int index = 0;
        for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator particle = particle_range.first;
             particle != particle_range.second; ++particle, ++index)
          {
            const Point<dim> location = particle->second.get_location();
            A(index,0) = 1;
            A(index,1) = location[0];
            A(index,2) = location[1];
            A(index,3) = location[0] * location[1];
            A(index,4) = location[0] * location[0];
            A(index,5) = location[1] * location[1];
          }

        dealii::FullMatrix<double> B(n_coefficients, n_coefficients);
        A.Tmmult(B, A, false);
        dealii::FullMatrix<double> B_inverse(B);
        B_inverse.gauss_jordan();

        index = 0;
        for (unsigned int i = 0; i < n_properties; ++i)
          {
            dealii::FullMatrix<double> r(6,1);
            r = 0;
            double max_value_for_particle_property=(particle_range.first)->second.get_properties()[i];
            double min_value_for_particle_property=(particle_range.first)->second.get_properties()[i];
            for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator particle = particle_range.first;
                 particle != particle_range.second; ++particle, ++index)
              {
                const double particle_property = particle->second.get_properties()[i];
                if (max_value_for_particle_property<particle_property)
                  max_value_for_particle_property=particle_property;
                if (min_value_for_particle_property>particle_property)
                  min_value_for_particle_property=particle_property;
                const Point<dim> position = particle->second.get_location();
                r(0,0) += particle_property;
                r(1,0) += particle_property * position[0];
                r(2,0) += particle_property * position[1];
                r(3,0) += particle_property * position[0] * position[1];
                r(4,0) += particle_property * position[0] * position[0];
                r(5,0) += particle_property * position[1] * position[1];
              }

            dealii::FullMatrix<double> c(6,1);
            c = 0;
            B_inverse.mmult(c, r);

            unsigned int index_positions = 0;
            for (typename std::vector<Point<dim> >::const_iterator itr = positions.begin(); itr != positions.end(); ++itr, ++index_positions)
              {
                Point<dim> support_point = *itr;
                double interpolated_value = c(0,0) + c(1,0)*(support_point[0]) + c(2,0)*(support_point[1]) + c(3,0)*(support_point[0] * support_point[1]) +  c(4,0)*(support_point[0] * support_point[0]) + c(5,0)*(support_point[1] * support_point[1]);
                if ( interpolated_value > max_value_for_particle_property)
                  interpolated_value  = max_value_for_particle_property;
                else if ( interpolated_value < min_value_for_particle_property)
                  interpolated_value  = min_value_for_particle_property;
                properties[index_positions].push_back(interpolated_value);
              }
          }
        return properties;
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(BiquadraticLeastSquares,
                                            "biquadratic",
                                            "Return the average of all tracer properties in the given cell.")
    }
  }
}
