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
        typename parallel::distributed::Triangulation<dim>::active_cell_iterator found_cell;

        if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
          {
            // We can not simply use one of the points as input for find_active_cell_around_point
            // because for vertices of mesh cells we might end up getting ghost_cells as return value
            // instead of the local active cell. So make sure we are well in the inside of a cell.
            Point<dim> approximated_cell_midpoint = positions[0];
            if (positions.size() > 1)
              {
                Tensor<1,dim> direction_to_center;
                for (unsigned int i = 1; i<positions.size()-1; ++i)
                  direction_to_center += positions[i] - positions[0];
                direction_to_center /= positions.size() - 1;
                approximated_cell_midpoint += direction_to_center;
              }

            found_cell =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                           this->get_triangulation(),
                                                           approximated_cell_midpoint)).first;
          }
        else
          found_cell = cell;

        const types::LevelInd cell_index = std::make_pair<unsigned int, unsigned int> (found_cell->level(),found_cell->index());
        const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator,
              typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator> particle_range = particles.equal_range(cell_index);

        const unsigned int n_particles = std::distance(particle_range.first,particle_range.second);
        const unsigned int n_properties = particles.begin()->second.get_properties().size();

        std::vector<std::vector<double> > properties;
        std::vector<double> cell_properties (n_properties,0.0);
        AssertThrow(n_particles != 0,
                    ExcMessage("At least one cell contained no particles. The 'constant "
                               "average' interpolation scheme does not support this case. "));

//        const std::vector<double> particle_properties = particle->second.get_properties();
        for (typename std::vector<Point<dim>>::const_iterator itr = positions.begin(); itr != positions.end(); itr++)
          {
            std::vector<double> cell_properties2 (n_properties,0.0);
            for (unsigned int i = 0; i < n_properties; ++i)
              {
                dealii::FullMatrix<double> A(6,6); // = dealii::FullMatrix<double>(3,3);
                dealii::FullMatrix<double> r(6,1); // = dealii::FullMatrix<double>(3,3);
                A = 0;
                r = 0;

                for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator particle = particle_range.first;
                     particle != particle_range.second; ++particle)
                  {

                    const double particle_property = particle->second.get_properties()[i];
                    const Point<dim> position = particle->second.get_location();
                    A(0,0) += 1;
                    A(0,1) += position[0];
                    A(0,2) += position[1];
                    A(0,3) += position[0] * position[1];
                    A(0,4) += position[0] * position[0];
                    A(0,5) += position[1] * position[1];
                    A(1,1) += position[0] * position[0];
                    A(1,2) += position[0] * position[1];
                    A(1,3) += position[0] * position[0] * position[1];
                    A(1,4) += position[0] * position[0] * position[0];
                    A(1,5) += position[0] * position[1] * position[1];
                    A(2,2) += position[1] * position[1];
                    A(2,3) += position[0] * position[1] * position[1];
                    A(2,4) += position[0] * position[0] * position[1];
                    A(2,5) += position[1] * position[1] * position[1];
                    A(3,3) += position[0] * position[0] * position[1] * position[1];
                    A(3,4) += position[0] * position[0] * position[0] * position[1];
                    A(3,5) += position[0] * position[1] * position[1] * position[1];
                    A(4,4) += position[0] * position[0] * position[0] * position[0];
                    A(4,5) += position[0] * position[0] * position[1] * position[1];
                    A(5,5) += position[1] * position[1] * position[1] * position[1];


                    r(0,0) += particle_property;
                    r(1,0) += particle_property * position[0];
                    r(2,0) += particle_property * position[1];
                    r(3,0) += particle_property * position[0] * position[1];
                    r(4,0) += particle_property * position[0] * position[0];
                    r(5,0) += particle_property * position[1] * position[1];
                  }

                A(1,0) = A(0,1);
                A(2,0) = A(0,2);
                A(2,1) = A(1,2);
                A(3,0) = A(0,3);
                A(3,1) = A(1,3);
                A(3,2) = A(2,3);
                A(4,0) = A(0,4);
                A(4,1) = A(1,4);
                A(4,2) = A(2,4);
                A(4,3) = A(3,4);
                A(5,0) = A(0,5);
                A(5,1) = A(1,5);
                A(5,2) = A(2,5);
                A(5,3) = A(3,5);
                A(5,4) = A(4,5);

                dealii::FullMatrix<double> c(6,1);
                c = 0;
                dealii::FullMatrix<double> A_inverse(A);
                A_inverse.gauss_jordan();
                A_inverse.mmult(c, r);

                Point<dim> support_point = *itr;
                cell_properties2[i] = c(0,0) + c(1,0)*(support_point[0]) + c(2,0)*(support_point[1]) + c(3,0)*(support_point[0] * support_point[1]) +  c(4,0)*(support_point[0] * support_point[0]) + c(5,0)*(support_point[1] * support_point[1]);
              }
            properties.push_back(cell_properties2);
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
