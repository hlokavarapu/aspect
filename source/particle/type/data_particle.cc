/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/particle/type/data_particle.h>

namespace aspect
{
  namespace Particle
  {
    namespace Type
    {
      template <int dim, int data_dim>
      DataParticle<dim,data_dim>::DataParticle()
      {
        for (unsigned int i=0; i<data_dim; ++i) val[i] = 0;
      }

      template <int dim, int data_dim>
      DataParticle<dim,data_dim>::DataParticle(const Point<dim> &new_loc, const double &new_id) : BaseParticle<dim>(new_loc, new_id)
          {
            for (unsigned int i=0; i<data_dim; ++i) val[i] = 0;
          }

      template <int dim, int data_dim>
          unsigned int
          DataParticle<dim,data_dim>::data_len()
          {
            return BaseParticle<dim>::data_len() + data_dim;
          }

      template <int dim, int data_dim>
          unsigned int
          DataParticle<dim,data_dim>::read_data(const std::vector<double> &data, const unsigned int &pos)
          {
            unsigned int  p;

            // Read the parent data first
            p = BaseParticle<dim>::read_data(data, pos);

            // Then read the DataParticle data
            for (unsigned i=0; i<data_dim; ++i)
              {
                val[i] = data.at(p++);
              }

            return p;
          }


      template <int dim, int data_dim>
          void
          DataParticle<dim,data_dim>::write_data(std::vector<double> &data) const
          {
            // Write the parent data first
            BaseParticle<dim>::write_data(data);

            // Then write the DataParticle data
            for (unsigned i=0; i<data_dim; ++i)
              {
                data.push_back(val[i]);
              }
          }

          /**
           * Sets the first dim components of val to the specified vector value
           *
           * @param [in] new_vec Vector to set the DataParticle data to
           */
      template <int dim, int data_dim>
          void
          DataParticle<dim,data_dim>::set_vector(Point<dim> new_vec)
          {
            AssertThrow(data_dim>=dim, std::out_of_range("set_vector"));
            for (unsigned int i=0; i<dim; ++i)
              val[i] = new_vec(i);
          }

          /**
           * Return a reference to an element of the DataParticle data
           *
           * @param [in] ind Index of the data array @return Reference to double
           * value at the requested index
           */
          template <int dim, int data_dim>
          double &
          DataParticle<dim,data_dim>::operator[](const unsigned int &ind)
          {
            AssertThrow(data_dim>ind, std::out_of_range("DataParticle[]"));
            return val[ind];
          }

          /**
           * Return the value of an element of the DataParticle data
           *
           * @param [in] ind Index of the data array @return Value at the
           * requested index
           */
          template <int dim, int data_dim>
          double
          DataParticle<dim,data_dim>::operator[](const unsigned int &ind) const
          {
            AssertThrow(data_dim>ind, std::out_of_range("DataParticle[]"));
            return val[ind];
          }

          /**
           * Set up the MPI data type information for the DataParticle type
           *
           * @param [in,out] data_info Vector to append MPIDataInfo objects to
           */
          template <int dim, int data_dim>
          void
          DataParticle<dim,data_dim>::add_mpi_types(std::vector<MPIDataInfo> &data_info)
          {
            // Set up the parent types first
            BaseParticle<dim>::add_mpi_types(data_info);

            // Then add our own
            data_info.push_back(aspect::Particle::MPIDataInfo("data", data_dim));
          }

      // A particle with associated values, such as scalars, vectors or tensors
      template <int dim, int data_dim>
      inline Point<dim>
      DataParticle<dim,data_dim>::get_vector () const
      {
        AssertThrow(data_dim >= dim, std::out_of_range ("get_vector"));
        Point < dim > p;
        for (unsigned int i = 0; i < dim; ++i)
          p (i) = val[i];
      }

    }
  }
}
