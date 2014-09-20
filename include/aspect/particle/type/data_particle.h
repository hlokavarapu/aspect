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

#ifndef __aspect__particle_type_data_particle_h
#define __aspect__particle_type_data_particle_h

#include <aspect/particle/type/base_particle.h>

namespace aspect
{
  namespace Particle
  {
    namespace Type
    {
    /**
     * DataParticle provides an example of how to extend the BaseParticle
     * class to include related particle data. This allows users to attach
     * scalars/vectors/tensors/etc to particles and ensure they are
     * transmitted correctly over MPI and written to output files.
     */
    template <int dim, int data_dim>
    class DataParticle : public BaseParticle<dim>
    {
      private:
        double      val[data_dim];

      private:
        DataParticle();

        DataParticle(const Point<dim> &new_loc, const double &new_id);

        static unsigned int data_len();

        virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos);


        virtual void write_data(std::vector<double> &data) const;

        /**
         * Returns a vector from the first dim components of val
         *
         * @return vector representation of first dim components of the
         * particle data
         */
        Point<dim>
        get_vector () const;

        /**
         * Sets the first dim components of val to the specified vector value
         *
         * @param [in] new_vec Vector to set the DataParticle data to
         */
        void set_vector(Point<dim> new_vec);

        /**
         * Return a reference to an element of the DataParticle data
         *
         * @param [in] ind Index of the data array @return Reference to double
         * value at the requested index
         */
        double &operator[](const unsigned int &ind);

        /**
         * Return the value of an element of the DataParticle data
         *
         * @param [in] ind Index of the data array @return Value at the
         * requested index
         */
        double operator[](const unsigned int &ind) const;

        /**
         * Set up the MPI data type information for the DataParticle type
         *
         * @param [in,out] data_info Vector to append MPIDataInfo objects to
         */
        static void add_mpi_types(std::vector<MPIDataInfo> &data_info);
    };
    }
  }
}

#endif

