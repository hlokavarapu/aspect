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

#ifndef __aspect__particle_output_vtu_h
#define __aspect__particle_output_vtu_h

#include <aspect/particle/output/interface.h>


namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      template <int dim, class T>
       class VTUOutput : public Interface<dim, T>
       {
         public:
           /**
            * Constructor.
            *
            * @param[in] The directory into which output files shall be placed.
            * @param[in] The MPI communicator that describes this simulation.
            */
           VTUOutput();

           /**
            * Write data about the particles specified in the first argument
            * to a file. If possible, encode the current simulation time
            * into this file using the data provided in the second argument.
            *
            * @param[in] particles The set of particles to generate a graphical
            *   representation for
            * @param[in] current_time Current time of the simulation, given as either
            *   years or seconds, as selected in the input file. In other words,
            *   output writers do not need to know the units in which time is
            *   described.
            * @return The name of the file that was written, or any other
            *   information that describes what output was produced if for example
            *   multiple files were created.
            */
           virtual
           std::string
           output_particle_data(const std::multimap<LevelInd, T> &particles,
                                std::vector<MPIDataInfo> &data_info,
                                const double &current_time);

         private:

           /**
            * A list of pairs (time, pvtu_filename) that have so far been written
            * and that we will pass to DataOutInterface::write_pvd_record
            * to create a master file that can make the association
            * between simulation time and corresponding file name (this
            * is done because there is no way to store the simulation
            * time inside the .pvtu or .vtu files).
            */
 //TODO: This needs to be serialized
           std::vector<std::pair<double,std::string> > times_and_pvtu_file_names;

           /**
            * Like the previous variable, but for the .visit file.
            */
           std::vector<std::string>                    vtu_file_names;
       };
    }
  }
}

#endif
