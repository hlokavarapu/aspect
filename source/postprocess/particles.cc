/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
#include <aspect/postprocess/particles.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <stdio.h>
#include <unistd.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    Particles<dim>::Particles ()
      :

      last_output_time(std::numeric_limits<double>::quiet_NaN())
    {}

    template <int dim>
    Particles<dim>::~Particles ()
    {}

    template <int dim>
    void
    Particles<dim>::initialize ()
    {}

    template <int dim>
    void
    Particles<dim>::generate_particles()
    {
      world.generate_particles();
    }

    template <int dim>
    void
    Particles<dim>::initialize_particles()
    {
      world.initialize_particles();
    }

    template <int dim>
    const Particle::World<dim> &
    Particles<dim>::get_particle_world() const
    {
      return world;
    }

    template <int dim>
    Particle::World<dim> &
    Particles<dim>::get_particle_world()
    {
      return world;
    }

    template <int dim>
    std::pair<std::string,std::string>
    Particles<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
        }

      // Do not advect the particles in the initial refinement stage
      const bool in_initial_refinement = (this->get_timestep_number() == 0)
                                         && (this->get_pre_refinement_step() < this->get_parameters().initial_adaptive_refinement);
      if (!in_initial_refinement)
        // Advance the particles in the world to the current time
        world.advance_timestep();

      statistics.add_value("Number of advected particles",world.n_global_particles());

      // If it's not time to generate an output file or we do not write output
      // return early with the number of particles that were advected
      if (this->get_time() < last_output_time + output_interval)
        return std::make_pair("Number of advected particles:",
                              Utilities::int_to_string(world.n_global_particles()));

      if (world.get_property_manager().need_update() == Particle::Property::update_output_step)
        world.update_particles();

      set_last_output_time (this->get_time());

      const std::string data_file_names = world.generate_output();

      // If we do not write output return early with the number of particles
      // that were advected
      if (data_file_names.empty())
        return std::make_pair("Number of advected particles:",
                              Utilities::int_to_string(world.n_global_particles()));

      statistics.add_value ("Particle file name(s)", data_file_names);
      return std::make_pair("Writing particle output:", data_file_names);
    }



    template <int dim>
    void
    Particles<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }


    template <int dim>
    template <class Archive>
    void Particles<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &last_output_time
      ;
    }


    template <int dim>
    void
    Particles<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;
      aspect::oarchive oa (os);

      world.save(os);
      oa << (*this);

      status_strings["Particles"] = os.str();
    }


    template <int dim>
    void
    Particles<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("Particles") != status_strings.end())
        {
          std::istringstream is (status_strings.find("Particles")->second);
          aspect::iarchive ia (is);

          // Load the particle world
          world.load(is);

          ia >> (*this);
        }
    }


    template <int dim>
    void
    Particles<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          prm.declare_entry ("Time between data output", "1e8",
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "output files. A value of zero indicates that "
                             "output should be generated every time step.\n\n"
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      Particle::World<dim>::declare_parameters(prm);
    }


    template <int dim>
    void
    Particles<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          output_interval = prm.get_double ("Time between data output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          AssertThrow(this->get_parameters().run_postprocessors_on_nonlinear_iterations == false,
                      ExcMessage("Postprocessing nonlinear iterations in models with "
                                 "particles is currently not supported."));

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      // Initialize the particle world
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&world))
        sim->initialize_simulator (this->get_simulator());
      world.parse_parameters(prm);
      world.initialize();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Particles,
                                  "particles",
                                  "A Postprocessor that creates particles that follow the "
                                  "velocity field of the simulation. The particles can be generated "
                                  "and propagated in various ways and they can carry a number of "
                                  "constant or time-varying properties. The postprocessor can write "
                                  "output positions and properties of all particles at chosen intervals, "
                                  "although this is not mandatory. It also allows other parts of the "
                                  "code to query the particles for information.")
  }
}
