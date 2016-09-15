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


#include <aspect/postprocess/richardson_extrapolation.h>
#include <aspect/global.h>
#include <deal.II/numerics/vector_tools.h>

#include <math.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    RichardsonExtrapolation<dim>::execute (TableHandler &)
    {
      if ( this->get_time() == end_time)
      {
          Triangulation<dim> refined_mesh;
          refined_mesh.copy_triangulation(this->get_triangulation());
          refined_mesh.refine_global(2);
          MappingQ1<dim> refined_mapping;
          DoFHandler<dim> refined_dof_handler;
          FiniteElement<dim> refined_fe(this->get_fe().base_element(introspection<dim>().base_elements.pressure));

          refined_dof_handler.initialize(refined_mesh, refined_fe);

          refined_dof_handler.distribute_dofs(this->get_fe());

          LinearAlgebra::Vector interpolated_solution(refined_dof_handler.n_dofs());
          intro

          VectorTools::interpolate_to_different_mesh(this->get_dof_handler(), this->get_solution().block(introspection<dim>().block_indices.pressure), refined_dof_handler, interpolated_solution);

          std::fstream interpolated_data;
          interpolated_data.open("interpolated_values.dat");
          interpolated_solution.print(interpolated_data, 14);
          interpolated_data.close();
      }
        return std::pair<std::string, std::string> ("Richardson Extrapolation: ",
                                                    "");
    }


    template <int dim>
    void
    RichardsonExtrapolation<dim>::declare_parameters (ParameterHandler &prm)
    {
        prm.enter_subsection("Postprocess");
        {
            prm.enter_subsection("Richardson extrapolation");
            {
                prm.declare_entry ("End time",
                        /* boost::lexical_cast<std::string>(std::numeric_limits<double>::max() /
                                             year_in_seconds) = */ "5.69e+300",
                                   Patterns::Double (),
                                   "The end time of the simulation. The default value is a number "
                                           "so that when converted from years to seconds it is approximately "
                                           "equal to the largest number representable in floating point "
                                           "arithmetic. For all practical purposes, this equals infinity. "
                                           "Units: Years if the "
                                           "'Use years in output instead of seconds' parameter is set; "
                                           "seconds otherwise.");
            }
            prm.leave_subsection();
        }
        prm.leave_subsection();
    }

//      template <int dim>
//      template <class Archive>
//      void RichardsonExtrapolation<dim>::serialize (Archive &ar, const unsigned int)
//      {
////          ar &end_time;
//      }
//
//
//      template <int dim>
//      void
//      RichardsonExtrapolation<dim>::save (std::map<std::string, std::string> &status_strings) const
//      {
////          std::ostringstream os;
////          aspect::oarchive oa (os);
////          oa << (*this);
//      }
//
//
//      template <int dim>
//      void
//      RichardsonExtrapolation<dim>::load (const std::map<std::string, std::string> &status_strings)
//      {
////              std::istringstream is (status_strings.find("PointValues")->second);
////              aspect::iarchive ia (is);
////              ia >> (*this);
//      }

    template <int dim>
    void
    RichardsonExtrapolation<dim>::parse_parameters (ParameterHandler &prm)
    {
        bool use_years = prm.get_bool ("Use years in output instead of seconds");
        prm.enter_subsection("Postprocess");
        {
            prm.enter_subsection("Richardson extrapolation");
            {
                // read end time from parameter file. if it is to be interpreted
                // in years rather than seconds, then do the conversion
                end_time = prm.get_double ("End time");
                if ( use_years == true)
                    end_time *= year_in_seconds;
            }
            prm.leave_subsection();
        }
        prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(RichardsonExtrapolation,
                                  "richardson extrapolation",
                                  "A postprocessor called Richardson Extrapolation.")
  }
}
