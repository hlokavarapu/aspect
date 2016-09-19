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
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <aspect/global.h>
//#include <deal.II/numerics/vector_tools.h>

#include <math.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>

#include <aspect/simulator.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Postprocess
  {

      template <int dim>
      void
      RichardsonExtrapolation<dim>::initialize() {
          output_file_name = this->get_output_directory()
                             + output_file_name + "_"
                             + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()))
                             + ".dat";
          input_file_name = this->get_output_directory()
                            + input_file_name + "_"
                            + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()))
                            + ".dat";
          if(Utilities::fexists(input_file_name))
            read_in_data();
      }

      template <int dim>
      void
      RichardsonExtrapolation<dim>::read_in_data()
      {
          /**
         * Assuming that the number of active cells * number of quadrature points = number of readin solution components,
         * we then allocate the input vectors.
         */
          quadrature_points_input = new std::vector<Point<dim>>;
          temperature_input = new std::vector<double>;
          pressure_input = new std::vector<double>;
          velocity_input = new std::vector<Tensor<1,dim>>;
          weight_input = new std::vector<double>;

          std::ifstream input(input_file_name);
          std::string line;

          while (std::getline(input, line))
          {
              Point<dim> q_point;
              Tensor<1,dim> velocity;
              double temperature, pressure, weight;
              std::istringstream iss(line);
              if (! (iss >> q_point[0] >> q_point[1] >> velocity[0] >> velocity[1] >> pressure >> temperature >> weight))
                  break;
              quadrature_points_input->push_back(q_point);
              velocity_input->push_back(velocity);
              pressure_input->push_back(pressure);
              temperature_input->push_back(temperature);
              weight_input->push_back(weight);
          }

          input.close();

          /**
         * For debugging purposes, we now write out what;s been read in.
         */

          std::ofstream debugging;
          debugging.open(output_file_name + ".g", std::ios_base::out);

          typename std::vector<double>::const_iterator itr_temperature = temperature_input->begin();
          typename std::vector<double>::const_iterator itr_pressure = pressure_input->begin();
          typename std::vector<Tensor<1,dim>>::const_iterator itr_velocity = velocity_input->begin();
//            typename std::vector<std::vector<double>>::const_iterator itr_compositional_fields = interpolated_compositional_fields.begin();

          typename std::vector<Point<dim>>::const_iterator itr_quadrature_points = quadrature_points_input->begin();
          typename std::vector<double>::const_iterator itr_weights = weight_input->begin();

          for (; itr_quadrature_points != quadrature_points_input->end(); itr_quadrature_points++, itr_velocity++, itr_pressure++, itr_temperature++, itr_weights++)
          {
              // Add metadata
              debugging << (*itr_quadrature_points)[0] << "\t" << (*itr_quadrature_points)[1]
                        << "\t" << (*itr_velocity)[0] << "\t" << (*itr_velocity)[1]
                        << "\t" << *itr_pressure << "\t" << *itr_temperature
                        << "\t" << *itr_weights << std::endl;
          }
          debugging.close();
      }

      template <int dim>
      void
      RichardsonExtrapolation<dim>::write_out_data()
      {
          std::ofstream interpolated_data_stream;
          interpolated_data_stream.open(output_file_name, std::ios_base::out);

          /**
          * Compute the Legendre gauss points at level 2 indirection.
          **/
          QGaussLobatto<1> base_quadrature(2);
          QIterated<dim> quadrature_rule (base_quadrature, 2);

          FEValues<dim> fe_values(this->get_mapping(),
                                  this->get_fe(),
                                  quadrature_rule,
                                  update_values |
                                  update_quadrature_points |
                                  update_JxW_values);

          const FEValuesExtractors::Scalar &extractor_pressure = this->introspection().extractors.pressure;
          const FEValuesExtractors::Scalar &extractor_temperature = this->introspection().extractors.temperature;
          const FEValuesExtractors::Vector &extractor_velocity = this->introspection().extractors.velocities;

          // iterate over all active cells on local mpi process
          for (typename DoFHandler<dim>::active_cell_iterator
                       cell = this->get_dof_handler().begin_active();
               cell != this->get_dof_handler().end();
               ++cell)
          {
              if (!(cell->is_locally_owned()))
                  continue;

              fe_values.reinit(cell);
              // Each vector is of the same length.
              const std::vector<Point<dim>>  quadrature_points = fe_values.get_quadrature_points();
              const std::vector<double>  jacobian_weight_points = fe_values.get_JxW_values();

              std::vector<double> interpolated_temperature(quadrature_points.size());
              std::vector<double> interpolated_pressure(quadrature_points.size());
              std::vector<Tensor<1,dim>> interpolated_velocity(quadrature_points.size());
//                std::vector<std::vector<double>> interpolated_compositional_fields(quadrature_points.size());

              fe_values[extractor_pressure].get_function_values(this->get_solution(), interpolated_pressure);
              fe_values[extractor_temperature].get_function_values(this->get_solution(), interpolated_temperature);
              fe_values[extractor_velocity].get_function_values(this->get_solution(), interpolated_velocity);
//                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
//                    fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
//                                                                                                            interpolated_compositional_fields[c]);

              typename std::vector<double>::const_iterator itr_temperature = interpolated_temperature.begin();
              typename std::vector<double>::const_iterator itr_pressure = interpolated_pressure.begin();
              typename std::vector<Tensor<1,dim>>::const_iterator itr_velocity = interpolated_velocity.begin();
//                typename std::vector<std::vector<double>>::const_iterator itr_compositional_fields = interpolated_compositional_fields.begin();

              typename std::vector<Point<dim>>::const_iterator itr_quadrature_points = quadrature_points.begin();
              typename std::vector<double>::const_iterator itr_weights = jacobian_weight_points.begin();

//                interpolated_data_stream << cell->active_cell_index() << std::endl;

              for (; itr_quadrature_points != quadrature_points.end(); itr_quadrature_points++, itr_velocity++, itr_pressure++, itr_temperature++, itr_weights++)
              {
                  // Add metadata
                  interpolated_data_stream << (*itr_quadrature_points)[0] << "\t" << (*itr_quadrature_points)[1]
                                           << "\t" << (*itr_velocity)[0] << "\t" << (*itr_velocity)[1]
                                           << "\t" << *itr_pressure << "\t" << *itr_temperature
                                           << "\t" << *itr_weights << std::endl;
              }
          }

          interpolated_data_stream.close();
      }

    template <int dim>
    std::pair<std::string,std::string>
    RichardsonExtrapolation<dim>::execute (TableHandler &)
    {
      std::string info = "";

      if ( this->get_time() == end_time)
      {
        if (Utilities::fexists(input_file_name)) {
            double velocity_l2_error = 0;
            double pressure_l2_error = 0;
//            double temperature_l2_error = 0;

            int n_points = 0;

            const QGauss<dim> quadrature_formula(this->get_parameters().stokes_velocity_degree);

            FEValues<dim> fe_values(this->get_mapping(),
                                    this->get_fe(),
                                    quadrature_formula,
                                    update_values |
                                    update_quadrature_points |
                                    update_JxW_values);

//          const FEValuesExtractors::Scalar &extractor_pressure = this->introspection().extractors.pressure;
//          const FEValuesExtractors::Scalar &extractor_temperature = this->introspection().extractors.temperature;
            const FEValuesExtractors::Vector &extractor_velocity = this->introspection().extractors.velocities;

            // iterate over all active cells on local mpi process
            typename std::vector<double>::const_iterator itr_temperature = temperature_input->begin();
            typename std::vector<double>::const_iterator itr_pressure = pressure_input->begin();
            typename std::vector<Tensor<1, dim>>::const_iterator itr_velocity = velocity_input->begin();
//            typename std::vector<std::vector<double>>::const_iterator itr_compositional_fields = interpolated_compositional_fields.begin();

            typename std::vector<Point<dim>>::const_iterator itr_quadrature_points = quadrature_points_input->begin();
            typename std::vector<double>::const_iterator itr_weights = weight_input->begin();

            for (; itr_quadrature_points !=
                   quadrature_points_input->end(); itr_quadrature_points++, itr_velocity++, itr_pressure++, itr_temperature++, itr_weights++) {
                  std::pair<const typename DoFHandler<dim>::active_cell_iterator,
                          Point<dim> > it = GridTools::find_active_cell_around_point(this->get_mapping(), this->get_dof_handler(), *itr_quadrature_points);
                  if (! it.first->is_locally_owned())
                      continue;

//                  std::vector<double> interpolated_temperature(fe_values.n_quadrature_points);
//                  std::vector<double> interpolated_pressure(fe_values.n_quadrature_points);
                  std::vector<Tensor<1,dim>> interpolated_velocity(fe_values.n_quadrature_points);

                  fe_values.reinit(it.first);

                  const std::vector<Point<dim>>  quadrature_points = fe_values.get_quadrature_points();
                  const std::vector<double>  jacobian_weight_points = fe_values.get_JxW_values();

//                  fe_values[extractor_pressure].get_function_values(this->get_solution(), interpolated_pressure);
//                  fe_values[extractor_temperature].get_function_values(this->get_solution(), interpolated_temperature);
                  fe_values[extractor_velocity].get_function_values(this->get_solution(), interpolated_velocity);

                  int index = 0;
                  for (int i=0; i<fe_values.n_quadrature_points; i++)
                  {
                      if ( *itr_quadrature_points == quadrature_points[i])
                          index = i;
                  }

                  velocity_l2_error += (interpolated_velocity[index] - (*itr_velocity))*(interpolated_velocity[index] - (*itr_velocity))*(*itr_weights)*(*itr_weights);
//                  pressure_l2_error += (interpolated_pressure[index] - (*itr_pressure))*(interpolated_pressure[index] - (*itr_pressure));
                  n_points++;
            }

          n_points = Utilities::MPI::sum(n_points, this->get_mpi_communicator());
          velocity_l2_error = std::sqrt(Utilities::MPI::sum(velocity_l2_error, this->get_mpi_communicator()));
//          pressure_l2_error = std::sqrt(Utilities::MPI::sum(pressure_l2_error, this->get_mpi_communicator())/n_points);

          if(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0) {
              info = "u_l2: " + std::to_string(velocity_l2_error);
              std::cout << std::setprecision(14) << "New: u_l2: " << velocity_l2_error << " p_l2: " << pressure_l2_error << std::endl;
          }
        }

        // Write out the current solution, interpolated at a higher resolved mesh.
        write_out_data();
      }

      return std::pair<std::string, std::string> ("Richardson Extrapolation: ",
                                                  info);
    }


    // Add function to compute the error between interpolated read in values with the current solution at current nodal points.

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
                prm.declare_entry("Input file name of interpolated data", "interpolated_data",
                                  Patterns::Anything (),
                                  "The file name containing the interpolated solution at the nodal values "
                                  "at the current grid resolution.");
                prm.declare_entry("Output file name of interpolated data", "interpolated_data",
                                  Patterns::Anything (),
                                  "The file name containing the interpolated solution at the nodal values "
                                  "at the current grid resolution * 4 cells.");
            }
            prm.leave_subsection();
        }
        prm.leave_subsection();
    }

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

                input_file_name = prm.get("Input file name of interpolated data");
                output_file_name = prm.get("Output file name of interpolated data");
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
