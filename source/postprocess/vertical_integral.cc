/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/postprocess/vertical_integral.h>

#include <aspect/simulator_access.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace
    {
      unsigned int
      get_index(const std::vector<unsigned int> &index,
                const std::vector<unsigned int> &size)
      {
        if ((index.size()  == 1)
            && (size.size() == 1))
          return index[0];
        else if ((index.size()  == 2)
            && (size.size() == 2))
          return index[1] * size[0] + index[0];
        else
          {
          Assert(false,
              ExcMessage("Called get_index with wrong size of arguments."));
          return 0;
          }
      }
    }

      template <int dim>
      std::pair<std::string,std::string>
      VerticalIntegral<dim>::execute (TableHandler &statistics)
      {
        // TODO: This criterion does not work for time_of_output == 0.0
        if ((this->get_time() < time_of_output)
              || (this->get_time() - this->get_old_timestep() >= time_of_output))
          return std::pair<std::string,std::string>();

        const GeometryModel::Box<dim> *
              geometry_model = dynamic_cast <const GeometryModel::Box<dim>*> (&this->get_geometry_model());
        AssertThrow (geometry_model,
                     ExcMessage("The Vertical integral postprocessor is only implemented for a box geometry."
                         "Please make sure you are using the right geometry or extend the postprocessor"));

        const Point<dim> model_origin = geometry_model->get_origin();
        const Point<dim> model_extent = geometry_model->get_extents();
        const unsigned int refinement = this->get_triangulation().n_levels() - 1;
        const unsigned int intervals = std::pow(2,refinement);
        std::vector<unsigned int> repetitions = geometry_model->get_repetitions();
        repetitions.resize(dim-1);

        // The surface grid will be centered on cell midpoints, therefore the
        // number of surface points is equal to the number of cells of the computational
        // grid
        std::vector<unsigned int>surface_grid_points(dim-1);
        unsigned int number_of_surface_points = 1;
        double surface_area = 1;
        Point<dim-1> grid_origin;
        Point<dim-1> grid_extent;

        for (unsigned int i=0; i < dim-1; i++)
          {
            grid_origin[i] = model_origin[i];
            grid_extent[i] = model_extent[i];
            surface_area  *= model_extent[i];
            surface_grid_points[i] = intervals * repetitions[i];
            number_of_surface_points *= intervals * repetitions[i];
          }

        const double surface_area_per_point = surface_area
            / number_of_surface_points;

        AssertThrow (this->introspection().compositional_name_exists(name_of_compositional_field),
                     ExcMessage("The compositional field " + name_of_compositional_field +
                                " you asked for is not used in the simulation."));
        const unsigned int compositional_index = this->introspection().compositional_index_for_name(name_of_compositional_field);

        std::vector<double> surface_grid(number_of_surface_points);

        const QMidpoint<dim> quadrature_formula;

        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_values |
                                 update_q_points |
                                 update_JxW_values);

        std::vector<double> composition(quadrature_formula.size());

        // loop over all of the surface cells and if one less than h/3 away from
        // the top surface, evaluate the stress at its center
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
              {
                fe_values.reinit (cell);

                // get the various components of the composition, then
                // evaluate the material properties there
                fe_values[this->introspection().extractors.compositional_fields[compositional_index]]
                .get_function_values(this->get_solution(),
                                     composition);

                // Compute the integral of the compositional field
                // over the entire cell, by looping over all quadrature points
                // (currently, there is only one, but the code is generic).
                for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                  {
                    const Point<dim> location = fe_values.quadrature_point(q);
                    std::vector<unsigned int> index(dim-1);
                    for (unsigned int i = 0; i<dim-1; i++)
                      index[i] = static_cast<unsigned int> (surface_grid_points[i] * (location[i]-grid_origin[i])/grid_extent[i]);

                    // JxW provides the volume quadrature weights. This is a general formulation
                    // necessary for when a quadrature formula is used that has more than one point.
                    surface_grid[get_index(index,surface_grid_points)] += composition[q] * fe_values.JxW(q);
                  }
              }

        Utilities::MPI::sum (surface_grid,
                             this->get_mpi_communicator(),
                             surface_grid);

        // Simple output routine
        const std::string filename = this->get_output_directory() +
                                     "vertically_integrated_" + name_of_compositional_field + "." +
                                     Utilities::int_to_string(this->get_timestep_number(), 5)
                                     + DataOutBase::default_suffix(output_format);
//
//        // on processor 0, collect all of the data the individual processors send
//        // and add and write them into one file
//        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
//          {
//            std::ofstream file (filename.c_str());
//            // output the field
//            for (unsigned int i = 0; i < surface_grid.size(); ++i)
//              file << surface_grid[i] << std::endl;
//          }

        // On the root process, write out the file. do this using the DataOutStack
        // class on a piecewise constant finite element space on
        // a 1d mesh with the correct subdivisions
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          {
            Triangulation<dim-1> mesh;
            GridGenerator::subdivided_hyper_rectangle (mesh,
                                                       surface_grid_points,
                                                       grid_origin,
                                                       grid_origin+grid_extent);

            FE_DGQ<dim-1> fe(0);
//            FEFaceValues<dim-1> fe_face_values (this->get_mapping(),
//                                                fe,
//                                                quadrature_formula_face,
//                                                update_JxW_values);
            DoFHandler<dim-1> dof_handler (mesh);
            dof_handler.distribute_dofs(fe);
            //Assert (dof_handler.n_dofs() == n_depth_zones, ExcInternalError());

            DataOut<dim-1> data_out;
            const std::string variables = "vertically_integrated_" + name_of_compositional_field;

            data_out.attach_dof_handler (dof_handler);

            Vector<double> tmp(number_of_surface_points);
            std::copy (surface_grid.begin(),
                       surface_grid.end(),
                       tmp.begin());
            tmp /= surface_area_per_point;
            data_out.add_data_vector (tmp,
                                      variables,
                                      DataOut<dim-1>::type_cell_data);
            data_out.build_patches ();

            std::ofstream f (filename.c_str());
            data_out.write (f, output_format);
          }

        return std::pair<std::string,std::string>("Writing vertical integral file: ",
                                                  filename);
      }

      template <int dim>
      void
      VerticalIntegral<dim>::
      declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
            prm.enter_subsection("Vertical integral");
            {
              prm.declare_entry ("Name of compositional field", "",
                                 Patterns::Anything (),
                                 "Option to remove the mean dynamic topography "
                                 "in the outputted data file (not visualization). "
                                 "'true' subtracts the mean, 'false' leaves "
                                 "the calculated dynamic topography as is. ");
              prm.declare_entry ("Time of output", "0.0",
                                 Patterns::Double (),
                                 "A parameter that denotes, at what time this "
                                 "postprocessor will be executed.");
              prm.declare_entry ("Output format", "vtu",
                                 Patterns::Selection(DataOutBase::get_output_format_names()),
                                 "The format in which the output shall be produced. The "
                                 "format in which the output is generated also determiens "
                                 "the extension of the file into which data is written.");
            }
            prm.leave_subsection();
          }
        prm.leave_subsection();
      }


      template <int dim>
      void
      VerticalIntegral<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
            prm.enter_subsection("Vertical integral");
            {
              name_of_compositional_field = prm.get("Name of compositional field");
              time_of_output         = prm.get_double("Time of output");
              output_format               = DataOutBase::parse_output_format(prm.get("Output format"));
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
      ASPECT_REGISTER_POSTPROCESSOR(VerticalIntegral,
                                    "vertical integral",
                                    "A visualization output object that integrates a given compositional "
                                    "field vertically and outputs the integrated value at the surface cells. "
                                    "The output is calculated as integrated volume per surface area, which "
                                    "is equal to the thickness of a layer containing all the material below "
                                    "each cell. The postprocessor is only implemented for a box geometry at "
                                    "the moment.")
  }
}
