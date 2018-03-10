#include "../time_dependent_annulus_01.h"


namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

      namespace AnalyticSolutions
      {

          template <int dim>
          void analytic_comp_solution(
                  double pos[],
                  double vel[],
                  double *pressure,
                  double *density,
                  std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > pressure_function,
                  std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > velocity_function,
                  std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > density_function)
          {
            /****************************************************************************************/
            /****************************************************************************************/
            /* Output */
            for (unsigned int i=0; i < dim; i++)
              vel[i] = velocity_function->value(Point<dim>(pos[0],pos[1]), i);

            (*pressure) = pressure_function->value(Point<dim>(pos[0], pos[1]));

            (*density) = density_function->value(Point<dim>(pos[0], pos[1]));
          }

          /**
           * The exact solution for the benchmark, given
           * density $\rho$.
           */
          template <int dim>
          class FunctionStreamlineComp : public FunctionStreamline<dim>
          {
          public:
              FunctionStreamlineComp (std_cxx11::shared_ptr<Functions::ParsedFunction<dim>> pressure,
                                  std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > velocity,
                                  std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > density)
                      :
                      FunctionStreamline<dim>(pressure, velocity),
                      density_function (density)
              {}

              virtual void vector_value (const Point< dim > &p,
                                         Vector< double >   &values) const
              {
                double pos[2]= {p(0),p(1)};

                AnalyticSolutions::analytic_comp_solution<dim>
                        (pos,
                         &values[0], &values[2], &values[4], this->pressure_function, this->velocity_function, density_function);
              }
          private:
              std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > density_function;
          };
      }

    template <int dim>
    class TimeDependentAnnulusCompositionalMaterialModel : public TimeDependentAnnulusMaterialModel<dim>
    {
      private:
        std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > density_function, pressure_function, velocity_function;

      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.position.size(); ++i)
            {
              out.densities[i] = in.composition[i][0];
              out.viscosities[i] = 1;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;

            }
        }
    };

      /**
        * A postprocessor that evaluates the accuracy of the solution.
        *
        */
      template <int dim>
      class TimeDependentAnnulusCompPostprocessor : public TimeDependentAnnulusPostprocessor<dim>
      {
      public:
          /**
           * Generate graphical output from the current solution.
           */
          virtual
          std::pair<std::string,std::string>
          execute (TableHandler &)
          {
            std_cxx1x::shared_ptr<Function<dim> > ref_func;
            if (dynamic_cast<const TimeDependentAnnulusCompositionalMaterialModel<dim> *>(&this->get_material_model()) != NULL)
            {
              const TimeDependentAnnulusCompositionalMaterialModel<dim> *
                      material_model
                      = dynamic_cast<const TimeDependentAnnulusCompositionalMaterialModel<dim> *>(&this->get_material_model());

              ref_func.reset (new AnalyticSolutions::FunctionStreamlineComp<dim>(material_model->get_pressure(), material_model->get_velocity(), material_model->get_density()));
            }
            else
            {
              AssertThrow(false,
                          ExcMessage("Postprocessor TimeDependentAnnulusCompPostprocessor only works with the material model TimeDependentAnnulusCompositionalMaterialModel."));
            }

            const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

            Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
            Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
            Vector<float> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
            Vector<float> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());
            Vector<float> cellwise_errors_rho (this->get_triangulation().n_active_cells());
            Vector<float> cellwise_errors_rhol2 (this->get_triangulation().n_active_cells());
            Vector<float> cellwise_errors_rholinf (this->get_triangulation().n_active_cells());


            ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                                this->get_fe().n_components());
            ComponentSelectFunction<dim> comp_p(dim, this->get_fe().n_components());
            ComponentSelectFunction<dim> comp_rho(dim+2, this->get_fe().n_components());

            VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                               this->get_solution(),
                                               *ref_func,
                                               cellwise_errors_u,
                                               quadrature_formula,
                                               VectorTools::L1_norm,
                                               &comp_u);
            VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                               this->get_solution(),
                                               *ref_func,
                                               cellwise_errors_p,
                                               quadrature_formula,
                                               VectorTools::L1_norm,
                                               &comp_p);
            VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                               this->get_solution(),
                                               *ref_func,
                                               cellwise_errors_ul2,
                                               quadrature_formula,
                                               VectorTools::L2_norm,
                                               &comp_u);
            VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                               this->get_solution(),
                                               *ref_func,
                                               cellwise_errors_pl2,
                                               quadrature_formula,
                                               VectorTools::L2_norm,
                                               &comp_p);
            VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                               this->get_solution(),
                                               *ref_func,
                                               cellwise_errors_rho,
                                               quadrature_formula,
                                               VectorTools::L1_norm,
                                               &comp_rho);
            VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                               this->get_solution(),
                                               *ref_func,
                                               cellwise_errors_rhol2,
                                               quadrature_formula,
                                               VectorTools::L2_norm,
                                               &comp_rho);

            VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                               this->get_solution(),
                                               *ref_func,
                                               cellwise_errors_rholinf,
                                               quadrature_formula,
                                               VectorTools::Linfty_norm,
                                               &comp_rho);



            const double u_l1 = Utilities::MPI::sum(cellwise_errors_u.l1_norm(),this->get_mpi_communicator());
            const double p_l1 = Utilities::MPI::sum(cellwise_errors_p.l1_norm(),this->get_mpi_communicator());
            const double rho_l1 = Utilities::MPI::sum(cellwise_errors_rho.l1_norm(), this->get_mpi_communicator());
            const double u_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_ul2.norm_sqr(),this->get_mpi_communicator()));
            const double p_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_pl2.norm_sqr(), this->get_mpi_communicator()));
            const double rho_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_rhol2.norm_sqr(), this->get_mpi_communicator()));
            const double rho_linf = Utilities::MPI::max(cellwise_errors_rholinf.linfty_norm(), this->get_mpi_communicator());

            std::ostringstream os;
            os << std::scientific << u_l1
               << ", " << p_l1
               << ", " << u_l2
               << ", " << p_l2
               << ", " << rho_l1
               << ", " << rho_l2
               << ", " << rho_linf;

            return std::make_pair("Errors u_L1, p_L1, u_L2, p_L2, rho_L1, rho_L2, rho_Linf:", os.str());
          }
      };
  }
}
