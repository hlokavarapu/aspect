#ifndef _aspect_solkz_h
#define _aspect_solkz_h

#include <aspect/material_model/simple.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  /**
  * This is the "Sol Kz" benchmark defined in the following paper:
  * @code
  *  @Article{DMGT11,
  *    author =       {T. Duretz and D. A. May and T. V. Gerya and P. J. Tackley},
  *    title =        {Discretization errors and free surface stabilization in the
  *                  finite difference and marker-in-cell method for applied
  *                  geodynamics: {A} numerical study},
  *    journal =      {Geochemistry Geophysics Geosystems},
  *    year =         2011,
  *    volume =       12,
  *    pages =        {Q07004/1--26}}
  * @endcode
  *
  * The results are published in Kronbichler, Heister and Bangerth paper.
  */
  namespace InclusionBenchmark
  {
    using namespace dealii;

    namespace AnalyticSolutions
    {
      /**
      * The exact solution for the SolKz benchmark.
      */
      template<int dim>
      class FunctionSolKz : public Function<dim>
      {
        public:
          FunctionSolKz() : Function<dim>() {}

          virtual void vector_value(const Point<dim> &p,
                                    Vector<double> &values) const;
      };
    }

    template <int dim>
    class SolKzMaterial : public MaterialModel::InterfaceCompatibility<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual inline double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        virtual inline double density (const double temperature,
                                const double pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         * Incompressibility does not necessarily imply that the density is
         * constant; rather, it may still depend on temperature or pressure.
         * In the current context, compressibility means whether we should
         * solve the contuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;
        /**
         * @}
         */

        void
        parse_parameters (ParameterHandler &prm);
    };

    /**
    * A postprocessor that evaluates the accuracy of the solution.
    *
    * The implementation of error evaluators that correspond to the
    * benchmarks defined in the paper Duretz et al. reference above.
    */
    template <int dim>
    class SolKzPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}
#endif
