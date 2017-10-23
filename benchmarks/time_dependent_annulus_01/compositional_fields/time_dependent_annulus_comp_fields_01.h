#include "../time_dependent_annulus_01.h"


namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

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
  }
}
