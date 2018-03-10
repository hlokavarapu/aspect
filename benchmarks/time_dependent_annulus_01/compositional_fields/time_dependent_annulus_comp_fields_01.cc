#include "time_dependent_annulus_comp_fields_01.h"

// explicit instantiations
namespace aspect
{
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(TimeDependentAnnulusCompositionalMaterialModel,
                                   "TimeDependentAnnulusCompositionalMaterialModel",
                                   "Material model for the time dependent annulus 01 benchmark.")


    ASPECT_REGISTER_POSTPROCESSOR(TimeDependentAnnulusCompPostprocessor,
                                  "TimeDependentAnnulusCompPostprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "derived analytical solution with the one computed by ASPECT "
                                  "and reports the error. Specifically, it can also compute the errors for "
                                  "the SolCx, SolKz and inclusion benchmarks. The postprocessor inquires "
                                  "which material model is currently being used and adjusts "
                                  "which exact solution to use accordingly.")


  }
}
