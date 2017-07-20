#include "simple_annulus.h"

// explicit instantiations
namespace aspect
{
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SimpleAnnulusMaterialModel,
                                   "SimpleAnnulusMaterial",
                                   "EGP and HVL benchmark material model.")


    ASPECT_REGISTER_POSTPROCESSOR(SimpleAnnulus2Postprocessor,
                                  "SimpleAnnulus2Postprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "derived analytical solution with the one computed by ASPECT "
                                  "and reports the error. Specifically, it can also compute the errors for "
                                  "the SolCx, SolKz and inclusion benchmarks. The postprocessor inquires "
                                  "which material model is currently being used and adjusts "
                                  "which exact solution to use accordingly.")


  }
}
