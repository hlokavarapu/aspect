/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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


#ifndef __aspect__model_modules_h
#define __aspect__model_modules_h

#include <aspect/global.h>

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/table.h>
#include <deal.II/base/function_lib.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;


    namespace Modules
    {

      class MaterialLookup
      {
        public:
          MaterialLookup(const std::string &filename,
                         const bool interpol);

          double
          specific_heat(double temperature,
                        double pressure) const;

          double
          density(double temperature,
                  double pressure) const;

          double
          thermal_expansivity(const double temperature,
                              const double pressure) const;

          double
          seismic_Vp(const double temperature,
                     const double pressure) const;

          double
          seismic_Vs(const double temperature,
                     const double pressure) const;

          double
          dHdT (const double temperature,
                const double pressure) const;

          double
          dHdp (const double temperature,
                const double pressure) const;

          double
          dRhodp (const double temperature,
                  const double pressure) const;

          double
          value (const double temperature,
                 const double pressure,
                 const Table<2,double> &values,
                 bool interpol) const;

        private:
          double get_nT(double temperature) const;
          double get_np(double pressure) const;

          Table<2,double> density_values;
          Table<2,double> thermal_expansivity_values;
          Table<2,double> specific_heat_values;
          Table<2,double> vp_values;
          Table<2,double> vs_values;
          Table<2,double> enthalpy_values;

          double delta_press;
          double min_press;
          double max_press;
          double delta_temp;
          double min_temp;
          double max_temp;
          unsigned int numtemp;
          unsigned int numpress;
          bool interpolation;
      };


      class LateralViscosityLookup
      {
        public:
          LateralViscosityLookup(const std::string &filename);

          double lateral_viscosity(double depth);

          int get_nslices() const;
        private:
          std::vector<double> values;
          double min_depth;
          double delta_depth;
          double max_depth;
      };


      class RadialViscosityLookup
      {
        public:
          RadialViscosityLookup(const std::string &filename);

          double radial_viscosity(const double depth) const;

        private:
          std_cxx1x::shared_ptr<Functions::InterpolatedTensorProductGridData<1> > data;
      };
  }
}
}

#endif
