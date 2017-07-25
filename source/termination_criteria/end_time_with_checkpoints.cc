/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/termination_criteria/end_time_with_checkpoints.h>


namespace aspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    bool
    EndTimeWithCheckpoints<dim>::execute()
    {
      return (this->get_time() > end_time);
    }

//    template <int dim>
//    double EndTimeWithCheckpoints<dim>::check_for_last_time_step (const double time_step) const
//    {
//      if ((this->get_time()<end_time)
//          &&
//          (this->get_time()+time_step > end_time))
//        return end_time - this->get_time();
//      else
//        return time_step;
//    }

    template <int dim>
    void
    EndTimeWithCheckpoints<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.declare_entry ("Checkpoint times",
                         /* boost::lexical_cast<std::string>(std::numeric_limits<double>::max() /
                                              year_in_seconds) = */ "5.69e+300",
                         Patterns::List(Patterns::Double(0)),
                         "A list of checkpoint times times of the simulation. The default value is a number "
                         "so that when converted from years to seconds it is approximately "
                         "equal to the largest number representable in floating point "
                         "arithmetic. For all practical purposes, this equals infinity. "
                         "Units: Years if the "
                         "'Use years in output instead of seconds' parameter is set; "
                         "seconds otherwise.");
    }


    template <int dim>
    void
    EndTimeWithCheckpoints<dim>::parse_parameters (ParameterHandler &prm)
    {
      // read end time from parameter file. if it is to be interpreted
      // in years rather than seconds, then do the conversion
        end_time = prm.get_double ("End time");

        prm.enter_subsection("Termination criteria");
        {
            prm.enter_subsection("End step with checkpoints");
            {
                std::vector<double> checkpoint_times = Utilities::string_to_double(
                        Utilities::split_string_list(
                                prm.get("Checkpoint times")));
            }
            prm.leave_subsection ();
        }
        prm.leave_subsection ();

        if (prm.get_bool ("Use years in output instead of seconds") == true)
        {
            end_time *= year_in_seconds;
            for (std::vector<double>::iterator itr = checkpoint_times.begin();
                    itr != checkpoint_times.end(); itr++)
              *itr *= year_in_seconds;
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace TerminationCriteria
  {
    ASPECT_REGISTER_TERMINATION_CRITERION(EndTimeWithCheckpoints,
                                          "end time vector",
                                          "Terminate the simulation once the end time "
                                          "specified in the input file has been reached. "
                                          "Unlike all other termination criteria, this "
                                          "criterion is \\textit{always} active, whether it "
                                          "has been explicitly selected or not in the input file "
                                          "(this is done to preserve historical behavior of "
                                          "\\aspect{}, but it also likely does not inconvenience "
                                          "anyone since it is what would be selected in most "
                                          "cases anyway).")
  }
}
