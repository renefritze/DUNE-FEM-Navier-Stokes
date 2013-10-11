
#include "main.hh"
#include <dune/navier/problems.hh>

/** \brief one single application of the discretisation and solver

  \param  mpicomm
      mostly useless atm, but mandatory
  \param  refine_level_factor
      integer to be multiplied by Dune::DGFGridInfo< GridType >::refineStepsForHalf()
      to get the used refine level for the constructed grid
  \param  stabil_coeff
      the set of coefficients to be used in the run. Default is used in all run types but StabRun().

**/
Stuff::RunInfoTimeMap singleRun(CollectiveCommunication& mpicomm, const int refine_level_factor, const int scheme_type);
//! output alert for neg. EOC
// void eocCheck( const RunInfoVector& runInfos );

bool setSchemeTypeFromString();

/**
 *  \brief  main function
 *
 *  ParameterContainer and Logger setup, select run type from parameter file and go
 *
 *  \param  argc
 *          number of arguments from command line
 *  \param  argv
 *          array of arguments from command line
 **/
int main(int argc, char** argv) {
  CollectiveCommunication mpicomm(init(argc, argv)); //( Dune::MPIManager::helper().getCommunicator() );

  if (setSchemeTypeFromString())
    Logger().Info() << "overrode scheme id from string" << std::endl;

  int err = 0;
  const unsigned int minref = Parameters().getParam("minref", 0, Dune::ValidateNotLess<int>(0));
  Stuff::RunInfoTimeMapMap rf;
  const int runtype = Parameters().getParam("runtype", 5);
  switch (runtype) {
    case 8: {
      Logger().Info() << "Reynolds runs\n";
      const int dt_steps = Parameters().getParam("dt_steps", 3, Dune::ValidateNotLess<int>(2));
      profiler().Reset(dt_steps - 1);
      int current_step = 0;
      Stuff::LoopTimer<int, Stuff::Logging::LogStream> loop_timer(current_step, dt_steps, Logger().Info());
      for (double viscosity = Parameters().getParam("viscosity", 0.1, Dune::ValidateNotLess<double>(0.0));
           dt_steps > current_step; ++loop_timer) {
        rf[current_step] = singleRun(mpicomm, minref, Parameters().getParam("scheme_type", 1, true));
        assert(rf.size());
        rf[current_step].begin()->second.refine_level = minref; // just in case the key changes from ref to sth else
        profiler().NextRun();
        viscosity /= 10.0f;
        Parameters().setParam("viscosity", viscosity);
      }
      break;
    }
    case 6: {
      Logger().Info() << "Time refine runs\n";
      const int dt_steps = Parameters().getParam("dt_steps", 3, Dune::ValidateNotLess<int>(2));
      profiler().Reset(dt_steps - 1);
      int current_step = 0;
      Stuff::LoopTimer<int, Stuff::Logging::LogStream, Stuff::QuadraticWeights> loop_timer(current_step, dt_steps,
                                                                                           Logger().Info());
      for (double dt = Parameters().getParam("fem.timeprovider.dt", 0.1, Dune::ValidateNotLess<double>(0.0));
           dt_steps > current_step; ++loop_timer) {
        rf[current_step] = singleRun(mpicomm, minref, Parameters().getParam("scheme_type", 1, true));
        assert(rf.size());
        rf[current_step].begin()->second.refine_level = minref; // just in case the key changes from ref to sth else
        profiler().NextRun();
        dt /= 2.0f;
        Parameters().setParam("fem.timeprovider.dt", dt);
      }
      break;
    }
    case 7: {
      Logger().Info() << "Scheme runs\n";
      profiler().Reset(4);
      int current_scheme = 2;
      Stuff::LoopTimer<int, Stuff::Logging::LogStream> loop_timer(current_scheme, 5, Logger().Info());
      for (; current_scheme < 6; ++loop_timer) {
        rf[current_scheme] = singleRun(mpicomm, minref, current_scheme);
        assert(rf.size());
        rf[current_scheme].begin()->second.refine_level = minref; // just in case the key changes from ref to sth else
        profiler().NextRun();
      }
      break;
    }
    case 5:
      Parameters().setParam("maxref", minref); // only one run with ref=minref
    case 0:
    default: {
      // ensures maxref>=minref
      const unsigned int maxref = Stuff::clamp(Parameters().getParam("maxref", (unsigned int)(0)), minref,
                                               Parameters().getParam("maxref", (unsigned int)(0)));
      profiler().Reset(maxref - minref + 1);
      Logger().Info() << "Grid refine runs\n";
      unsigned int ref = minref;
      Stuff::LoopTimer<unsigned int, Stuff::Logging::LogStream, Stuff::LinearWeights> loop_timer(
          ref, maxref - minref + 1, Logger().Info());
      for (; ref <= maxref; ++loop_timer) {
        rf[ref] = singleRun(mpicomm, ref, Parameters().getParam("scheme_type", 1, true));
        rf[ref].begin()->second.refine_level = ref; // just in case the key changes from ref to sth else
        profiler().NextRun();
      }
      break;
    }
  }
  profiler().OutputMap(mpicomm, rf);

  if (NAVIER_DATA_NAMESPACE::hasExactSolution && Parameters().getParam("calculate_errors", true)) {
    Stuff::TimeSeriesOutput out(rf);
    out.writeTex(Parameters().getParam("fem.io.datadir", std::string(".")) + std::string("/timeseries"));
  }

  Logger().Dbg() << "\nRun from: " << commit_string << std::endl;
  return err;
}

Stuff::RunInfoTimeMap singleRun(CollectiveCommunication& mpicomm, const int refine_level_factor,
                                const int scheme_type) {
  Stuff::Profiler::ScopedTiming pf_t("SingleRun");
  Stuff::Logging::LogStream& infoStream = Logger().Info();
  Stuff::Logging::LogStream& debugStream = Logger().Dbg();

  infoStream << "\n- initialising grid" << std::endl;
  const int gridDim = GridType::dimensionworld;
  Dune::GridPtr<GridType> gridPtr(Parameters().DgfFilename(gridDim));
  const int refine_level = (refine_level_factor) * Dune::DGFGridInfo<GridType>::refineStepsForHalf();
  gridPtr->globalRefine(refine_level);

  const int polOrder = POLORDER;
  debugStream << "  - polOrder: " << polOrder << std::endl;
  //    const double grid_width = Dune::GridWidth::calcGridWidth( gridPart );
  //    infoStream << (boost::format("  - max grid width: %f\n") % grid_width) << std::endl;

  try {
    return ThetaschemeRunner<GridType, CollectiveCommunication>(*gridPtr, mpicomm).run(scheme_type);
  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  }
  catch (std::bad_alloc& b) {
    std::cerr << "Memory allocation failed: " << b.what();
    Logger().Info().Resume();
    Stuff::meminfo(Logger().Info());
  }
  catch (assert_exception& a) {
    std::cerr << "Exception thrown at:\n" << a.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
  return Stuff::RunInfoTimeMap();
}

bool setSchemeTypeFromString() {
  bool changed = false;
  if (Dune::Parameter::exists("scheme_type_string")) {
    const std::vector<std::string>& scheme_names = Dune::NavierStokes::ThetaSchemeDescription<0>::scheme_names;
    Stuff::ValidateInList<std::string> validator(scheme_names);
    const std::string scheme_string = Parameters().getParam("scheme_type_string", scheme_names.back(), validator);
    const int scheme_id = Stuff::getIdx(scheme_names, scheme_string);
    Parameters().setParam("scheme_type", scheme_id);
    changed = true;
  }
  return changed;
}

/** Copyright (c) 2012, Rene Milk
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the FreeBSD Project.
**/
