#include "main.hh"

#include "oseen.hh"

/** \brief one single application of the discretisation and solver

  \param  mpicomm
      mostly useless atm, but mandatory
  \param  refine_level_factor
      integer to be multiplied by Dune::DGFGridInfo< GridType >::refineStepsForHalf()
      to get the used refine level for the constructed grid
  \param  stabil_coeff
      the set of coefficients to be used in the run. Default is used in all run types but StabRun().

**/
Stuff::RunInfoVector singleRun(CollectiveCommunication& mpicomm, int refine_level_factor);

//! output alert for neg. EOC
void eocCheck(const Stuff::RunInfoVector& runInfos);

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
#ifdef NDEBUG
  try
#endif
  {
    CollectiveCommunication mpicomm(init(argc, argv));

    int err = 0;

    const int minref = Parameters().getParam("minref", 0);
    // ensures maxref>=minref
    const int maxref = Stuff::clamp(Parameters().getParam("maxref", 0), minref, Parameters().getParam("maxref", 0));
    profiler().Reset(maxref - minref + 1);
    for (int ref = minref; ref <= maxref; ++ref) {
      singleRun(mpicomm, ref);
      profiler().NextRun();
    }

    Logger().Dbg() << "\nRun from: " << commit_string << std::endl;
    return err;
  }

#ifdef NDEBUG
  catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (std::bad_alloc& b) {
    std::cerr << "Memory allocation failed: " << b.what();
    Logger().Info().Resume();
    Stuff::meminfo(Logger().Info());
  }
  catch (assert_exception& a) {
    std::cerr << "Exception thrown at:\n" << a.what() << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception thrown at:\n" << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
#endif
}

Stuff::RunInfoVector singleRun(CollectiveCommunication& mpicomm, int refine_level_factor) {
  profiler().StartTiming("SingleRun");
  Stuff::Logging::LogStream& infoStream = Logger().Info();
  Stuff::Logging::LogStream& debugStream = Logger().Dbg();
  Stuff::RunInfoVector runInfoVector;

  /* ********************************************************************** *
   * initialize the grid                                                    *
   * ********************************************************************** */
  infoStream << "\n- initialising grid" << std::endl;
  const int gridDim = GridType::dimensionworld;
  Dune::GridPtr<GridType> gridPtr(Parameters().DgfFilename(gridDim));
  static bool firstRun = true;
  int refine_level = (refine_level_factor) * Dune::DGFGridInfo<GridType>::refineStepsForHalf();
  if (firstRun && refine_level_factor > 0) {
    refine_level = (refine_level_factor) * Dune::DGFGridInfo<GridType>::refineStepsForHalf();
    gridPtr->globalRefine(refine_level);
  }

  typedef Dune::AdaptiveLeafGridPart<GridType> GridPartType;
  GridPartType gridPart(*gridPtr);

  infoStream << "\n- initialising problem" << std::endl;

  const int polOrder = POLORDER;
  debugStream << "  - polOrder: " << polOrder << std::endl;

  const double grid_width = Dune::GridWidth::calcGridWidth(gridPart);
  infoStream << "  - max grid width: " << grid_width << std::endl;

  //	Dune::CompileTimeChecker< ( VELOCITY_POLORDER >= 2 ) > RHS_ADAPTER_CRAPS_OUT_WITH_VELOCITY_POLORDER_LESS_THAN_2;

  const double reynolds = Parameters().getParam("reynolds", 1.0);
  //	const double theta_ = 1.0;
  const double d_t = 1.0;
  //	const double operator_weight_beta_ = 1.0;
  //	const double operator_weight_alpha_ = 1.0;
  //	const double oseen_alpha = Parameters().getParam( "alpha", 1.0 );
  const double oseen_viscosity = Parameters().getParam("viscosity", 1.0);
  const double lambda = (reynolds * 0.5) - std::sqrt((std::pow(reynolds, 2) * 0.25) + (4 * std::pow(M_PI, 2)));
  //	const double pressure_C = ( std::exp( 3 * lambda ) - std::exp(-1  * lambda ) ) / ( - 8 * lambda );

  //	const double lambda = - 8 *M_PI * M_PI / ( reynolds + std::sqrt(reynolds*reynolds + 64 * M_PI * M_PI));

  Parameters().setParam("lambda", lambda);
  Parameters().setParam("viscosity", oseen_viscosity);
  Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
  stab_coeff.FactorFromParams("D12", 0);
  stab_coeff.FactorFromParams("C12", 0);
  stab_coeff.Add("E12", 0.5);
  stab_coeff.FactorFromParams("E12", 0.5);

  typedef Dune::Oseen::Traits<CollectiveCommunication, GridPartType, gridDim, polOrder, VELOCITY_POLORDER,
                              PRESSURE_POLORDER> OseenTraits;

  //	CollectiveCommunication comm;// = Dune::MPIManager::helper().getCommunicator();
  OseenTraits::TimeProviderType timeprovider_(OseenTraits::SchemeDescriptionType::crank_nicholson(0.5), mpicomm);
  OseenTraits::OseenModelTraits::DiscreteOseenFunctionSpaceWrapperType functionSpaceWrapper(gridPart);

  typedef OseenTraits::OseenModelTraits::DiscreteOseenFunctionWrapperType DiscreteOseenFunctionWrapperType;
  DiscreteOseenFunctionWrapperType currentFunctions("current_", functionSpaceWrapper, gridPart);
  DiscreteOseenFunctionWrapperType nextFunctions("next_", functionSpaceWrapper, gridPart);
  DiscreteOseenFunctionWrapperType errorFunctions("error_", functionSpaceWrapper, gridPart);
  OseenTraits::ExactSolutionType exactSolution(timeprovider_, gridPart, functionSpaceWrapper);
  exactSolution.project();
  //	exactSolution.exactPressure().setShift( pressure_C );
  //	Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();

  //	currentFunctions.assign( exactSolution );
  currentFunctions.clear();
  nextFunctions.clear();

  typedef OseenTraits::OseenModelTraits::PressureFunctionSpaceType PressureFunctionSpaceType;
  //	PressureFunctionSpaceType pressureFunctionSpace;
  //	Stuff::VolumeDiffFunction<PressureFunctionSpaceType> vol(pressureFunctionSpace, -12.8430582392842 );

  //	Dune::BetterL2Projection
  //		::project( 0.0, vol, nextFunctions.discretePressure() );
  //	double meanPressure_exact = Stuff::integralAndVolume( exactSolution.exactPressure(),
  //nextFunctions.discretePressure().space() ).first;
  double meanPressure_discrete =
      Stuff::meanValue(nextFunctions.discretePressure(), nextFunctions.discretePressure().space());
  //	double GD = Stuff::boundaryIntegral( stokesDirichletData, nextFunctions.discreteVelocity().space() );

  //	Dune::L2Norm< GridPartType > l2_Error( gridPart );

  Logger().Info().Resume();
  Logger().Info()
      //					<< "L2-Error Pressure (abs|rel): " << std::setw(8) << l2_error_pressure << " | " <<
      //relative_l2_error_pressure << "\n"
      //					<< "L2-Error Velocity (abs|rel): " << std::setw(8) << l2_error_velocity << " | " <<
      //relative_l2_error_velocity << "\n"
      << "Mean pressure (discrete): " << meanPressure_discrete << std::endl << std::endl;

  return runInfoVector;
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
