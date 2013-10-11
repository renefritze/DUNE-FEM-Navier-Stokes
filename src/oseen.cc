#include "main.hh"

#define STOKES_CONV_ONLY 1

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

#include "oseen.hh"

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
  catch (std::exception& e) {
    std::cerr << "Exception thrown at:\n" << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
#endif
}

Stuff::RunInfoVector singleRun(CollectiveCommunication& comm, int refine_level_factor) {
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

  const int polOrder = POLORDER;
  typedef Dune::Oseen::Traits<CollectiveCommunication, GridType, gridDim, polOrder, VELOCITY_POLORDER,
                              PRESSURE_POLORDER> OseenTraits;

  typedef typename OseenTraits::GridPartType GridPartType;
  GridPartType gridPart(*gridPtr);

  infoStream << "\n- initialising problem" << std::endl;
  debugStream << "  - polOrder: " << polOrder << std::endl;
  const double grid_width = Dune::GridWidth::calcGridWidth(gridPart);
  infoStream << "  - max grid width: " << grid_width << std::endl;

  //	Dune::CompileTimeChecker< ( VELOCITY_POLORDER >= 2 ) > RHS_ADAPTER_CRAPS_OUT_WITH_VELOCITY_POLORDER_LESS_THAN_2;

  //	const double reynolds = Parameters().getParam( "reynolds", 1.0 );
  //	const double theta_ = 1.0;
  //	const double d_t = 1.0;
  //	const double operator_weight_beta_ = 1.0;
  //	const double operator_weight_alpha_ = 1.0;
  const double oseen_alpha = Parameters().getParam("alpha", 1.0);
  const double oseen_viscosity = Parameters().getParam("viscosity", 1.0);
  //	const double lambda = ( reynolds * 0.5 )
  //						  - std::sqrt(
  //								  ( std::pow( reynolds, 2 ) * 0.25 )
  //								  + ( 4 * std::pow( M_PI, 2 ) )
  //									  ) ;
  //	const double pressure_C = ( std::exp( 3 * lambda ) - std::exp(-1  * lambda ) ) / ( - 8 * lambda );

  //	const double lambda = - 8 *M_PI * M_PI / ( reynolds + std::sqrt(reynolds*reynolds + 64 * M_PI * M_PI));

  //	Parameters().setParam( "lambda", lambda );
  Parameters().setParam("viscosity", oseen_viscosity);
  Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
  //	stab_coeff.FactorFromParams( "D12", 0 );
  //	stab_coeff.FactorFromParams( "C12", 0 );
  //	stab_coeff.Add( "E12", 0.5 );
  //	stab_coeff.FactorFromParams( "E12", 0.5 );

  OseenTraits::TimeProviderType timeprovider_(OseenTraits::SchemeDescriptionType::crank_nicholson(0.5), comm);
  OseenTraits::OseenModelTraits::DiscreteOseenFunctionSpaceWrapperType functionSpaceWrapper(gridPart);

  typedef OseenTraits::OseenModelTraits::DiscreteOseenFunctionWrapperType DiscreteOseenFunctionWrapperType;
  DiscreteOseenFunctionWrapperType currentFunctions("current_", functionSpaceWrapper, gridPart);
  DiscreteOseenFunctionWrapperType nextFunctions("next_", functionSpaceWrapper, gridPart);
  DiscreteOseenFunctionWrapperType tmpFunctions("tmp_", functionSpaceWrapper, gridPart);
  DiscreteOseenFunctionWrapperType errorFunctions("error_", functionSpaceWrapper, gridPart);
  OseenTraits::ExactSolutionType exactSolution(timeprovider_, gridPart, functionSpaceWrapper);
  exactSolution.project();
  //	exactSolution.exactPressure().setShift( pressure_C );
  Stuff::Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
  Stuff::printDiscreteFunctionMatlabStyle(exactSolution.discretePressure(), "p_exakt", matlabLogStream);
  Stuff::printDiscreteFunctionMatlabStyle(exactSolution.discreteVelocity(), "u_exakt", matlabLogStream);

  OseenTraits::OseenModelTraits::AnalyticalDirichletDataType stokesDirichletData(timeprovider_, functionSpaceWrapper);

  //	OseenTraits::OseenModelTraits::PressureFunctionSpaceType
  //			continousPressureSpace;
  OseenTraits::OseenModelTraits::VelocityFunctionSpaceType continousVelocitySpace;

  OseenTraits::OseenModelTraits::AnalyticalForceFunctionType force(timeprovider_, continousVelocitySpace,
                                                                   oseen_viscosity, oseen_alpha);
  OseenTraits::OseenModelType stokesModel(stab_coeff, force, stokesDirichletData, oseen_viscosity, /*viscosity*/
                                          oseen_alpha,                                             /*alpha*/
                                          1, // Parameters().getParam( "cscale", 1.0 ),/*convection_scale_factor*/
                                          1 /*pressure_gradient_scale_factor*/);

  //	currentFunctions.assign( exactSolution );
  currentFunctions.clear();
  nextFunctions.clear();

  OseenTraits::ConvectionType convection(timeprovider_, continousVelocitySpace);
  DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType discrete_convection(
      "convetion", currentFunctions.discreteVelocity().space());
  Dune::L2Projection<double, double, OseenTraits::ConvectionType,
                     DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType>()(convection, discrete_convection);

  OseenTraits::OseenModelTraits::DiscreteSigmaFunctionSpaceType sigma_space(gridPart);
  Dune::Oseen::RhsDatacontainer<typename OseenTraits::OseenModelTraits> rhs_container(
      currentFunctions.discreteVelocity().space(), sigma_space);
  auto beta = exactSolution.discreteVelocity();
  typedef OSEEN_DATA_NAMESPACE::Beta<OseenTraits::OseenModelTraits::VelocityFunctionSpaceType,
                                     OseenTraits::TimeProviderType> ContBetaType;
  ContBetaType beta_cont(timeprovider_, continousVelocitySpace);
  Dune::BetterL2Projection::project(timeprovider_, beta_cont, beta);

  OseenTraits::OseenPassType oseenPass(stokesModel, gridPart, functionSpaceWrapper, beta, true /*do_osseen_disc*/);
  nextFunctions.clear();
  currentFunctions.clear();
  oseenPass.printInfo();
  oseenPass.apply(currentFunctions, nextFunctions, &rhs_container);

  errorFunctions.discretePressure().assign(exactSolution.discretePressure());
  errorFunctions.discretePressure() -= nextFunctions.discretePressure();
  errorFunctions.discreteVelocity().assign(exactSolution.discreteVelocity());
  errorFunctions.discreteVelocity() -= nextFunctions.discreteVelocity();

  double meanPressure_exact =
      Stuff::integralAndVolume(exactSolution.exactPressure(), nextFunctions.discretePressure().space()).first;
  double meanPressure_discrete =
      Stuff::meanValue(nextFunctions.discretePressure(), nextFunctions.discretePressure().space());
  //	typedef OseenTraits::OseenModelTraits::PressureFunctionSpaceType
  //			PressureFunctionSpaceType;
  //	PressureFunctionSpaceType pressureFunctionSpace;
  //	Stuff::ConstantFunction<PressureFunctionSpaceType> vol(pressureFunctionSpace, meanPressure_discrete );
  //	Dune::BetterL2Projection
  //		::project( 0.0, vol, tmpFunctions.discretePressure() );
  //	nextFunctions.discretePressure() -= tmpFunctions.discretePressure();
  double meanPressure_discrete_after =
      Stuff::meanValue(nextFunctions.discretePressure(), nextFunctions.discretePressure().space());

  double GD = Stuff::boundaryIntegral(stokesDirichletData, nextFunctions.discreteVelocity().space());

  Dune::L2Norm<GridPartType> l2_Error(gridPart);

  const double l2_error_pressure = l2_Error.norm(errorFunctions.discretePressure());
  const double l2_error_velocity = l2_Error.norm(errorFunctions.discreteVelocity());
  const double relative_l2_error_pressure = l2_error_pressure / l2_Error.norm(exactSolution.discretePressure());
  const double relative_l2_error_velocity = l2_error_velocity / l2_Error.norm(exactSolution.discreteVelocity());

  //	errorFunctions.discretePressure().assign( exactSolution.discretePressure() );
  //	errorFunctions.discretePressure() -= nextFunctions.discretePressure();

  const double l2_error_pressure_after = l2_Error.norm(errorFunctions.discretePressure());
  const double relative_l2_error_pressure_after =
      l2_error_pressure_after / l2_Error.norm(exactSolution.discretePressure());

  Logger().Info().Resume();
  Logger().Info() << "L2-Error Pressure (abs|rel): " << std::setw(8) << l2_error_pressure << " | "
                  << relative_l2_error_pressure << "\n"
                  << "L2-Error Pressure after (abs|rel): " << std::setw(8) << l2_error_pressure_after << " | "
                  << relative_l2_error_pressure_after << "\n"
                  << "L2-Error Velocity (abs|rel): " << std::setw(8) << l2_error_velocity << " | "
                  << relative_l2_error_velocity << "\n"
                  << "Mean pressure (exact|discrete|after): " << meanPressure_exact << " | " << meanPressure_discrete
                  << " | " << meanPressure_discrete_after << std::endl << "GD: " << GD << "\n"
      //					<< "lambda: " << lambda
                  << "current time: " << timeprovider_.time() << std::endl;

  typedef Dune::tuple<const DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType*,
                      const DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType*,
                      const DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType*,
                      const DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType*,
                      const DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType*,
                      const DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType*,
                      const DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType*> OutputTupleType;
  typedef Dune::NavierStokes::TimeAwareDataWriter<OseenTraits::TimeProviderType, GridPartType::GridType,
                                                  OutputTupleType> DataWriterType;
  OutputTupleType out(&nextFunctions.discreteVelocity(), &nextFunctions.discretePressure(),
                      &exactSolution.discreteVelocity(), &exactSolution.discretePressure(),
                      &errorFunctions.discreteVelocity(), &errorFunctions.discretePressure(), &discrete_convection);

  DataWriterType dt(timeprovider_, gridPart.grid(), out);
  dt.write();
  errorFunctions.discreteVelocity().write_ascii("test.ascii");
  errorFunctions.discreteVelocity().write_xdr("test.xdr");

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
