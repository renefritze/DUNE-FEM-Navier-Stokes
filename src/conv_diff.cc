#include "main.hh"
#include <dune/navier/global_defines.hh>

#include <cstdio>
#if defined(USE_PARDG_ODE_SOLVER) && defined(USE_BFG_CG_SCHEME)
#warning("USE_PARDG_ODE_SOLVER enabled, might conflict with custom solvers")
#endif

#if defined(UGGRID) && defined(DEBUG)
#warning("UGGRID in debug mode is likely to produce a segfault")
#endif

#define USE_GRPAE_VISUALISATION (HAVE_GRAPE && !defined(AORTA_PROBLEM))

#include <vector>
#include <string>

#include <iostream>
#include <cmath>
#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/common/exceptions.hh>   // We use exceptions
#include <dune/grid/common/capabilities.hh>

#include "conv_diff.hh"

#if ENABLE_MPI
typedef Dune::CollectiveCommunication<MPI_Comm> CollectiveCommunication;
#else
typedef Dune::CollectiveCommunication<double> CollectiveCommunication;
#endif

//! the strings used for column headers in tex output
typedef std::vector<std::string> ColumnHeaders;

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
  Stuff::Signals::installSignalHandler(SIGINT);
#ifdef NDEBUG
  try
#endif
  {

    Dune::MPIManager::initialize(argc, argv);
    // assert( Dune::Capabilities::isParallel< GridType >::v );
    CollectiveCommunication mpicomm(Dune::MPIManager::helper().getCommunicator());

    /* ********************************************************************** *
     * initialize all the stuff we need                                       *
     * ********************************************************************** */
    if (argc < 2) {
      std::cerr << "\nUsage: " << argv[0] << " parameterfile \n"
                << "\n\t --- OR --- \n";
      std::cerr << "\nUsage: " << argv[0] << " paramfile:"
                << "file"
                << " more-opts:val ..." << std::endl;
      std::cerr << "\nUsage: " << argv[0] << " -d paramfile "
                << "\n\t(for displaying solutions in grape) " << std::endl;
      Parameters().PrintParameterSpecs(std::cerr);
      std::cerr << std::endl;
      return 2;
    }
#if USE_GRPAE_VISUALISATION
    if (!strcmp(argv[1], "-d") || !strcmp(argv[1], "-r")) {
      return display(argc, argv);
    }
#endif
    if (!(Parameters().ReadCommandLine(argc, argv))) {
      return 1;
    }

    // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
    //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
    const bool useLogger = false;
    Logger().Create(Parameters().getParam("loglevel", 62, useLogger),
                    Parameters().getParam("logfile", std::string("dune_stokes"), useLogger),
                    Parameters().getParam("fem.io.logdir", std::string(), useLogger),
                    Parameters().getParam("fem.io.logdir", std::string(), useLogger));

    int err = 0;

    const int minref = Parameters().getParam("minref", 0);
    // ensures maxref>=minref
    const int maxref = Stuff::clamp(Parameters().getParam("maxref", 0), minref, Parameters().getParam("maxref", 0));
    profiler().Reset(maxref - minref + 1);
    Stuff::RunInfoVectorMap rf;
    for (unsigned int ref = minref; ref <= maxref; ++ref) {
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
  Logging::LogStream& infoStream = Logger().Info();
  Logging::LogStream& debugStream = Logger().Dbg();
  Stuff::RunInfoVector runInfoVector;

  /* ********************************************************************** *
   * initialize the grid                                                    *
   * ********************************************************************** */
  infoStream << "\n- initialising grid" << std::endl;
  const int gridDim = GridType::dimensionworld;
  Dune::GridPtr<GridType> gridPtr(Parameters().DgfFilename(gridDim));
  int refine_level = (refine_level_factor) * Dune::DGFGridInfo<GridType>::refineStepsForHalf();
  gridPtr->globalRefine(refine_level);

  typedef Dune::ConvDiff::Traits<CollectiveCommunication, GridType, gridDim, polOrder, VELOCITY_POLORDER,
                                 PRESSURE_POLORDER> ConvDiffTraits;

  typedef typename ConvDiffTraits::GridPartType<GridType> GridPartType;
  GridPartType gridPart(*gridPtr);

  /* ********************************************************************** *
   * initialize problem                                                     *
   * ********************************************************************** */
  infoStream << "\n- initialising problem" << std::endl;

  const int polOrder = POLORDER;
  debugStream << "  - polOrder: " << polOrder << std::endl;

  Parameters().setParam("reduced_oseen_solver", true);

  const double alpha = Parameters().getParam("alpha", 1.0);
  const double viscosity = Parameters().getParam("viscosity", 1.0);

  Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
  stab_coeff.FactorFromParams("D12", 0);
  stab_coeff.FactorFromParams("C12", 0);
  stab_coeff.Add("E12", 0.5);

  CollectiveCommunication comm = Dune::MPIManager::helper().getCommunicator();

  ConvDiffTraits::TimeProviderType timeprovider_(ConvDiffTraits::SchemeDescriptionType::crank_nicholson(0.5), comm);
  ConvDiffTraits::OseenModelTraits::DiscreteOseenFunctionSpaceWrapperType functionSpaceWrapper(gridPart);

  typedef ConvDiffTraits::OseenModelTraits::DiscreteOseenFunctionWrapperType DiscreteOseenFunctionWrapperType;
  DiscreteOseenFunctionWrapperType currentFunctions("current_", functionSpaceWrapper, gridPart);
  DiscreteOseenFunctionWrapperType nextFunctions("next_", functionSpaceWrapper, gridPart);
  DiscreteOseenFunctionWrapperType errorFunctions("error_", functionSpaceWrapper, gridPart);
  ConvDiffTraits::ExactSolutionType exactSolution(timeprovider_, gridPart, functionSpaceWrapper);
  exactSolution.project();

  Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
  //	Stuff::printDiscreteFunctionMatlabStyle( exactSolution.discretePressure(), "p_exakt", matlabLogStream );
  //	Stuff::printDiscreteFunctionMatlabStyle( exactSolution.discreteVelocity(), "u_exakt", matlabLogStream );

  ConvDiffTraits::StartPassType startPass;
  ConvDiffTraits::OseenModelTraits::AnalyticalDirichletDataType stokesDirichletData =
      ConvDiffTraits::OseenModelTraits::AnalyticalDirichletDataTraitsImplementation::getInstance(timeprovider_,
                                                                                                 functionSpaceWrapper);

  ConvDiffTraits::OseenModelTraits::VelocityFunctionSpaceType continousVelocitySpace;

  ConvDiffTraits::OseenModelTraits::AnalyticalForceFunctionType force(viscosity, continousVelocitySpace, alpha);
  ConvDiffTraits::OseenModelType stokesModel(stab_coeff, force, stokesDirichletData, viscosity, alpha,
                                             1.0, /*convection_scale_factor*/
                                             0.0 /*pressure_gradient_scale_factor*/);
  currentFunctions.assign(exactSolution);

  ConvDiffTraits::OseenModelTraits::DiscreteSigmaFunctionSpaceType sigma_space(gridPart);
  ConvDiffTraits::OseenModelTraits::DiscreteSigmaFunctionType discrete_velocityGradient("velocityGradient",
                                                                                        sigma_space);
  ConvDiffTraits::OseenModelTraits::SigmaFunctionSpaceType cont_sigma_space;
  ConvDiffTraits::VelocityGradientType velocityGradient(timeprovider_, cont_sigma_space);

  ConvDiffTraits::ConvectionType convection(timeprovider_, continousVelocitySpace);
  DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType discrete_convection(
      "convection", currentFunctions.discreteVelocity().space());
  DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType discrete_exactConvection(
      "exact_convection", currentFunctions.discreteVelocity().space());
  DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType convection_diff(
      "convection_diff", currentFunctions.discreteVelocity().space());
  Dune::L2Projection<double, double, ConvDiffTraits::ConvectionType,
                     DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType>()(convection, discrete_convection);

  Dune::RhsDatacontainer<typename ConvDiffTraits::OseenModelTraits> rhs_container(
      currentFunctions.discreteVelocity().space(), sigma_space);
  ConvDiffTraits::OseenPassType oseenPass(startPass, stokesModel, gridPart, functionSpaceWrapper,
                                          exactSolution.discreteVelocity(), true);
  oseenPass.apply(currentFunctions, nextFunctions, &rhs_container, &velocityGradient);

  ConvDiffTraits::ExactConvectionType exactConvection(timeprovider_, continousVelocitySpace);
  Dune::L2Projection<double, double, ConvDiffTraits::ExactConvectionType,
                     DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType>()(exactConvection,
                                                                                       discrete_exactConvection);

  typedef Stuff::GradientSplitterFunction<DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType,
                                          ConvDiffTraits::OseenModelTraits::DiscreteSigmaFunctionType>
  GradientSplitterFunctionType;
  GradientSplitterFunctionType gradient_splitter(functionSpaceWrapper.discreteVelocitySpace(),
                                                 rhs_container.velocity_gradient);

  ConvDiffTraits::VelocityGradientYType velocityGradientY(timeprovider_, continousVelocitySpace);
  DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType discrete_velocityGradientY(
      "velocityGradientY", currentFunctions.discreteVelocity().space());
  Dune::L2Projection<double, double, ConvDiffTraits::VelocityGradientYType,
                     DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType>()(velocityGradientY,
                                                                                       discrete_velocityGradientY);

  ConvDiffTraits::VelocityGradientXType velocityGradientX(timeprovider_, continousVelocitySpace);
  DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType discrete_velocityGradientX(
      "velocityGradientX", currentFunctions.discreteVelocity().space());
  Dune::L2Projection<double, double, ConvDiffTraits::VelocityGradientXType,
                     DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType>()(velocityGradientX,
                                                                                       discrete_velocityGradientX);

  Dune::BetterL2Projection::project(0.0, velocityGradient, discrete_velocityGradient);

  GradientSplitterFunctionType exact_gradient_splitter(functionSpaceWrapper.discreteVelocitySpace(),
                                                       discrete_velocityGradient);

  ConvDiffTraits::VelocityLaplaceType velocityLaplace(timeprovider_, continousVelocitySpace);
  DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType discrete_velocityLaplace(
      "exact_laplace", currentFunctions.discreteVelocity().space());
  Dune::L2Projection<double, double, ConvDiffTraits::VelocityLaplaceType,
                     DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType>()(velocityLaplace,
                                                                                       discrete_velocityLaplace);

  typedef Stuff::L2Error<GridPartType> L2ErrorType;
  L2ErrorType l2Error(gridPart);
  L2ErrorType::Errors errors_convection =
      l2Error.get(rhs_container.convection, discrete_exactConvection, convection_diff);
  L2ErrorType::Errors errors_velocity = l2Error.get(nextFunctions.discreteVelocity(), exactSolution.discreteVelocity(),
                                                    errorFunctions.discreteVelocity());

  L2ErrorType::Errors errors_gradient = l2Error.get(rhs_container.velocity_gradient, discrete_velocityGradient);
  L2ErrorType::Errors errors_laplace = l2Error.get(rhs_container.velocity_laplace, discrete_velocityLaplace);
  double GD = Stuff::boundaryIntegral(stokesDirichletData, nextFunctions.discreteVelocity().space());

  Logger().Info().Resume();
  Logger().Info()
      //			<< "L2-Error Pressure (abs|rel): " << std::setw(8) << l2_error_pressure << " | " <<
      //relative_l2_error_pressure << "\n"
      << errors_convection.str() << errors_velocity.str() << errors_gradient.str() << errors_laplace.str()
      << "current time: " << timeprovider_.time()
                             //					<< "Mean pressure (exact|discrete): " << meanPressure_exact << " | " <<
                             //meanPressure_discrete << std::endl
      << "\nGD: " << GD << std::endl;

  typedef Stuff::FullTuple<const DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType*> OutputTupleType;
  typedef Dune::TimeAwareDataWriter<ConvDiffTraits::TimeProviderType, GridPartType::GridType, OutputTupleType>
  DataWriterType;
  OutputTupleType out(&nextFunctions.discreteVelocity(), &exactSolution.discreteVelocity(),
                      &errorFunctions.discreteVelocity(), &rhs_container.convection, &discrete_exactConvection,
                      gradient_splitter[0].get(), gradient_splitter[1].get(), exact_gradient_splitter[0].get(),
                      exact_gradient_splitter[1].get());

  DataWriterType dt(timeprovider_, gridPart.grid(), out);
  dt.write();

  OutputTupleType out2(&discrete_velocityLaplace, &rhs_container.velocity_laplace, 0, 0, 0, 0, 0, 0, 0);

  DataWriterType dt2(timeprovider_, gridPart.grid(), out2);
  dt2.write();

  return runInfoVector;
}

void eocCheck(const Stuff::RunInfoVector& runInfos) {
  bool ups = false;
  Stuff::RunInfoVector::const_iterator it = runInfos.begin();
  Stuff::RunInfo last = *it;
  ++it;
  for (; it != runInfos.end(); ++it) {
    ups = (last.L2Errors[0] < it->L2Errors[0] || last.L2Errors[1] < it->L2Errors[1]);
    last = *it;
  }
  if (ups) {
    Logger().Err() << "----------------------------------------------------------\n"
                   << "-                                                        -\n"
                   << "-                  negative EOC                          -\n"
                   << "-                                                        -\n"
                   << "----------------------------------------------------------\n" << std::endl;
  }
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
