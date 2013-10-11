#ifndef THETASCHEMEBASE_H
#define THETASCHEMEBASE_H

#include <dune/navier/exactsolution.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/stuff/tuple.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/stuff/error.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/profiler.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>
#include <boost/scoped_ptr.hpp>
#include <algorithm>

#include <dune/navier/thetascheme_traits.hh>
#include <dune/navier/fractionaldatawriter.hh>
#include <dune/navier/global_defines.hh>
#include <dune/oseen/pass.hh>

namespace Dune {
namespace NavierStokes {
template <class TraitsImp>
class ThetaSchemeBase {
protected:
  typedef TraitsImp Traits;
  typedef typename Traits::CommunicatorType CommunicatorType;
  typedef typename Traits::ExactSolutionType ExactSolutionType;
  typedef Stuff::TupleSerializer<typename Traits::DiscreteOseenFunctionWrapperType,
                                 typename Traits::DiscreteOseenFunctionWrapperType, ExactSolutionType,
                                 typename Traits::DiscreteOseenFunctionWrapperType> TupleSerializerType1;
  typedef typename TupleSerializerType1::TupleType OutputTupleType1;
  typedef TimeAwareDataWriter<typename Traits::TimeProviderType, typename Traits::GridPartType::GridType,
                              OutputTupleType1> DataWriterType1;
  typedef CheckPointer<typename Traits::GridPartType::GridType, OutputTupleType1> CheckPointerType;
  typedef Stuff::TupleSerializer<typename Traits::DiscreteOseenFunctionWrapperType> TupleSerializerType2;
  typedef typename TupleSerializerType2::TupleType OutputTupleType2;
  typedef TimeAwareDataWriter<typename Traits::TimeProviderType, typename Traits::GridPartType::GridType,
                              OutputTupleType2> DataWriterType2;
  typedef Dune::Oseen::RhsDatacontainer<typename Traits::OseenModelTraits> DataContainerType;
  typedef typename Traits::DiscreteOseenFunctionWrapperType::DiscreteVelocityFunctionType DiscreteVelocityFunctionType;
  typedef typename Traits::DiscreteOseenFunctionWrapperType::DiscretePressureFunctionType DiscretePressureFunctionType;

  mutable typename Traits::GridPartType gridPart_;
  const typename Traits::ThetaSchemeDescriptionType& scheme_params_;

protected:
  CommunicatorType communicator_;
  typename Traits::TimeProviderType timeprovider_;
  typename Traits::DiscreteOseenFunctionSpaceWrapperType functionSpaceWrapper_;
  mutable typename Traits::DiscreteOseenFunctionWrapperType currentFunctions_;
  mutable typename Traits::DiscreteOseenFunctionWrapperType nextFunctions_;
  typename Traits::DiscreteOseenFunctionWrapperType errorFunctions_;
  ExactSolutionType exactSolution_;
  mutable typename Traits::DiscreteOseenFunctionWrapperType dummyFunctions_;
  mutable typename Traits::DiscreteOseenFunctionWrapperType updateFunctions_;
  mutable typename Traits::DiscreteOseenFunctionWrapperType rhsFunctions_;
  OutputTupleType1& data_tuple_1;
  DataWriterType1 dataWriter1_;
  CheckPointerType check_pointer_;
  DataWriterType2 dataWriter2_;
  const typename Traits::OseenPassType::Traits::DiscreteSigmaFunctionSpaceType sigma_space_;
  mutable DataContainerType rhsDatacontainer_;
  mutable typename Traits::DiscreteOseenFunctionWrapperType lastFunctions_;

  typedef Stuff::L2Error<typename Traits::GridPartType> L2ErrorType;
  L2ErrorType l2Error_;

public:
  const double viscosity_;
  const double d_t_;
  const double reynolds_;
  double current_max_gridwidth_;

public:
  virtual ~ThetaSchemeBase() {}

  ThetaSchemeBase(typename Traits::GridPartType gridPart,
                  const typename Traits::ThetaSchemeDescriptionType& scheme_params,
                  CommunicatorType comm = Dune::MPIManager::helper().getCommunicator())
    : gridPart_(gridPart)
    , scheme_params_(scheme_params)
    , communicator_(comm)
    , timeprovider_(scheme_params_, communicator_)
    , functionSpaceWrapper_(gridPart_)
    , currentFunctions_("current_", functionSpaceWrapper_, gridPart_)
    , nextFunctions_("next_", functionSpaceWrapper_, gridPart_)
    , errorFunctions_("error_", functionSpaceWrapper_, gridPart_)
    , exactSolution_(timeprovider_, gridPart_, functionSpaceWrapper_)
    , dummyFunctions_("dummy", functionSpaceWrapper_, gridPart_)
    , updateFunctions_("updates", functionSpaceWrapper_, gridPart_)
    , rhsFunctions_("rhs-adapter", functionSpaceWrapper_, gridPart_)
    , data_tuple_1(TupleSerializerType1::getTuple(currentFunctions_, errorFunctions_, exactSolution_, dummyFunctions_))
    , dataWriter1_(timeprovider_, gridPart_.grid(), data_tuple_1)
    , check_pointer_(gridPart_.grid(), data_tuple_1, timeprovider_)
    , dataWriter2_(timeprovider_, gridPart_.grid(), TupleSerializerType2::getTuple(updateFunctions_, rhsFunctions_))
    , sigma_space_(gridPart_)
    , rhsDatacontainer_(currentFunctions_.discreteVelocity().space(), sigma_space_)
    , lastFunctions_("last", functionSpaceWrapper_, gridPart_)
    , l2Error_(gridPart)
    , viscosity_(Parameters().getParam("viscosity", 1.0, Dune::ValidateNotLess<double>(0.0)))
    , d_t_(timeprovider_.deltaT())
    , reynolds_(1.0 / viscosity_)
    , current_max_gridwidth_(Dune::GridWidth::calcGridWidth(gridPart_)) {
    Logger().Info() << scheme_params_;
    NAVIER_DATA_NAMESPACE::SetupCheck check;
    if (!check(this, gridPart_, scheme_params_, timeprovider_, functionSpaceWrapper_))
      DUNE_THROW(InvalidStateException, check.error());
  }

  void nextStep(const int step, Stuff::RunInfo& info) {
    current_max_gridwidth_ = Dune::GridWidth::calcGridWidth(gridPart_);
    lastFunctions_.assign(currentFunctions_);
    currentFunctions_.assign(nextFunctions_);
    exactSolution_.project();
    const bool last_substep = (step == (Traits::ThetaSchemeDescriptionType::numberOfSteps_ - 1));

    // error calc
    if (NAVIER_DATA_NAMESPACE::hasExactSolution && Parameters().getParam("calculate_errors", true)) {
      Stuff::Profiler::ScopedTiming error_time("error_calc");

      errorFunctions_.discretePressure().assign(exactSolution_.discretePressure());
      errorFunctions_.discretePressure() -= currentFunctions_.discretePressure();
      errorFunctions_.discreteVelocity().assign(exactSolution_.discreteVelocity());
      errorFunctions_.discreteVelocity() -= currentFunctions_.discreteVelocity();

      double meanPressure_exact =
          Stuff::integralAndVolume(exactSolution_.exactPressure(), currentFunctions_.discretePressure().space()).first;
      double meanPressure_discrete = Stuff::integralAndVolume(currentFunctions_.discretePressure(),
                                                              currentFunctions_.discretePressure().space()).first;

      Dune::L2Norm<typename Traits::GridPartType> l2_Error(gridPart_);
      Dune::H1Norm<typename Traits::GridPartType> h1_Error(gridPart_);

      //						if ( Parameters().getParam( "error_scaling", false ) ) {
      //								const double scale		= 1 / std::sqrt( viscosity_ );
      //								errorFunctions_.discretePressure() *= scale;
      //								errorFunctions_.discreteVelocity() *= scale;
      //						}

      const double l2_error_pressure_ = l2_Error.norm(errorFunctions_.discretePressure());
      const double l2_error_velocity_ = l2_Error.norm(errorFunctions_.discreteVelocity());
      const double h1_error_pressure_ = h1_Error.norm(errorFunctions_.discretePressure());
      const double h1_error_velocity_ = h1_Error.norm(errorFunctions_.discreteVelocity());
      const double relative_l2_error_pressure_ = l2_error_pressure_ / l2_Error.norm(exactSolution_.discretePressure());
      const double relative_l2_error_velocity_ = l2_error_velocity_ / l2_Error.norm(exactSolution_.discreteVelocity());
      const double relative_h1_error_velocity_ = h1_error_velocity_ / h1_Error.norm(exactSolution_.discreteVelocity());
      std::vector<double> error_vector;
      error_vector.push_back(l2_error_velocity_);
      error_vector.push_back(l2_error_pressure_);
      std::vector<double> h1_error_vector;
      h1_error_vector.push_back(h1_error_velocity_);
      h1_error_vector.push_back(h1_error_pressure_);

#ifdef NDEBUG
      if (last_substep) // no need to be so verbose otherwise
#endif
      {
        Logger().Info().Resume();
        if (Parameters().getParam("parabolic", false))
          Logger().Info() << boost::format("L2-Error Velocity (abs|rel): %e | %e") % l2_error_velocity_ %
                                 relative_l2_error_velocity_;
        else
          Logger().Info() << boost::format("L2-Error Pressure (abs|rel): %e | %e \t Velocity (abs|rel): %e | %e") %
                                 l2_error_pressure_ % relative_l2_error_pressure_ % l2_error_velocity_ %
                                 relative_l2_error_velocity_;
#ifndef NDEBUG
        Logger().Info() << boost::format(
                               "\nH1-Error Velocity (abs|rel): %e | %e \t Mean pressure (exact|discrete) %e | %e") %
                               h1_error_velocity_ % relative_h1_error_velocity_ % meanPressure_exact %
                               meanPressure_discrete;
#endif
        Logger().Info() << std::endl;
      }
      const double max_l2_error = Parameters().getParam("max_error", 1e2, Dune::ValidateGreater<double>(0.0));
      info.L2Errors = error_vector;
      info.H1Errors = h1_error_vector;
      if (l2_error_velocity_ > max_l2_error ||
          (!Parameters().getParam("parabolic", false) && l2_error_pressure_ > max_l2_error))
        throw Stuff::singlerun_abort_exception("Aborted, L2 error above " + Stuff::toString(max_l2_error));
      if (!Parameters().getParam("parabolic", false) &&
          (std::isnan(l2_error_velocity_) || std::isnan(l2_error_pressure_)))
        throw Stuff::singlerun_abort_exception("L2 error is Nan");
    }
    // end error calc

    if (last_substep) {
      typedef Dune::StabilizationCoefficients::ValueType Pair;
      Dune::StabilizationCoefficients stabil_coeff =
          Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();

      info.codim0 = gridPart_.grid().size(0);
      info.grid_width = current_max_gridwidth_;
      info.run_time = profiler().GetTiming("full_step");
      info.delta_t = timeprovider_.deltaT();
      info.current_time = timeprovider_.subTime();
      info.viscosity = viscosity_;
      info.reynolds = reynolds_;

      info.c11 = Pair(stabil_coeff.Power("C11"), stabil_coeff.Factor("C11"));
      info.c12 = Pair(stabil_coeff.Power("C12"), stabil_coeff.Factor("C12"));
      info.d11 = Pair(stabil_coeff.Power("D11"), stabil_coeff.Factor("D11"));
      info.d12 = Pair(stabil_coeff.Power("D12"), stabil_coeff.Factor("D12"));
      info.bfg = Parameters().getParam("do-bfg", true);
      // TODO gridname
      //						info.gridname		= gridPart_.grid().name();
      info.refine_level = Parameters().getParam("minref", 0, Dune::ValidateNotLess<int>(0));

      info.polorder_pressure = Traits::OseenModelTraits::pressureSpaceOrder;
      info.polorder_sigma = Traits::OseenModelTraits::sigmaSpaceOrder;
      info.polorder_velocity = Traits::OseenModelTraits::velocitySpaceOrder;

      info.solver_accuracy = Parameters().getParam("absLimit", 1e-4);
      info.inner_solver_accuracy = Parameters().getParam("inner_absLimit", 1e-4);
      info.bfg_tau = Parameters().getParam("bfg-tau", 0.1);

      info.problemIdentifier = NAVIER_DATA_NAMESPACE::identifier;
      info.algo_id = scheme_params_.algo_id;
      info.extra_info = (boost::format("%s on %s") % COMMIT % std::getenv("HOSTNAME")).str();

      Logger().Info() << boost::format("current time (substep %d ): %f (%f)\n") % step % timeprovider_.subTime() %
                             timeprovider_.previousSubTime();
    }

    if (last_substep || !Parameters().getParam("write_fulltimestep_only", false))
      writeData();
    timeprovider_.nextFractional();
  }

  void Init() {
    typename Traits::TimeProviderType::StepZeroGuard step0(timeprovider_.stepZeroGuard(d_t_));
    // initial flow field at t = 0
    exactSolution_.project();
    currentFunctions_.assign(exactSolution_);
    nextFunctions_.assign(exactSolution_);
    writeData();
    // the guard dtor sets current time to t_0 + dt_k
  }

  Stuff::RunInfoTimeMap run() {
    Stuff::RunInfoTimeMap runInfoMap;
    Init();

    for (; timeprovider_.time() <= timeprovider_.endTime();) {
      assert(timeprovider_.time() > 0.0);
      Stuff::RunInfo info = full_timestep();
      const double real_time = timeprovider_.subTime();
      try {
        nextStep(Traits::substep_count - 1, info);
      }
      catch (Stuff::singlerun_abort_exception& e) {
        Logger().Err() << e.what() << std::endl;
        // fill up the map with dummy data so it can still be used in output
        runInfoMap[real_time] = info;
        for (; timeprovider_.time() <= timeprovider_.endTime();) {
          timeprovider_.nextFractional();
          runInfoMap[timeprovider_.subTime()] = Stuff::RunInfo::dummy();
        }
        return runInfoMap;
      }
      timeprovider_.printRemainderEstimate(Logger().Info());
      runInfoMap[real_time] = info;
    }
    assert(runInfoMap.size() > 0);
    return runInfoMap;
  }

  virtual Stuff::RunInfo full_timestep() = 0;

  void setUpdateFunctions() const {
    updateFunctions_.assign(nextFunctions_);
    updateFunctions_ -= currentFunctions_;
  }

  void writeData() {
    Stuff::Profiler::ScopedTiming io_time("IO");
    dataWriter1_.write();
    dataWriter2_.write();
    //					check_pointer_.write( timeprovider_.time(), timeprovider_.timeStep() );
  }

  DataContainerType& rhsDatacontainer() { return rhsDatacontainer_; }

  const ExactSolutionType& exactSolution() const { return exactSolution_; }

  const typename Traits::DiscreteOseenFunctionWrapperType& currentFunctions() const { return currentFunctions_; }

  const typename Traits::TimeProviderType& timeprovider() const { return timeprovider_; }

  void cheatRHS() {
    typedef typename DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType VelocityFunctionSpaceType;
    VelocityFunctionSpaceType continousVelocitySpace_;

    // ----
    typedef NAVIER_DATA_NAMESPACE::VelocityLaplace<VelocityFunctionSpaceType, typename Traits::TimeProviderType>
    VelocityLaplace;
    VelocityLaplace velocity_laplace(timeprovider_, continousVelocitySpace_);
    Dune::BetterL2Projection::project(timeprovider_.previousSubTime(), velocity_laplace,
                                      rhsDatacontainer_.velocity_laplace);
    currentFunctions_.discreteVelocity().assign(exactSolution_.discreteVelocity());

    typedef NAVIER_DATA_NAMESPACE::VelocityConvection<VelocityFunctionSpaceType, typename Traits::TimeProviderType>
    VelocityConvection;
    VelocityConvection velocity_convection(timeprovider_, continousVelocitySpace_);
    typedef NAVIER_DATA_NAMESPACE::PressureGradient<VelocityFunctionSpaceType, typename Traits::TimeProviderType>
    PressureGradient;
    PressureGradient pressure_gradient(timeprovider_, continousVelocitySpace_);

    Dune::BetterL2Projection::project(timeprovider_.previousSubTime(), pressure_gradient,
                                      rhsDatacontainer_.pressure_gradient);
    Dune::BetterL2Projection // we need evals from the _previous_ (t_{k-1}) step
        ::project(timeprovider_.previousSubTime(), velocity_convection, rhsDatacontainer_.convection);
  }

  struct DiscretizationWeights {
    const double theta, alpha, beta, theta_times_delta_t, viscosity, one_neg_two_theta_dt;
    DiscretizationWeights(const double d_t, const double visc)
      : theta(1.0 - (std::sqrt(2) / 2.0f))
      , alpha((1.0 - 2 * theta) / (1.0 - theta))
      , beta(1.0 - alpha)
      , theta_times_delta_t(theta * d_t)
      , viscosity(visc)
      , one_neg_two_theta_dt((1. - 2. * theta) * d_t) {}
  };
};
template <class Traits>
class DataOnlyScheme : public ThetaSchemeBase<Traits> {
protected:
  typedef ThetaSchemeBase<Traits> BaseType;
  using BaseType::timeprovider_;

public:
  DataOnlyScheme(typename Traits::GridPartType gridPart,
                 const typename Traits::ThetaSchemeDescriptionType& scheme_params,
                 typename BaseType::CommunicatorType comm = typename BaseType::CommunicatorType())
    : BaseType(gridPart, scheme_params, comm) {}

  virtual Stuff::RunInfo full_timestep() {
    timeprovider_.printRemainderEstimate(Logger().Info());
    return Stuff::RunInfo();
  }
};
} // end namespace NavierStokes
} // end namespace Dune

#endif // THETASCHEMEBASE_H

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
