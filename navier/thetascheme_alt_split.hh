#ifndef THETASCHEME_ALT_SPLIT_HH
#define THETASCHEME_ALT_SPLIT_HH

#include <dune/navier/thetascheme_base.hh>

namespace Dune {
namespace NavierStokes {
template <class Traits>
class ThetaSchemeAltSplitting : public ThetaSchemeBase<Traits> {
protected:
  typedef ThetaSchemeBase<Traits> BaseType;
  using BaseType::gridPart_;
  using BaseType::scheme_params_;

protected:
  using BaseType::timeprovider_;
  using BaseType::functionSpaceWrapper_;
  using BaseType::currentFunctions_;
  using BaseType::nextFunctions_;
  using BaseType::exactSolution_;
  using BaseType::dummyFunctions_;
  using BaseType::rhsFunctions_;
  using BaseType::rhsDatacontainer_;
  using BaseType::lastFunctions_;
  using BaseType::l2Error_;

public:
  using BaseType::viscosity_;
  using BaseType::reynolds_;

public:
  ThetaSchemeAltSplitting(typename Traits::GridPartType gridPart,
                          const typename Traits::ThetaSchemeDescriptionType& scheme_params,
                          typename BaseType::CommunicatorType comm = typename BaseType::CommunicatorType())
    : BaseType(gridPart, scheme_params, comm) {}

  virtual Stuff::RunInfo full_timestep() {
    Stuff::RunInfo info;
    {
      Stuff::Profiler::ScopedTiming fullstep_time("full_step");
      Stuff::RunInfo info_dummy;
      // stokes step A
      typename BaseType::DiscreteVelocityFunctionType u_n("u_n", dummyFunctions_.discreteVelocity().space());
      u_n.assign(currentFunctions_.discreteVelocity());
      stokesStep(scheme_params_.step_sizes_[0], scheme_params_.thetas_[0]);
      BaseType::nextStep(0, info_dummy);

      Parameters().setParam("reduced_oseen_solver", true);
      // Nonlinear step
      nonlinearStep(scheme_params_.step_sizes_[1], scheme_params_.thetas_[1], u_n);
      BaseType::nextStep(1, info_dummy);
      Parameters().setParam("reduced_oseen_solver", false);

      // stokes step B
      info = stokesStep(scheme_params_.step_sizes_[2], scheme_params_.thetas_[2]);
    }
    BaseType::nextStep(2, info);

    return info;
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

  Stuff::RunInfo
  stokesStep(const double /*dt_k*/,
             const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& /*theta_values*/) const {
    DiscretizationWeights discretization_weights(BaseType::d_t_, viscosity_);

    if (Parameters().getParam("silent_stokes", true))
      Logger().Suspend(Stuff::Logging::LogStream::default_suspend_priority + 1);

    const bool first_stokes_step = timeprovider_.timeStep() <= 1;
    const typename Traits::AnalyticalForceType force(timeprovider_, currentFunctions_.discreteVelocity().space(),
                                                     viscosity_, 0.0 /*stokes alpha*/);

    boost::scoped_ptr<typename Traits::StokesForceAdapterType> ptr_stokesForce_vanilla(
        first_stokes_step
            ? new typename Traits::StokesForceAdapterType(timeprovider_, currentFunctions_.discreteVelocity(), force,
                                                          discretization_weights)
            : new typename Traits::StokesForceAdapterType(timeprovider_, currentFunctions_.discreteVelocity(), force,
                                                          discretization_weights, rhsDatacontainer_));

    typedef Stuff::L2Error<typename Traits::GridPartType> L2ErrorType;
    L2ErrorType l2Error(gridPart_);

    // CHEAT (projecting the anaylitcal evals into the container filled by last pass
    const bool do_cheat = Parameters().getParam("rhs_cheat", false) && !first_stokes_step;
    dummyFunctions_.discreteVelocity().assign(currentFunctions_.discreteVelocity());
    //					if ( do_cheat ) //do cheat rhs assembly unconditionally, below we'll choose according to do_cheat which
    //rhs to put into the model
    {
      typedef typename BaseType::DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType
      VelocityFunctionSpaceType;
      VelocityFunctionSpaceType continousVelocitySpace_;
      typedef NAVIER_DATA_NAMESPACE::VelocityConvection<VelocityFunctionSpaceType, typename Traits::TimeProviderType>
      VelocityConvection;
      VelocityConvection velocity_convection(timeprovider_, continousVelocitySpace_);
      Dune::BetterL2Projection // we need evals from the _previous_ (t_0) step
          ::project(timeprovider_.previousSubTime(), velocity_convection, rhsDatacontainer_.convection);
      //						// ----
      typedef NAVIER_DATA_NAMESPACE::VelocityLaplace<VelocityFunctionSpaceType, typename Traits::TimeProviderType>
      VelocityLaplace;
      VelocityLaplace velocity_laplace(timeprovider_, continousVelocitySpace_);
      Dune::BetterL2Projection // this seems currently inconsequential to the produced error
          ::project(timeprovider_.previousSubTime(), velocity_laplace, rhsDatacontainer_.velocity_laplace);

      //						typename L2ErrorType::Errors errors_convection = l2Error.get(	exactSolution_.discreteVelocity() ,
      //																			currentFunctions_.discreteVelocity(),
      //																			dummyFunctions_.discreteVelocity() );
      //						std::cerr << "BLAH " << errors_convection.str();

      currentFunctions_.discreteVelocity().assign(exactSolution_.discreteVelocity());
    } // END CHEAT

    boost::scoped_ptr<typename Traits::StokesForceAdapterType> ptr_stokesForce(
        first_stokes_step
            ? new typename Traits::StokesForceAdapterType(timeprovider_, currentFunctions_.discreteVelocity(), force,
                                                          discretization_weights)
            : new typename Traits::StokesForceAdapterType(timeprovider_, currentFunctions_.discreteVelocity(), force,
                                                          discretization_weights, rhsDatacontainer_));
    typename L2ErrorType::Errors errors_rhs =
        l2Error.get(static_cast<typename Traits::StokesForceAdapterType::BaseType>(*ptr_stokesForce),
                    static_cast<typename Traits::StokesForceAdapterType::BaseType>(*ptr_stokesForce_vanilla),
                    dummyFunctions_.discreteVelocity());
    std::cerr << "RHS " << errors_rhs.str();

    Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();

    {
      //					if ( Parameters().getParam( "stab_coeff_visc_scale", true ) ) {
      //						stab_coeff.Factor( "D11", ( 1 / stokes_viscosity ) );
      //						stab_coeff.Factor( "C11", stokes_viscosity );
      //					}
      //					else {
      //						stab_coeff.FactorFromParams("D11");
      //						stab_coeff.FactorFromParams("C11");
      //					}
      //					stab_coeff.FactorFromParams("D12");
      //					stab_coeff.FactorFromParams("C12");
      //					stab_coeff.Add( "E12", 0.5 );
      //					stab_coeff.print( Logger().Info() );
    }

    typename Traits::AnalyticalDirichletDataType stokesDirichletData(timeprovider_, functionSpaceWrapper_);

    typename Traits::StokesModelType stokesModel(
        stab_coeff, do_cheat ? *ptr_stokesForce : *ptr_stokesForce_vanilla, stokesDirichletData,
        discretization_weights.alpha * discretization_weights.theta_times_delta_t, /*viscosity*/
        1.0,                                                                       /*alpha*/
        0.0,                                                                       /*convection_scale_factor*/
        discretization_weights.theta_times_delta_t /*pressure_gradient_scale_factor*/);
    typename Traits::StokesPassType stokesPass(stokesModel, gridPart_, functionSpaceWrapper_,
                                               dummyFunctions_.discreteVelocity(), false);

    stokesPass.apply(currentFunctions_, nextFunctions_, &rhsDatacontainer_);
    BaseType::setUpdateFunctions();
    Stuff::RunInfo info;
    stokesPass.getRuninfo(info);
    if (Parameters().getParam("silent_stokes", true))
      Logger().Resume(Stuff::Logging::LogStream::default_suspend_priority + 1);
    return info;
  }

  void nonlinearStep(const double /*dt_k*/,
                     const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& /*theta_values*/,
                     const typename BaseType::DiscreteVelocityFunctionType& u_n) {
    DiscretizationWeights discretization_weights(BaseType::d_t_, viscosity_);

    const typename Traits::AnalyticalForceType force(timeprovider_, currentFunctions_.discreteVelocity().space(),
                                                     viscosity_, 0.0 /*stokes alpha*/);

    // CHEAT (projecting the anaylitcal evals into the container filled by last pass
    if (Parameters().getParam("rhs_cheat", false)) {
      typedef typename BaseType::DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType
      VelocityFunctionSpaceType;
      VelocityFunctionSpaceType continousVelocitySpace_;

      typedef NAVIER_DATA_NAMESPACE::PressureGradient<VelocityFunctionSpaceType, typename Traits::TimeProviderType>
      PressureGradient;
      PressureGradient pressure_gradient(timeprovider_, continousVelocitySpace_);
      Dune::BetterL2Projection // we need evals from the _previous_ (t_0) step
          ::project(timeprovider_.previousSubTime(), pressure_gradient, rhsDatacontainer_.pressure_gradient);
      // ----
      typedef NAVIER_DATA_NAMESPACE::VelocityLaplace<VelocityFunctionSpaceType, typename Traits::TimeProviderType>
      VelocityLaplace;
      VelocityLaplace velocity_laplace(timeprovider_, continousVelocitySpace_);
      Dune::BetterL2Projection::project(timeprovider_.previousSubTime(), velocity_laplace,
                                        rhsDatacontainer_.velocity_laplace);
      currentFunctions_.discreteVelocity().assign(exactSolution_.discreteVelocity());
    } // END CHEAT

    typename Traits::NonlinearForceAdapterType nonlinearForce(timeprovider_, currentFunctions_.discreteVelocity(),
                                                              force, discretization_weights, rhsDatacontainer_);

    rhsFunctions_.discreteVelocity().assign(nonlinearForce);
    unsigned int oseen_iterations = Parameters().getParam("oseen_iterations", (unsigned int)(1));
    assert(oseen_iterations > 0);
    nonlinearStepSingle(nonlinearForce, discretization_weights, u_n);
  }

  template <class T>
  void nonlinearStepSingle(const T& nonlinearForce, const DiscretizationWeights& discretization_weights,
                           const typename BaseType::DiscreteVelocityFunctionType& u_n) {
    typename Traits::StokesStartPassType stokesStartPass;

    typename Traits::AnalyticalDirichletDataType stokesDirichletData(timeprovider_, functionSpaceWrapper_);
    Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
    //					if ( Parameters().getParam( "stab_coeff_visc_scale", true ) ) {
    //						stab_coeff.Factor( "D11", ( 1 / oseen_viscosity )  );
    //						stab_coeff.Factor( "C11", oseen_viscosity );
    //					}
    //					else {
    //						stab_coeff.FactorFromParams("D11");
    //						stab_coeff.FactorFromParams("C11");
    //					}
    //					stab_coeff.FactorFromParams("D12");
    //					stab_coeff.FactorFromParams("C12");
    //					stab_coeff.Add( "E12", 0.5 );

    //					stab_coeff.print( Logger().Info() );

    typename BaseType::DiscreteVelocityFunctionType beta("beta", dummyFunctions_.discreteVelocity().space());
    typename BaseType::DiscreteVelocityFunctionType tmp("tmp", dummyFunctions_.discreteVelocity().space());
    tmp.assign(u_n);
    const double theta = discretization_weights.theta;
    tmp *= (2.0 * theta) / (1.0 - theta);
    beta.assign(currentFunctions_.discreteVelocity());
    beta *= theta / (1.0 - theta);
    beta += tmp;

    typename Traits::NonlinearModelType stokesModel(
        stab_coeff, nonlinearForce, stokesDirichletData,
        discretization_weights.beta * discretization_weights.one_neg_two_theta_dt, /*viscosity*/
        1.0,                                                                       /*alpha*/
        discretization_weights.one_neg_two_theta_dt,                               /*convection_scale_factor*/
        0.0 /*pressure_gradient_scale_factor*/);
    typename Traits::NonlinearPassType oseenPass(stokesModel, gridPart_, functionSpaceWrapper_, beta, true);
    oseenPass.apply(currentFunctions_, nextFunctions_, &rhsDatacontainer_);
  }
};
} // end namespace NavierStokes
} // end namespace Dune

#endif // THETASCHEME_ALT_SPLIT_HH

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
