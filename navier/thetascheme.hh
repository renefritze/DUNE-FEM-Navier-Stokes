#ifndef METADATA_HH
#define METADATA_HH

#include <dune/navier/thetascheme_base.hh>

namespace Dune {
namespace NavierStokes {
template <class Traits>
class ThetaScheme : public ThetaSchemeBase<Traits> {
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
  ThetaScheme(typename Traits::GridPartType gridPart, const typename Traits::ThetaSchemeDescriptionType& scheme_params,
              typename BaseType::CommunicatorType comm = typename BaseType::CommunicatorType())
    : BaseType(gridPart, scheme_params, comm) {}

  virtual Stuff::RunInfo full_timestep() {
    Stuff::Profiler::ScopedTiming fullstep_time("full_step");
    Stuff::RunInfo info;
    for (int i = 0; i < Traits::substep_count; ++i) {
      const double dt_k = scheme_params_.step_sizes_[i];
      substep(dt_k, scheme_params_.thetas_[i]);
      if (i != Traits::substep_count - 1)
        // the last step increase is done after one call level up
        BaseType::nextStep(i, info);
    }
    return info;
  }

  boost::shared_ptr<typename Traits::OseenForceAdapterFunctionType>
  prepare_rhs(const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& theta_values) {
    const bool first_step = timeprovider_.timeStep() <= 2;
    const typename Traits::AnalyticalForceType force(timeprovider_, currentFunctions_.discreteVelocity().space(),
                                                     viscosity_, 0.0 /*stokes alpha*/);
    const bool do_cheat = Parameters().getParam("rhs_cheat", false);

    if (!Parameters().getParam("parabolic", false) &&
        (scheme_params_.algo_id == Traits::ThetaSchemeDescriptionType::scheme_names[3] /*CN*/)) {
      // reconstruct the prev convection term
      auto beta = currentFunctions_.discreteVelocity();
      beta *= 1.5;
      auto dummy = lastFunctions_.discreteVelocity();
      dummy *= 0.5;
      beta -= dummy;
      Dune::BruteForceReconstruction<typename Traits::OseenModelType>::getConvection(
          beta, rhsDatacontainer_.velocity_gradient, rhsDatacontainer_.convection);
    }
    auto null_f(exactSolution_.discreteVelocity());
    null_f *= 0.0;
    boost::shared_ptr<typename Traits::OseenForceAdapterFunctionType> ptr_oseenForceVanilla(
        first_step // in our very first step no previous computed data is avail. in rhs_container
            ? new typename Traits::OseenForceAdapterFunctionType(timeprovider_, exactSolution_.discreteVelocity(),
                                                                 force, reynolds_, theta_values)
            : new typename Traits::OseenForceAdapterFunctionType(timeprovider_, currentFunctions_.discreteVelocity(),
                                                                 force, reynolds_, theta_values, rhsDatacontainer_));
    //					if ( do_cheat )
    BaseType::cheatRHS();
    boost::shared_ptr<typename Traits::OseenForceAdapterFunctionType> ptr_oseenForce(
        first_step // in our very first step no previous computed data is avail. in rhs_container
            ? new typename Traits::OseenForceAdapterFunctionType(timeprovider_, exactSolution_.discreteVelocity(),
                                                                 force, reynolds_, theta_values)
            : new typename Traits::OseenForceAdapterFunctionType(timeprovider_, currentFunctions_.discreteVelocity(),
                                                                 force, reynolds_, theta_values, rhsDatacontainer_));
    typename BaseType::L2ErrorType::Errors errors_rhs =
        l2Error_.get(static_cast<typename Traits::StokesForceAdapterType::BaseType>(*ptr_oseenForce),
                     static_cast<typename Traits::StokesForceAdapterType::BaseType>(*ptr_oseenForceVanilla));
    Logger().Dbg().Resume(9000);
    Logger().Dbg() << "RHS " << errors_rhs.str() << std::endl;

    rhsFunctions_.discreteVelocity().assign(*ptr_oseenForce);
    return do_cheat ? ptr_oseenForce : ptr_oseenForceVanilla;
  }

  typename Traits::OseenPassType
  prepare_pass(const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& theta_values) {
    const auto rhs = prepare_rhs(theta_values);
    const double dt_n = timeprovider_.deltaT();
    const bool do_convection_disc =
        !(Parameters().getParam("navier_no_convection", false) || Parameters().getParam("parabolic", false));
    auto beta = currentFunctions_.discreteVelocity(); //=u^n = bwe linearization
    if (do_convection_disc && (scheme_params_.algo_id == Traits::ThetaSchemeDescriptionType::scheme_names[3] /*CN*/)) {
      // linearization: 1.5u^n-0.5u^{n-1}
      beta *= 1.5;
      auto dummy = lastFunctions_.discreteVelocity();
      dummy *= 0.5;
      beta -= dummy;
    } else if (!do_convection_disc)
      beta.clear();
    Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
    stab_coeff.FactorFromParams("C11");
    stab_coeff.FactorFromParams("C12");
    stab_coeff.FactorFromParams("D11");
    stab_coeff.FactorFromParams("D12");
    typename Traits::AnalyticalDirichletDataType oseenDirichletData(timeprovider_, functionSpaceWrapper_,
                                                                    theta_values[0], 1 - theta_values[0]);
    typename Traits::OseenModelType oseenModel(
        stab_coeff, *rhs, oseenDirichletData, theta_values[0] / reynolds_, /*viscosity*/
        1.0f / dt_n,                                                       /*alpha*/
        //											do_convection_disc ? theta_values[0] * dt_n : 0.0, /*convection_scale_factor*/
        theta_values[0], /*convection_scale_factor*/
        theta_values[0]  /*pressure_gradient_scale_factor*/
        );

    return typename Traits::OseenPassType(oseenModel, gridPart_, functionSpaceWrapper_, beta /*beta*/,
                                          do_convection_disc /*do_oseen_disc*/);
  }

  void substep(const double /*dt_k*/,
               const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& theta_values) {
    {
      typename Traits::AnalyticalDirichletDataType oseenDirichletData(timeprovider_, functionSpaceWrapper_);
      Dune::BetterL2Projection::project(timeprovider_, oseenDirichletData, dummyFunctions_.discreteVelocity());
      const double boundaryInt =
          Stuff::boundaryIntegral(oseenDirichletData, BaseType::currentFunctions().discreteVelocity().space());
      Logger().Dbg() << boost::format("discrete Boundary integral: %e\n") % boundaryInt;
    }
    auto oseenPass = prepare_pass(theta_values);
    if (timeprovider_.timeStep() <= 2)
      oseenPass.printInfo();
    if (Parameters().getParam("silent_stokes", true))
      Logger().Info().Suspend(Stuff::Logging::LogStream::default_suspend_priority + 10);
    if (Parameters().getParam("clear_u", false))
      nextFunctions_.discreteVelocity().clear();
    if (Parameters().getParam("clear_p", false))
      nextFunctions_.discretePressure().clear();
    oseenPass.apply(currentFunctions_, nextFunctions_, &rhsDatacontainer_);
    Logger().Info().Resume(Stuff::Logging::LogStream::default_suspend_priority + 10);
    BaseType::setUpdateFunctions();
    currentFunctions_.assign(nextFunctions_);
  }
};
} // end namespace NavierStokes
} // end namespace Dune

#endif // METADATA_HH

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
