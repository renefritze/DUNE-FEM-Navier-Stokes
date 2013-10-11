#ifndef THETASCHEME_TRAITS_HH
#define THETASCHEME_TRAITS_HH

#include <dune/oseen/modelinterface.hh>
#include <dune/oseen/pass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/stuff/functions.hh>
#include <dune/stuff/misc.hh>

namespace Dune {
namespace NavierStokes {
namespace {
const std::string scheme_names_array[] = {"N.A.", "FWE", "BWE", "CN", "FS0", "FS1"};
}
//! for each step keep a set of theta values and one value for dt
// thetas_[stepnumber][theta_subscript_index]
template <int numberOfSteps>
struct ThetaSchemeDescription {
  static const std::vector<std::string> scheme_names;
  typedef ThetaSchemeDescription<numberOfSteps> ThisType;
  static const int numberOfSteps_ = numberOfSteps;
  typedef Dune::array<double, 4> ThetaValueArray;
  typedef Dune::array<ThetaValueArray, numberOfSteps> ThetaArray;
  typedef Stuff::wraparound_array<double, numberOfSteps> TimestepArray;
  ThetaArray thetas_;
  TimestepArray step_sizes_;
  std::string algo_id;

private:
  ThetaSchemeDescription() : algo_id("N.A.") {}

  ThetaSchemeDescription(const ThetaArray& a, const TimestepArray& s, std::string id)
    : thetas_(a)
    , step_sizes_(s)
    , algo_id(id) {}

  ThetaSchemeDescription(const ThetaArray& a, double dt, std::string id)
    : thetas_(a)
    , algo_id(id) {
    Stuff::fill_entirely(step_sizes_, dt);
  }

public:
  static ThetaSchemeDescription<1> crank_nicholson(double delta_t) {
    ThetaValueArray c;
    Stuff::fill_entirely(c, 0.5f);
    ThetaArray a;
    Stuff::fill_entirely(a, c);
    return ThisType(a, delta_t, scheme_names[3]);
  }
  static ThetaSchemeDescription<1> forward_euler(double delta_t) {
    // double braces to silence gcc warnigns (throughout this file)
    ThetaValueArray c = {{0.0f, 1.0f, 1.0f, 0.0f}};
    ThetaArray a;
    Stuff::fill_entirely(a, c);
    return ThetaSchemeDescription<1>(a, delta_t, scheme_names[1]);
  }
  static ThetaSchemeDescription<1> backward_euler(double delta_t) {
    ThetaValueArray c = {{1.0f, 0.0f, 0.0f, 1.0f}};
    ThetaArray a;
    Stuff::fill_entirely(a, c);
    return ThetaSchemeDescription<1>(a, delta_t, scheme_names[2]);
  }
  static ThetaSchemeDescription<3> fs0(double delta_t) {
    const double theta = 1.0 - (std::sqrt(2.0) / 2.0f);
    const double theta_squigly = 1.0 - (2.0 * theta);
    const double tau = theta_squigly / (1.0 - theta);
    const double eta = 1.0 - tau;
    typedef ThetaSchemeDescription<3> ReturnType;
    ReturnType::ThetaValueArray step_one = {{tau * theta, eta * theta, eta * theta, tau * theta}};
    ReturnType::ThetaValueArray step_two = {
        {eta * theta_squigly, tau * theta_squigly, tau * theta_squigly, eta * theta_squigly}};
    ReturnType::ThetaValueArray step_three = {{tau * theta, eta * theta, eta * theta, tau * theta}};
    ReturnType::ThetaArray a;
    a[0] = step_one;
    a[1] = step_two;
    a[2] = step_three;
    ReturnType::TimestepArray c;
    c[0] = theta * delta_t;
    c[1] = theta_squigly * delta_t;
    c[2] = theta * delta_t;
    return ReturnType(a, c, scheme_names[4]);
  }
  static ThetaSchemeDescription<3> fs1(double delta_t) {
    const double theta = 1.0 - (std::sqrt(2.0) / 2.0f);
    const double theta_squigly = 1.0 - (2.0 * theta);
    const double tau = theta_squigly / (1.0 - theta);
    const double eta = 1.0 - tau;
    typedef ThetaSchemeDescription<3> ReturnType;
    ReturnType::ThetaValueArray step_one = {{tau * theta, eta * theta, theta, 0}};
    ReturnType::ThetaValueArray step_two = {{eta * theta_squigly, tau * theta_squigly, 0, theta_squigly}};
    ReturnType::ThetaValueArray step_three = {{tau * theta, eta * theta, theta, 0}};
    ReturnType::ThetaArray a;
    a[0] = step_one;
    a[1] = step_two;
    a[2] = step_three;
    ReturnType::TimestepArray c;
    c[0] = theta * delta_t;
    c[1] = theta_squigly * delta_t;
    c[2] = theta * delta_t;
    return ReturnType(a, c, scheme_names[5]);
  }
};

template <int N>
const std::vector<std::string> ThetaSchemeDescription<N>::scheme_names(scheme_names_array, scheme_names_array + 6);

template <class Stream, int N>
inline Stream& operator<<(Stream& s, ThetaSchemeDescription<N> desc) {
  s << boost::format("%s (%d-step) scheme description:\n") % desc.algo_id % N;
  for (int i = 0; i < N; ++i) {
    s << boost::format("dt_k: %e\ttheta1: %e\ttheta2: %e\ttheta3: %e\ttheta4: %e\n") // ideally this could be re-used,
                                                                                     // but that result in exception
                                                                                     // 'too-few-args"..
             %
             desc.step_sizes_[i] % desc.thetas_[i][0] % desc.thetas_[i][1] % desc.thetas_[i][2] % desc.thetas_[i][3];
  }
  return s;
}

template <class CommunicatorImp, class GridImp, template <class, class> class AnalyticalForceImp,
          template <class, class> class AnalyticalDirichletDataImp, template <class, class> class ExactPressureImp,
          template <class, class> class ExactVelocityImp, int subStepCount, int gridDim, int sigmaOrder,
          int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder>
struct ThetaSchemeTraits {
  typedef ThetaSchemeTraits<CommunicatorImp, GridImp, AnalyticalForceImp, AnalyticalDirichletDataImp, ExactPressureImp,
                            ExactVelocityImp, subStepCount, gridDim, sigmaOrder, velocityOrder, pressureOrder> ThisType;
  typedef GridImp GridType;

  static const int substep_count = subStepCount;
  typedef ThetaSchemeDescription<subStepCount> ThetaSchemeDescriptionType;
  typedef FractionalTimeProvider<ThetaSchemeDescriptionType, CommunicatorImp> TimeProviderType;

  typedef NonlinearStep::DiscreteOseenModelTraits<TimeProviderType, GridType, AnalyticalForceImp,
                                                  OseenStep::ForceAdapterFunction, AnalyticalDirichletDataImp,
                                                  typename ThetaSchemeDescriptionType::ThetaValueArray, gridDim,
                                                  sigmaOrder, velocityOrder, pressureOrder> OseenModelTraits;
  typedef typename OseenModelTraits::GridPartType GridPartType;
  typedef NonlinearStep::DiscreteOseenModelTraits<TimeProviderType, GridType, AnalyticalForceImp,
                                                  StokesStep::ForceAdapterFunction, AnalyticalDirichletDataImp,
                                                  typename ThetaSchemeDescriptionType::ThetaValueArray, gridDim,
                                                  sigmaOrder, velocityOrder, pressureOrder> StokesModelTraits;

  typedef NonlinearStep::DiscreteOseenModelTraits<TimeProviderType, GridType, AnalyticalForceImp,
                                                  NonlinearStep::ForceAdapterFunction, AnalyticalDirichletDataImp,
                                                  typename ThetaSchemeDescriptionType::ThetaValueArray, gridDim,
                                                  sigmaOrder, velocityOrder, pressureOrder> NonlinearModelTraits;

  typedef NonlinearStep::DiscreteOseenModelTraits<TimeProviderType, GridType, AnalyticalForceImp,
                                                  OseenStep::DummyForceAdapterFunction, AnalyticalDirichletDataImp,
                                                  typename ThetaSchemeDescriptionType::ThetaValueArray, gridDim,
                                                  sigmaOrder, velocityOrder, pressureOrder> OseenModelAltRhsTraits;

  typedef typename OseenModelTraits::ForceAdatperType OseenForceAdapterFunctionType;
  typedef typename OseenModelAltRhsTraits::ForceAdatperType OseenAltRhsForceAdapterFunctionType;
  typedef typename StokesModelTraits::ForceAdatperType StokesForceAdapterType;
  typedef typename NonlinearModelTraits::ForceAdatperType NonlinearForceAdapterType;

  typedef ExactPressureImp<typename OseenModelTraits::PressureFunctionSpaceType, TimeProviderType> ExactPressureType;
  typedef ExactVelocityImp<typename OseenModelTraits::VelocityFunctionSpaceType, TimeProviderType> ExactVelocityType;

  typedef ExactSolution<ThisType> ExactSolutionType;

  typedef typename OseenModelTraits::PressureFunctionSpaceType PressureFunctionSpaceType;
  typedef typename OseenModelTraits::VelocityFunctionSpaceType VelocityFunctionSpaceType;

  typedef typename OseenModelTraits::DiscreteOseenFunctionSpaceWrapperType DiscreteOseenFunctionSpaceWrapperType;

  typedef typename OseenModelTraits::DiscreteOseenFunctionWrapperType DiscreteOseenFunctionWrapperType;
  typedef typename OseenModelTraits::RealAnalyticalForceType AnalyticalForceType;
  typedef typename OseenModelTraits::AnalyticalDirichletDataType AnalyticalDirichletDataType;

  typedef Dune::StartPass<DiscreteOseenFunctionWrapperType, -1> StokesStartPassType;

  typedef CommunicatorImp CommunicatorType;

  typedef Dune::DiscreteOseenModelDefault<OseenModelTraits> OseenModelType;
  typedef Dune::DiscreteOseenModelDefault<StokesModelTraits> StokesModelType;
  typedef Dune::DiscreteOseenModelDefault<NonlinearModelTraits> NonlinearModelType;
  typedef Dune::DiscreteOseenModelDefault<OseenModelAltRhsTraits> OseenModelAltRhsType;

  typedef Dune::OseenPass<OseenModelType> OseenPassType;
  typedef Dune::OseenPass<StokesModelType> StokesPassType;
  typedef Dune::OseenPass<NonlinearModelType> NonlinearPassType;
  typedef Dune::OseenPass<OseenModelAltRhsType> OseenPassAltRhsType;
};
}
}

#endif // THETASCHEME_TRAITS_HH

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
