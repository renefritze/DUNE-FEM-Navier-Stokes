#ifndef RHSADAPTER_HH
#define RHSADAPTER_HH

#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/stuff/printing.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/navier/problems.hh>

namespace Dune {
namespace NavierStokes {
namespace OseenStep {
/** \brief take previous step solution \f$u_{k-1}\f$ and analytical RHS to form function to be passed to either
  StokesStep
  * given analytical force \f$f_{ana}\f$ and discrete function \f$u\f$ representing previous time step's velocity
  solution,
  *	this calculates new right hand side \f$F := \theta_{3} \delta{t_n} f_{k} + \theta_{2} \delta{t_n}f_{k-1}
              + u_{k-1} + \frac{\theta_{1} \delta{t_n}}{Re} \Delta u_{k-1}
               + \theta_{1} \delta{t_n} - \left( _{k-1} \cdot \nabla \right) _{k-1} \f$
  * \note theta value array is 0 based, so all indices have a -1 offset to the paper
  */
template <class TimeProviderType, class AnalyticalForceType, class DiscreteVelocityFunctionType, class ThetaValuesType>
class ForceAdapterFunction : public DiscreteVelocityFunctionType {
protected:
  typedef ForceAdapterFunction<TimeProviderType, AnalyticalForceType, DiscreteVelocityFunctionType, ThetaValuesType>
  ThisType;
  const TimeProviderType& timeProvider_;
  const AnalyticalForceType& force_;
  const double reynolds_;
  const ThetaValuesType& theta_values_;

public:
  typedef DiscreteVelocityFunctionType BaseType;

  //! this sginature is used in the first stokes where we have analytical data to derive from
  ForceAdapterFunction(const TimeProviderType& timeProvider, const DiscreteVelocityFunctionType& velocity,
                       const AnalyticalForceType& force, const double reynolds, const ThetaValuesType& theta_values,
                       int /*polOrd*/ = -1)
    : BaseType("stokes-ana-rhsdapater", velocity.space())
    , timeProvider_(timeProvider)
    , force_(force)
    , reynolds_(reynolds)
    , theta_values_(theta_values) {

    //						NAVIER_DATA_NAMESPACE
    typedef typename DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType VelocityFunctionSpaceType;
    VelocityFunctionSpaceType continousVelocitySpace_;

    typedef NAVIER_DATA_NAMESPACE::VelocityLaplace<VelocityFunctionSpaceType, TimeProviderType> VelocityLaplace;
    VelocityLaplace velocity_laplace(timeProvider_, continousVelocitySpace_);
    typedef NAVIER_DATA_NAMESPACE::VelocityConvection<VelocityFunctionSpaceType, TimeProviderType> VelocityConvection;
    VelocityConvection velocity_convection(timeProvider_, continousVelocitySpace_);
    typedef NAVIER_DATA_NAMESPACE::PressureGradient<VelocityFunctionSpaceType, TimeProviderType> PressureGradient;
    PressureGradient pressure_gradient(timeProvider_, continousVelocitySpace_);

    DiscreteVelocityFunctionType velocity_convection_discrete("velocity_convection_discrete", velocity.space());
    DiscreteVelocityFunctionType velocity_laplace_discrete("velocity_laplace_discrete", velocity.space());
    DiscreteVelocityFunctionType pressure_gradient_discrete("velocity_laplace_discrete", velocity.space());

    Dune::BetterL2Projection // we need evals from the _previous_ (t_{k-1}) step
        ::project(timeProvider_.previousSubTime(), velocity_convection, velocity_convection_discrete);
    Dune::BetterL2Projection::project(timeProvider_.previousSubTime(), velocity_laplace, velocity_laplace_discrete);
    Dune::BetterL2Projection::project(timeProvider_.previousSubTime(), pressure_gradient, pressure_gradient_discrete);

    AddCommon(velocity, velocity_convection_discrete, velocity_laplace_discrete, pressure_gradient_discrete);
  }

  //! this signature is used in all other stokes steps where we get the data from the previous step's discretisation
  template <class RhsContainerType>
  ForceAdapterFunction(const TimeProviderType& timeProvider, const DiscreteVelocityFunctionType& velocity,
                       const AnalyticalForceType& force, const double reynolds, const ThetaValuesType& theta_values,
                       const RhsContainerType& rhs_container, int /*polOrd*/ = -1)
    : BaseType("stokes-ana-rhsdapater", velocity.space())
    , timeProvider_(timeProvider)
    , force_(force)
    , reynolds_(reynolds)
    , theta_values_(theta_values) {
    AddCommon(velocity, rhs_container.convection, rhs_container.velocity_laplace, rhs_container.pressure_gradient);
  }

protected:
  /** \f$F := \theta_{3} \delta{t_n} f_{k} + \theta_{2} \delta{t_n}f_{k-1}
  + u_{k-1} + \frac{\theta_{1} \delta{t_n}}{Re} \Delta u_{k-1}
   + \theta_{1} \delta{t_n} - \left( _{k-1} \cdot \nabla \right) _{k-1} \f$
  \note theta value array is 0 based, so all indices have a -1 offset to the paper**/
  void AddCommon(const DiscreteVelocityFunctionType& velocity, const DiscreteVelocityFunctionType& convection,
                 const DiscreteVelocityFunctionType& velocity_laplace,
                 const DiscreteVelocityFunctionType& pressure_gradient) {
    const double dt_n = timeProvider_.deltaT();
    this->clear();

    DiscreteVelocityFunctionType tmp("rhs-ana-tmp", velocity.space());
    Dune::BetterL2Projection::project(timeProvider_, force_, tmp);
    tmp *= (theta_values_[3]);
    *this += tmp;

    Dune::BetterL2Projection::project(timeProvider_.previousSubTime(), force_, tmp);
    tmp *= (theta_values_[2]);
    *this += tmp;

    tmp.assign(velocity_laplace);
    tmp *= (theta_values_[1]) / reynolds_;
    *this += tmp; // this = f + beta_re_qoutient * laplace

    tmp.assign(convection);
    tmp *= (theta_values_[1]);
    *this -= tmp;

    tmp.assign(pressure_gradient);
    tmp *= (theta_values_[1]);
    *this -= tmp;
    //*this *= dt_n ;
    tmp.assign(velocity);
    tmp /= dt_n;

    *this += tmp;
  }
};

//! this is just a type match wrapper for a velocity discrete function
template <class TimeProviderType, class AnalyticalForceType, class DiscreteVelocityFunctionType, class ThetaValuesType>
class DummyForceAdapterFunction : public DiscreteVelocityFunctionType {
public:
  DummyForceAdapterFunction(const DiscreteVelocityFunctionType& other) : DiscreteVelocityFunctionType(other) {}
};

} // end namespace OseenStep
namespace StokesStep {
/** \brief take previous step solution and analytical RHS to form function to be passed to either StokesStep
  * given analytical force \f$f_{ana}\f$ and discrete function \f$u\f$ representing previous time step's velocity
  *solution,
  *	this calculates new right hand side \f$f := f_{ana} + frac{1}{\theta \tau}u + \frac{\beta}{Re} \Delta u - \left( u
  *\cdot \nabla \right) u \f$
  *
  */
template <class TimeProviderType, class AnalyticalForceType, class DiscreteVelocityFunctionType, class ThetaValuesType>
class ForceAdapterFunction : public DiscreteVelocityFunctionType {
protected:
  typedef ForceAdapterFunction<TimeProviderType, AnalyticalForceType, DiscreteVelocityFunctionType, ThetaValuesType>
  ThisType;
  const TimeProviderType& timeProvider_;
  const AnalyticalForceType& force_;

public:
  typedef DiscreteVelocityFunctionType BaseType;

  //! this sginature is used in the first stokes where we have analytical data to derive from
  template <class DiscretizationWeightsType>
  ForceAdapterFunction(const TimeProviderType& timeProvider, const DiscreteVelocityFunctionType& velocity,
                       const AnalyticalForceType& force, const DiscretizationWeightsType& weights, int /*polOrd*/ = -1)
    : BaseType("stokes-ana-rhsdapater", velocity.space())
    , timeProvider_(timeProvider)
    , force_(force) {

    //						NAVIER_DATA_NAMESPACE
    typedef typename DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType VelocityFunctionSpaceType;
    VelocityFunctionSpaceType continousVelocitySpace_;

    typedef NAVIER_DATA_NAMESPACE::VelocityLaplace<VelocityFunctionSpaceType, TimeProviderType> VelocityLaplace;
    VelocityLaplace velocity_laplace(timeProvider_, continousVelocitySpace_);
    typedef NAVIER_DATA_NAMESPACE::VelocityConvection<VelocityFunctionSpaceType, TimeProviderType> VelocityConvection;
    VelocityConvection velocity_convection(timeProvider_, continousVelocitySpace_);

    DiscreteVelocityFunctionType velocity_convection_discrete("velocity_convection_discrete", velocity.space());
    DiscreteVelocityFunctionType velocity_laplace_discrete("velocity_laplace_discrete", velocity.space());

    Dune::BetterL2Projection // we need evals from the _previous_ (t_0) step
        ::project(timeProvider_.previousSubTime(), velocity_convection, velocity_convection_discrete);
    Dune::BetterL2Projection::project(timeProvider_.previousSubTime(), velocity_laplace, velocity_laplace_discrete);

    AddCommon(velocity, velocity_convection_discrete, velocity_laplace_discrete, weights);
  }

  //! this signature is used in all other stokes steps where we get the data from the previous step's discretisation
  template <class RhsContainerType, class DiscretizationWeightsType>
  ForceAdapterFunction(const TimeProviderType& timeProvider, const DiscreteVelocityFunctionType& velocity,
                       const AnalyticalForceType& force, const DiscretizationWeightsType& weights,
                       const RhsContainerType& rhs_container, int /*polOrd*/ = -1)
    : BaseType("stokes-ana-rhsdapater", velocity.space())
    , timeProvider_(timeProvider)
    , force_(force) {
    AddCommon(velocity, rhs_container.convection, rhs_container.velocity_laplace, weights);
  }

protected:
  //! F = alpha*f_{n+theta} beta*+f_{n} \beta / \Re * laplace u + ( 1/(theta * tau) ) u - ( u * nable ) u
  template <class DiscretizationWeightsType>
  void AddCommon(const DiscreteVelocityFunctionType& velocity, const DiscreteVelocityFunctionType& convection,
                 const DiscreteVelocityFunctionType& velocity_laplace, const DiscretizationWeightsType& weights) {
    Dune::BetterL2Projection::project(timeProvider_, force_, *this); // this = f_{n+theta}
    *this *= weights.alpha;

    DiscreteVelocityFunctionType tmp("rhs-ana-tmp", velocity.space());
    Dune::BetterL2Projection::project(timeProvider_.previousSubTime(), force_, tmp); // tmp = f_{n+theta}
    tmp *= weights.beta;
    *this += tmp;

    tmp.assign(velocity_laplace);
    tmp *= weights.beta * weights.viscosity;
    *this += tmp; // this = f + beta_re_qoutient * laplace

    *this -= convection;

    *this *= weights.theta_times_delta_t;
    *this += velocity;
  }
};

template <class TimeProviderType, class AnalyticalDirichletType>
class DirichletAdapterFunction
    : public Fem::Function<typename AnalyticalDirichletType::FunctionSpaceType,
                           DirichletAdapterFunction<TimeProviderType, AnalyticalDirichletType>> {
protected:
  typedef DirichletAdapterFunction<TimeProviderType, AnalyticalDirichletType> ThisType;
  typedef Fem::Function<typename AnalyticalDirichletType::FunctionSpaceType, ThisType> BaseType;
  const TimeProviderType& timeProvider_;
  const AnalyticalDirichletType gd_;

public:
  DirichletAdapterFunction(const TimeProviderType& timeProvider,
                           const typename AnalyticalDirichletType::FunctionSpaceType space)
    : BaseType(space)
    , timeProvider_(timeProvider)
    , gd_(space) {}

  template <class IntersectionIteratorType>
  inline void evaluate(const typename AnalyticalDirichletType::DomainType& arg,
                       typename AnalyticalDirichletType::RangeType& ret, const IntersectionIteratorType intIt) const {
    const double time = timeProvider_.subTime();
    gd_.evaluate(time, arg, ret, intIt);
  }

  inline void evaluate(const typename AnalyticalDirichletType::DomainType& /*arg*/,
                       typename AnalyticalDirichletType::RangeType& /*ret*/) const {
    NEEDS_IMPLEMENTATION
  }
};
template <template <class> class DiricheltDataImp, class TimeProviderType>
struct DirichletAdapterFunctionTraits {

  template <class FunctionSpaceImp, class GridPartImp>
  struct Implementation {
    typedef DirichletAdapterFunction<TimeProviderType, DiricheltDataImp<FunctionSpaceImp>> AnalyticalDirichletDataType;

    template <class DiscreteOseenFunctionWrapper>
    static AnalyticalDirichletDataType getInstance(const TimeProviderType& timeProvider,
                                                   const DiscreteOseenFunctionWrapper& wrapper) {
      return AnalyticalDirichletDataType(timeProvider, wrapper.discreteVelocitySpace());
    }
  };
};

} // namespace StokesStep
namespace NonlinearStep {
/** \brief take previous step solution and analytical RHS to form function to be passed to localdg code
  * given analytical force \f$f_{ana}\f$ and discrete functions \f$u,p\f$ representing previous time step's velocity and
  *pressure solution,
  *	this calculates new right hand side \f$f := f_{ana} + \frac{\alpha}{Re} \Delta u - \nabla p +
  *\frac{1}{(1-2\theta)dt} u\f$
  *
  */
template <class TimeProviderType, class AnalyticalForceType, class DiscreteVelocityFunctionType, class ThetaValuesType>
class ForceAdapterFunction : public DiscreteVelocityFunctionType {
protected:
  typedef ForceAdapterFunction<TimeProviderType, AnalyticalForceType, DiscreteVelocityFunctionType, ThetaValuesType>
  ThisType;
  typedef DiscreteVelocityFunctionType BaseType;
  const TimeProviderType& timeProvider_;

public:
  template <class RhsContainerType, class DiscretizationWeightsType>
  ForceAdapterFunction(const TimeProviderType& timeProvider, const DiscreteVelocityFunctionType& velocity,
                       const AnalyticalForceType& force, const DiscretizationWeightsType& weights,
                       const RhsContainerType& rhs_container, int /*polOrd */ = -1)
    : BaseType("nonlinear-rhsdapater", velocity.space())
    , timeProvider_(timeProvider) {
    // F = f + \alpha \Re \delta u - \nabla p + ( 1/(1-2 \theta) ) * u
    Dune::BetterL2Projection::project(timeProvider_, force, *this); // this = f_{n+theta}
    *this *= weights.beta;

    DiscreteVelocityFunctionType tmp("rhs-ana-tmp", velocity.space());
    Dune::BetterL2Projection::project(timeProvider_.previousSubTime(), force, tmp); // tmp = f_{n+theta}
    tmp *= weights.alpha;
    *this += tmp;

    tmp.assign(rhs_container.velocity_laplace);
    tmp *= weights.alpha * weights.viscosity;
    *this += tmp; // this = f + beta_re_qoutient * laplace

    tmp.assign(rhs_container.pressure_gradient);
    tmp *= weights.theta_times_delta_t;
    *this -= tmp;

    *this *= weights.one_neg_two_theta_dt;
    *this += velocity;
  }
};
} // end namespace NonlinearStep

} // end namespace NavierStokes
} // end namespace Dune

#endif // RHSADAPTER_HH

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
