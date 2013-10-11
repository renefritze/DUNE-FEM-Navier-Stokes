#ifndef NAVIER_PROBLEMS_TAYLOR_HH
#define NAVIER_PROBLEMS_TAYLOR_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include "common.hh"

namespace NavierProblems {
namespace Taylor {

static const std::string identifier = "Taylor";
static const bool hasExactSolution = true;
ALLGOOD_SETUPCHECK;
static const double P = M_PI; // pi_factor;

template <class FunctionSpaceImp, class TimeProviderImp>
class Force : public Dune::TimeFunction<FunctionSpaceImp, Force<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Force<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
    *  \brief  constructor
    *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
    **/
  Force(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double viscosity = 0.0,
        const double alpha = 0.0)
    : BaseType(timeprovider, space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  /**
    *  \brief  destructor
    *  doing nothing
    **/
  ~Force() {}

  /**
    *  \brief  evaluates the force
    *  \param  arg
    *          point to evaluate at
    *  \param  ret
    *          value of force at given point
    **/
  inline void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    const double x = arg[0];
    const double y = arg[1];
    const double v = viscosity_;
    const double E = std::exp(-2 * std::pow(P, 2) * viscosity_ * time);
    const double F = std::exp(-4 * std::pow(P, 2) * viscosity_ * time);
    const double S_x = std::sin(P * x);
    const double S_y = std::sin(P * y);
    const double S_2x = std::sin(2 * P * x);
    const double S_2y = std::sin(2 * P * y);
    const double C_x = std::cos(P * x);
    const double C_y = std::cos(P * y);
    //					  ret[0] = - C_x * E * P * ( S_x * E + v * S_y * P )	+ 0.5 * P * F * S_2x;
    //					  ret[1] = - C_y * E * P * ( S_y * E - v * S_x * P )	+ 0.5 * P * F * S_2y;

    RangeType u_eval;
    VelocityEvaluate(0, time, arg, u_eval);
    ret = RangeType(0);
    //					  ret += alpha * u_eval;

    // laplace
    ret[0] -= +2 * P * P * C_x * S_y * E;
    ret[1] -= +2 * P * P * S_x * C_y * E;
    // conv
    //					  ret[0] += - P * E * E * S_x * C_x ;
    //					  ret[1] += - P * E * E * S_y * C_y ;
    // p-grad
    ret[0] += 0.5 * P * S_2x * F;
    ret[1] += 0.5 * P * S_2y * F;
  }

private:
  const double viscosity_;
  const double alpha_;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double /*lambda*/, const double time, const DomainType& arg, RangeType& ret) {
  const double x = arg[0];
  const double y = arg[1];
  const double v = Parameters().getParam("viscosity", 1.0);
  const double E = std::exp(-2 * std::pow(P, 2) * v * time);
  const double S_x = std::sin(P * x);
  const double S_y = std::sin(P * y);
  const double C_x = std::cos(P * x);
  const double C_y = std::cos(P * y);

  ret[0] = (-1 / v) * C_x * S_y * E;
  ret[1] = (1 / v) * S_x * C_y * E;
}

/**
 *  \brief  describes the dirichlet boundary data
 *
 *  \tparam DirichletTraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template <class FunctionSpaceImp, class TimeProviderImp>
class DirichletData : public Dune::IntersectionTimeFunction<
    FunctionSpaceImp, DirichletData<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef DirichletData<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::IntersectionTimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
    *  \brief  constructor
    *  \param  viscosity,alpha   dummies
    **/
  DirichletData(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double /*viscosity*/ = 0.0,
                const double /*alpha*/ = 0.0)
    : BaseType(timeprovider, space) {}

  /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
  ~DirichletData() {}

  template <class IntersectionType>
  void evaluateTime(const double time, const DomainType& arg, RangeType& ret,
                    const IntersectionType& /*intersection */) const {
    dune_static_assert(FunctionSpaceImp::dimDomain == 2, "__CLASS__ evaluate not implemented for world dimension");
    VelocityEvaluate(0.0, time, arg, ret);
  }
  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    VelocityEvaluate(0.0, time, arg, ret);
  }
};
template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityConvection : public Dune::TimeFunction<
    FunctionSpaceImp, VelocityConvection<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityConvection<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *
   *  doing nothing besides Base init
   **/
  VelocityConvection(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                     const double /*parameter_a*/ = M_PI / 2.0, const double /*parameter_d*/ = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
  ~VelocityConvection() {}

  template <class IntersectionType>
  void evaluateTime(const double time, const DomainType& arg, RangeType& ret,
                    const IntersectionType& /*intersection */) const {
    dune_static_assert(FunctionSpaceImp::dimDomain == 2, "__CLASS__ evaluate not implemented for world dimension");
    VelocityEvaluate(lambda_, time, arg, ret);
  }
  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    VelocityEvaluate(lambda_, time, arg, ret);
  }

private:
  const double lambda_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class Velocity
    : public Dune::TimeFunction<FunctionSpaceImp, Velocity<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Velocity<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *
   *  doing nothing besides Base init
   **/
  Velocity(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
           const double /*parameter_a*/ = M_PI / 2.0, const double /*parameter_d*/ = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
  ~Velocity() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(FunctionSpaceImp::dimDomain == 2, "__CLASS__ evaluate not implemented for world dimension");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

  /**
   * \brief  evaluates the dirichlet data
   * \param  arg
   *         point to evaluate at
   * \param  ret
   *         value of dirichlet boundary data at given point
   **/
  //					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double lambda_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class Pressure
    : public Dune::TimeFunction<FunctionSpaceImp, Pressure<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Pressure<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *
   *  doing nothing besides Base init
   **/
  Pressure(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
           const double /*parameter_a*/ = M_PI / 2.0, const double /*parameter_d*/ = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0))
    , shift_(0.0) {}

  /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
  ~Pressure() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(FunctionSpaceImp::dimDomain == 2, "__CLASS__ evaluate not implemented for world dimension");
    const double x = arg[0];
    const double y = arg[1];
    const double v = Parameters().getParam("viscosity", 1.0);
    const double F = std::exp(-4 * std::pow(P, 2) * v * time);
    const double C_2x = std::cos(2 * P * x);
    const double C_2y = std::cos(2 * P * y);
    ret = (-1 / (4 * v)) * (C_2x + C_2y) * F;
  }

  template <class DiscreteFunctionSpace>
  void setShift(const DiscreteFunctionSpace& /*space*/) {
    //					shift_ = -1 * Stuff::meanValue( *this, space );
  }

  /**
   * \brief  evaluates the dirichlet data
   * \param  arg
   *         point to evaluate at
   * \param  ret
   *         value of dirichlet boundary data at given point
   **/
  //					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double lambda_;
  double shift_;
};

NULLFUNCTION_TP(VelocityLaplace)
NULLFUNCTION_TP(PressureGradient)
} // end ns
} // end ns

#endif // NAVIER_PROBLEMS_TAYLOR_HH

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
