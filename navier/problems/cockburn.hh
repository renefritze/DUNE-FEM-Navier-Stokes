#ifndef NAVIER_PROBLEMS_COCKBURN_HH
#define NAVIER_PROBLEMS_COCKBURN_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/grid.hh>
#include <dune/stuff/math.hh>
#include <dune/stuff/parametercontainer.hh>
#include "common.hh"

namespace NavierProblems {
namespace Cockburn {

static const std::string identifier = "Cockburn";
static const bool hasExactSolution = true;

static const double P = M_PI; // pi_factor;

struct SetupCheck {
  std::string err;
  template <class Scheme, class GridPart, class... Rest>
  bool operator()(Scheme* /*scheme*/, const GridPart& gridPart, const Rest&... /*rest*/) {
    Stuff::GridDimensions<typename GridPart::GridType> grid_dim(gridPart.grid());
    bool ok = Stuff::aboutEqual(grid_dim.coord_limits[0].min(), -1.) &&
              Stuff::aboutEqual(grid_dim.coord_limits[1].min(), -1.) &&
              Stuff::aboutEqual(grid_dim.coord_limits[0].max(), 1.) &&
              Stuff::aboutEqual(grid_dim.coord_limits[1].max(), 1.);
    err = (boost::format("\n******\nSetupCheck Failed!\ngrid dimension %f,%f - %f,%f\n") %
           grid_dim.coord_limits[0].min() % grid_dim.coord_limits[1].min() % grid_dim.coord_limits[0].max() %
           grid_dim.coord_limits[1].max()).str();
    if (!ok)
      return false;
    const double v = Parameters().getParam("viscosity", -10.0);
    ok = Stuff::aboutEqual(v, 1.0);
    err = (boost::format("viscosity %f\n") % v).str();
    return ok;
  }
  std::string error() { return err; }
};

NULLFUNCTION_TP(VelocityLaplace)

template <class DomainType, class RangeType>
void VelocityEvaluate(const double /*lambda*/, const double /*time*/, const DomainType& arg, RangeType& ret) {
  const double x1 = arg[0];
  const double x2 = arg[1];
  const double exp_of_x1 = std::exp(x1);
  const double sin_of_x2 = std::sin(x2);
  const double cos_of_x2 = std::cos(x2);
  // return
  ret[0] = x2 * cos_of_x2;
  ret[0] += sin_of_x2;
  ret[0] *= -1.0 * exp_of_x1;
  ret[1] = exp_of_x1 * x2 * sin_of_x2;
}

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
    //					  ret[0] = - C_x * E * P * ( S_x * E + v * S_y * P )	+ 0.5 * P * F * S_2x;
    //					  ret[1] = - C_y * E * P * ( S_y * E - v * S_x * P )	+ 0.5 * P * F * S_2y;

    RangeType u_eval;
    VelocityEvaluate(0, time, arg, u_eval);
    ret = RangeType(0);
    //
    ret[0] = u_eval[0] * u_eval[0];
    ret[0] += u_eval[1] * std::exp(x) * (-2 * std::cos(y) + (y * std::sin(y)));
    ret[1] = u_eval[0] * u_eval[1];
    ret[1] += u_eval[1] * (-1) * u_eval[0];
  }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class PressureGradient : public Dune::TimeFunction<
    FunctionSpaceImp, PressureGradient<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef PressureGradient<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  PressureGradient(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                   const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~PressureGradient() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    const double x = arg[0];
    const double y = arg[1];
    ret[0] = 2 * std::exp(x) * std::sin(y);
    ret[1] = 2 * std::exp(x) * std::cos(y);
  }

private:
  const double parameter_a_;
  const double parameter_d_;
};

//! gd
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

  //    inline void evaluateTime( const DomainType& arg, RangeType& ret ) const {VelocityEvaluate( lambda_, 0, arg,
  // ret);}

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
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
  Velocity(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
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
    ret = 2 * std::exp(x) * std::sin(y);
  }

  template <class DiscreteFunctionSpace>
  void setShift(const DiscreteFunctionSpace& space) {
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

} // end ns
} // end ns

#endif // NAVIER_PROBLEMS_COCKBURN_HH

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
