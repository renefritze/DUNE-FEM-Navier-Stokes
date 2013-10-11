#ifndef NAVIER_PROBLEMS_TRIVIAL_HH
#define NAVIER_PROBLEMS_TRIVIAL_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include "common.hh"

namespace NavierProblems {
namespace Trivial {

static const std::string identifier = "Trvial";
static const bool hasExactSolution = true;
ALLGOOD_SETUPCHECK;

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

  ~Force() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    const double x = arg[0];
    const double y = arg[1];
    ret = RangeType(0);
    // zeitableitung u
    ret[0] += 2;
    ret[1] += 0;
    // pressure gradient
    ret[0] += -1 * time;
    ret[1] += 0;
    // conv
    if (!Parameters().getParam("navier_no_convection", false)) {
      assert(false);
      ret[0] += -x;
      ret[1] += -y;
    }
  }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double /*lambda*/, const double time, const DomainType& /*arg*/, RangeType& ret) {
  //	const double x				= arg[0];
  //	const double y				= arg[1];
  ret[0] = 2 * time;
  ret[1] = 0;
}

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
   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
   **/
  VelocityConvection(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double viscosity = 0.0,
                     const double alpha = 0.0)
    : BaseType(timeprovider, space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  ~VelocityConvection() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    VelocityEvaluate(0.0f, time, arg, ret);
  }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

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
  *
  *  doing nothing besides Base init
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
    dune_static_assert(dim_ == 2, "__CLASS__ evaluate not implemented for world dimension");
    VelocityEvaluate(0.0, time, arg, ret);
  }

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    VelocityEvaluate(0.0, time, arg, ret);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
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
           const double /*parameter_a*/ = M_PI / 2.0, const double /*parameter_d */ = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~Velocity() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
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
    dune_static_assert(dim_ == 2, "__CLASS__ evaluate not implemented for world dimension");
    const double x = arg[0];
    //			const double y				= arg[1];
    //			const double v				= Parameters().getParam( "viscosity", 1.0 );
    //			const double F				= std::exp( -16 * std::pow( M_PI, 2 ) * time );
    //			const double C1				= std::cos(4*M_PI* ( x + 0.25 ) );
    //			const double C2				= std::cos(4*M_PI* ( y + 0.5 ) );

    ret = (0.5 - x) * time;
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

template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityLaplace
    : public Dune::TimeFunction<FunctionSpaceImp, VelocityLaplace<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityLaplace<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  VelocityLaplace(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                  const double /*parameter_a*/ = M_PI / 2.0, const double /*parameter_d*/ = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~VelocityLaplace() {}

  void evaluateTime(const double /*time*/, const DomainType& /*arg*/, RangeType& ret) const {
    dune_static_assert(dim_ == 2, "__CLASS__ evaluate not implemented for world dimension");
    ret = RangeType(0);
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
class PressureGradient : public Dune::TimeFunction<
    FunctionSpaceImp, PressureGradient<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef PressureGradient<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  PressureGradient(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                   const double /*parameter_a*/ = M_PI / 2.0, const double /*parameter_d*/ = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~PressureGradient() {}

  void evaluateTime(const double time, const DomainType& /*arg*/, RangeType& ret) const {
    dune_static_assert(dim_ == 2, "__CLASS__ evaluate not implemented for world dimension");
    //			const double x				= arg[0];
    //			const double y				= arg[1];
    ret[0] = -1 * time;
    ret[1] = 0;
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

} // end ns
} // end ns

#endif // NAVIER_PROBLEMS_TRIVIAL_HH

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
