#ifndef REAL_BVP_HH
#define REAL_BVP_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/oseen/boundarydata.hh>
#include "common.hh"

namespace NavierProblems {

namespace Null {
static const std::string identifier = "NULL";
static const bool hasExactSolution = true;
ALLGOOD_SETUPCHECK;

NULLFUNCTION_TP(PressureGradient)
NULLFUNCTION_TP(Pressure)
NULLFUNCTION_TP(VelocityConvection)
NULLFUNCTION_TP(VelocityLaplace)
NULLFUNCTION_TP(Velocity)
NULLFUNCTION_TP(Force)
NULLFUNCTION_TP_BOUNDARY(DirichletData)
}

namespace BVP_A {

static const std::string identifier = "BVP";
static const bool hasExactSolution = false;
ALLGOOD_SETUPCHECK;

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
    : BaseType(timeprovider, space)
    , z_max(Parameters().getParam("z_max", 3.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~DirichletData() {}
  void evaluateTime(const double /*time*/, const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }

  template <class IntersectionType>
  void evaluateTime(const double time, const DomainType& arg, RangeType& ret,
                    const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 3, "__CLASS__ evaluate not implemented for world dimension");
    DomainType normal(0);
    normal[2] = 1;
    const double gd_factor = time * Parameters().getParam("gd_factor", 1.0);
    if (arg[2] > 0.0 && arg[2] < z_max) {
      normal *= 0.0;
    }
    ret = normal;
    ret *= gd_factor;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double z_max;
};

NULLFUNCTION_TP(PressureGradient)
NULLFUNCTION_TP(Pressure)
NULLFUNCTION_TP(VelocityConvection)
NULLFUNCTION_TP(VelocityLaplace)
NULLFUNCTION_TP(Velocity)
NULLFUNCTION_TP(Force)
}

namespace BVP {

static const std::string identifier = "BVP";
static const bool hasExactSolution = false;
// currently not possible, see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=45114
//    template < class FunctionSpaceImp, class TimeProviderImp >
//    using DirichletData = Stuff::InstationaryBoundaryFluxFunction<FunctionSpaceImp, TimeProviderImp>;
//    template < class FunctionSpaceImp, class TimeProviderImp >
//    class DirichletData : public Stuff::InstationaryBoundaryFluxFunction<FunctionSpaceImp, TimeProviderImp>
//    {
//	    typedef Stuff::InstationaryBoundaryFluxFunction<FunctionSpaceImp, TimeProviderImp>
//		BaseType;
//	public:
//	    DirichletData( const TimeProviderImp& timeprovider,
//			   const FunctionSpaceImp& space,
//			   const double /*constant*/ = 0.0 )
//	    : BaseType ( timeprovider, space )
//	    {}
//    };

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
    : BaseType(timeprovider, space)
    , z_max(Parameters().getParam("z_max", 3.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~DirichletData() {}
  void evaluateTime(const double /*time*/, const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }

  template <class IntersectionType>
  void evaluateTime(const double time, const DomainType& arg, RangeType& ret,
                    const IntersectionType& intersection) const {
    const int id = intersection.boundaryId();
    typedef Dune::FieldVector<typename IntersectionType::ctype, IntersectionType::dimension - 1> LocalVectorType;

    LocalVectorType center = Stuff::getBarycenterLocal(intersection.intersectionSelfLocal());
    RangeType normal = intersection.unitOuterNormal(center);
    ret = normal;
    double factor = Parameters().getParam("gd_factor", 1.0) * time;
    switch (id) {
      case 1: {
        factor = 0;
        break;
      }
      case 2: {
        factor *= -1;
        break;
      }
      case 6:
      case 5:
      case 4:
      case 3: {
        factor *= 1;
        break;
      }
      default:
        assert(false);
    }
    ret *= factor * time;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double z_max;
};

NULLFUNCTION_TP(PressureGradient)
NULLFUNCTION_TP(Pressure)
NULLFUNCTION_TP(VelocityConvection)
NULLFUNCTION_TP(VelocityLaplace)

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

  void evaluateTime(const double /*time*/, const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }

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

NULLFUNCTION_TP(Force)
} // end ns
} // end ns

#endif // REAL_BVP_HH

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
