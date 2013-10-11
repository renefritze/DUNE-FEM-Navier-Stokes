#ifndef NAVIER_PROBLEMS_RUNTIME_HH
#define NAVIER_PROBLEMS_RUNTIME_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/runtimefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include "common.hh"

namespace NavierProblems {
namespace Runtime {

static const std::string identifier = "Runtime";
static const bool hasExactSolution = true;
ALLGOOD_SETUPCHECK;

NV_RUNTIME_FUNC(Force);
NV_RUNTIME_FUNC(VelocityConvection);
NV_RUNTIME_FUNC(Velocity);
NV_RUNTIME_FUNC(Pressure);
NV_RUNTIME_FUNC(VelocityLaplace);
NV_RUNTIME_FUNC(PressureGradient);
NV_RUNTIME_FUNC(Beta);

template <class FunctionSpaceImp, class TimeProviderImp>
class DirichletData : public Dune::IntersectionTimeFunction<
    FunctionSpaceImp, DirichletData<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef DirichletData<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::IntersectionTimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef Velocity<FunctionSpaceImp, TimeProviderImp> VelocityType;
  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  DirichletData(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double /*viscosity*/ = 0.0,
                const double /*alpha*/ = 0.0)
    : BaseType(timeprovider, space)
    , velocity_(timeprovider, space) {}

  ~DirichletData() {}

  template <class IntersectionType>
  void evaluateTime(const double time, const DomainType& arg, RangeType& ret,
                    const IntersectionType& /*intersection */) const {
    velocity_.evaluateTime(time, arg, ret);
  }

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    velocity_.evaluateTime(time, arg, ret);
  }

private:
  const VelocityType velocity_;
};

} // end ns

} // end ns

#endif // NAVIER_PROBLEMS_RUNTIME_HH

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
