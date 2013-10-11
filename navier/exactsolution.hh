#ifndef EXACTSOLUTION_HH
#define EXACTSOLUTION_HH

#include <dune/stuff/customprojection.hh>

namespace Dune {
namespace NavierStokes {

template <class TraitsImp>
class ExactSolution : public TraitsImp::DiscreteOseenFunctionWrapperType {

  typedef TraitsImp TraitsType;
  typedef typename TraitsType::DiscreteOseenFunctionWrapperType BaseType;

  const typename TraitsType::TimeProviderType& timeprovider_;
  typename TraitsType::PressureFunctionSpaceType continousPressureSpace_;
  typename TraitsType::VelocityFunctionSpaceType continousVelocitySpace_;
  typename TraitsType::ExactVelocityType velocity_;
  typename TraitsType::ExactPressureType pressure_;

public:
  ExactSolution(const typename TraitsType::TimeProviderType& timeprovider, typename TraitsType::GridPartType& gridPart,
                typename TraitsType::DiscreteOseenFunctionSpaceWrapperType& space_wrapper)
    : BaseType("exact", space_wrapper, gridPart)
    , timeprovider_(timeprovider)
    , velocity_(timeprovider_, continousVelocitySpace_)
    , pressure_(timeprovider_, continousPressureSpace_) {
    project();
  }

  void project() { BaseType::projectInto(velocity_, pressure_); }

  const typename TraitsType::ExactVelocityType& exactVelocity() const { return velocity_; }

  const typename TraitsType::ExactPressureType& exactPressure() const { return pressure_; }

  typename TraitsType::ExactVelocityType& exactVelocity() { return velocity_; }

  typename TraitsType::ExactPressureType& exactPressure() { return pressure_; }

  void atTime(const double time, BaseType& dest) const {
    Dune::BetterL2Projection::project(time, pressure_, dest.discretePressure());

    Dune::BetterL2Projection::project(time, velocity_, dest.discreteVelocity());
  }

public:
  typedef typename BaseType::Traits::FunctionTupleType FunctionTupleType;
};

} // end namespace NavierStokes
} // end namespace Dune

#endif // EXACTSOLUTION_HH

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
