#ifndef WEIGHED_FUNCTION_HH
#define WEIGHED_FUNCTION_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>

namespace Dune {
namespace NavierStokes {

/** simple wrapper that allows us to create functions
 *
 */
template <class FunctionSpaceImp, class TimeProviderImp, class FunctionImp>
class WeighedIntersectionFunction : public Dune::IntersectionTimeFunction<
    FunctionSpaceImp, WeighedIntersectionFunction<FunctionSpaceImp, TimeProviderImp, FunctionImp>, TimeProviderImp> {
public:
  typedef WeighedIntersectionFunction<FunctionSpaceImp, TimeProviderImp, FunctionImp> ThisType;
  typedef Dune::IntersectionTimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef FunctionImp FunctionType;
  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  WeighedIntersectionFunction(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                              const double weight_a = 1.0, const double weight_b = 0.0)
    : BaseType(timeprovider, space)
    , weight_a_(weight_a)
    , weight_b_(weight_b)
    , function_(timeprovider, space) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~WeighedIntersectionFunction() {}
  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    RangeType b;
    function_.evaluateTime(time, arg, ret);
    function_.evaluateTime(time - BaseType::timeProvider_.deltaT(), arg, b);
    ret *= weight_a_;
    b *= weight_b_;
    ret += b;
  }

  template <class IntersectionType>
  void evaluateTime(const double time, const DomainType& arg, RangeType& ret,
                    const IntersectionType& intersection) const {
    RangeType b;
    function_.evaluateTime(time, arg, ret, intersection);
    function_.evaluateTime(time - BaseType::timeProvider_.deltaT(), arg, b, intersection);
    ret *= weight_a_;
    b *= weight_b_;
    ret += b;
  }

private:
  const double weight_a_;
  const double weight_b_;
  FunctionType function_;
};

} // namespace Dune
} // namespace NavierStokes
#endif // WEIGHED_FUNCTION_HH

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
