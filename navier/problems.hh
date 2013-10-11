#ifndef NAVIER_PROBLEMS_HH
#define NAVIER_PROBLEMS_HH

#ifndef NAVIER_DATA_NAMESPACE
#define NAVIER_DATA_NAMESPACE NavierProblems::Trivial
#endif
#ifndef OSEEN_DATA_NAMESPACE
#define OSEEN_DATA_NAMESPACE NavierProblems::Trivial
#endif

#include "problems/timedisc.hh"
#include "problems/trivial.hh"
#include "problems/taylor.hh"
#include "problems/heat.hh"
#include "problems/testcase2d.hh"
#include "problems/testcase3d.hh"
#include "problems/cockburn.hh"
#include "problems/real_bvp.hh"
#include "problems/2dtube.hh"
#include "problems/runtime.hh"
#include "problems/damped.hh"

namespace NavierProblems {
namespace {
//! a dummy interface class that throws errors if we've forgotten to implement data funcs in given namespace
template <class FunctionSpaceImp, class TimeProviderImp>
struct Interface {
  typedef OSEEN_DATA_NAMESPACE::PressureGradient<FunctionSpaceImp, TimeProviderImp> PressureGradient;
  typedef OSEEN_DATA_NAMESPACE::VelocityConvection<FunctionSpaceImp, TimeProviderImp> VelocityConvection;
  typedef OSEEN_DATA_NAMESPACE::VelocityLaplace<FunctionSpaceImp, TimeProviderImp> VelocityLaplace;
  typedef OSEEN_DATA_NAMESPACE::Velocity<FunctionSpaceImp, TimeProviderImp> Velocity;
  typedef OSEEN_DATA_NAMESPACE::Force<FunctionSpaceImp, TimeProviderImp> Force;
  typedef OSEEN_DATA_NAMESPACE::DirichletData<FunctionSpaceImp, TimeProviderImp> DirichletData;
};
}
}

#endif // NAVIER_PROBLEMS_HH

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
