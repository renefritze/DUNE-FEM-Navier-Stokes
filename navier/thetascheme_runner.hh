#ifndef THETASCHEME_RUNNER_HH
#define THETASCHEME_RUNNER_HH

#include <dune/navier/thetascheme_traits.hh>
#include <dune/navier/thetascheme.hh>
#include <dune/navier/thetascheme_alt_split.hh>
#include <dune/navier/global_defines.hh>
#include <dune/common/static_assert.hh>

template <class GridType, class CollectiveCommunicationType>
class ThetaschemeRunner {
private:
  typedef Dune::NavierStokes::ThetaSchemeTraits<
      CollectiveCommunicationType, GridType, NAVIER_DATA_NAMESPACE::Force, NAVIER_DATA_NAMESPACE::DirichletData,
      NAVIER_DATA_NAMESPACE::Pressure, NAVIER_DATA_NAMESPACE::Velocity, 1, // number of substeps
      GridType::dimensionworld, POLORDER, VELOCITY_POLORDER, PRESSURE_POLORDER> OneStepThetaSchemeTraitsType;
  typedef Dune::NavierStokes::ThetaScheme<OneStepThetaSchemeTraitsType> OneStepThetaSchemeType;
  typedef typename OneStepThetaSchemeTraitsType::ThetaSchemeDescriptionType OneStepThetaSchemeDescriptionType;
  typedef Dune::NavierStokes::DataOnlyScheme<OneStepThetaSchemeTraitsType> DataOnlySchemeType;
  typedef typename OneStepThetaSchemeTraitsType::ThetaSchemeDescriptionType DataOnlySchemeDescriptionType;
  typedef Dune::NavierStokes::ThetaSchemeTraits<
      CollectiveCommunicationType, GridType, NAVIER_DATA_NAMESPACE::Force, NAVIER_DATA_NAMESPACE::DirichletData,
      NAVIER_DATA_NAMESPACE::Pressure, NAVIER_DATA_NAMESPACE::Velocity, 3, // number of substeps
      GridType::dimensionworld, POLORDER, VELOCITY_POLORDER, PRESSURE_POLORDER> ThreeStepThetaSchemeTraitsType;
  typedef Dune::NavierStokes::ThetaScheme<ThreeStepThetaSchemeTraitsType> ThreeStepThetaSchemeType;
  typedef Dune::NavierStokes::ThetaSchemeAltSplitting<ThreeStepThetaSchemeTraitsType>
  ThreeStepThetaSchemeAltSplittingType;
  typedef typename ThreeStepThetaSchemeTraitsType::ThetaSchemeDescriptionType ThreeStepThetaSchemeDescriptionType;

public:
  ThetaschemeRunner(GridType& grid, CollectiveCommunicationType& comm)
    : grid_part_(grid)
    , comm_(comm) // create gridpart from passed pointer
  {
    typedef boost::is_same<typename ThreeStepThetaSchemeTraitsType::GridPartType,
                           typename OneStepThetaSchemeTraitsType::GridPartType> IdenticalGridParts;
    dune_static_assert(IdenticalGridParts::value,
                       "both traits classes need to produce the same gridpart type, otherwise you're doing it wrong");
  }

  Stuff::RunInfoTimeMap run(const int scheme_type) {
    const double dt_ = Parameters().getParam("fem.timeprovider.dt", double(0.1), Dune::ValidateGreater<double>(0.0));
    switch (scheme_type) {
      case 1:
        return OneStepThetaSchemeType(grid_part_, OneStepThetaSchemeDescriptionType::forward_euler(dt_)).run();
      case 2:
        return OneStepThetaSchemeType(grid_part_, OneStepThetaSchemeDescriptionType::backward_euler(dt_)).run();
      case 3:
        return OneStepThetaSchemeType(grid_part_, OneStepThetaSchemeDescriptionType::crank_nicholson(dt_)).run();
      case 4:
        if (Parameters().getParam("old_timestep", false))
          return ThreeStepThetaSchemeAltSplittingType(grid_part_, ThreeStepThetaSchemeDescriptionType::fs0(dt_)).run();
        else
          return ThreeStepThetaSchemeType(grid_part_, ThreeStepThetaSchemeDescriptionType::fs0(dt_)).run();
      default:
        Logger().Info() << "Using default value for theta scheme type\n";
      case 5:
        if (Parameters().getParam("old_timestep", false))
          return ThreeStepThetaSchemeAltSplittingType(grid_part_, ThreeStepThetaSchemeDescriptionType::fs1(dt_)).run();
        else
          return ThreeStepThetaSchemeType(grid_part_, ThreeStepThetaSchemeDescriptionType::fs1(dt_)).run();
      case 0:
        return DataOnlySchemeType(grid_part_, DataOnlySchemeDescriptionType::crank_nicholson(dt_)).run();
    }
  }

private:
  typename ThreeStepThetaSchemeTraitsType::GridPartType grid_part_;
  CollectiveCommunicationType& comm_;
};

#endif // THETASCHEME_RUNNER_HH

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
