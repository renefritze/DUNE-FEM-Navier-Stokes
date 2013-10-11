#ifndef OSEEN_HH
#define OSEEN_HH

#include <dune/navier/global_defines.hh>
#include <dune/navier/problems.hh>

#include <dune/oseen/modelinterface.hh>
#include <dune/oseen/pass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/navier/thetascheme_traits.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/navier/fractionaldatawriter.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <cmath>

namespace Dune {
namespace Oseen {

template <class TimeProviderType, class GridImp, template <class, class> class ForceFuntionType,
          template <class, class> class AnalyticalDirichletDataImp, int gridDim, int sigmaOrder,
          int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder,
          template <class, class, int, template <class> class BaseFunctionStorageImp = Dune::CachingStorage>
          class GalerkinSpaceImp = Dune::DiscontinuousGalerkinSpace>
class DiscreteModelTraits {
public:
  //! using DGAdaptiveLeafGridPart is mandated by DUNE-FEM, but not in any way checked...
  typedef Dune::DGAdaptiveLeafGridPart<GridImp> GridPartType;
  //! for CRTP trick
  typedef DiscreteOseenModelDefault<DiscreteModelTraits> DiscreteModelType;

  //! we use caching quadratures for the entities
  typedef Dune::CachingQuadrature<GridPartType, 0> VolumeQuadratureType;

  //! we use caching quadratures for the faces
  typedef Dune::CachingQuadrature<GridPartType, 1> FaceQuadratureType;

  //! polynomial order for the discrete sigma function space
  static const int sigmaSpaceOrder = sigmaOrder;
  //! polynomial order for the discrete velocity function space
  static const int velocitySpaceOrder = velocityOrder;
  //! polynomial order for the discrete pressure function space
  static const int pressureSpaceOrder = pressureOrder;

  //    private:

  //! function space type for the velocity
  typedef Dune::FunctionSpace<double, double, gridDim, gridDim> VelocityFunctionSpaceType;

  //! discrete function space type for the velocity
  typedef GalerkinSpaceImp<VelocityFunctionSpaceType, GridPartType, velocitySpaceOrder>
  DiscreteVelocityFunctionSpaceType;

  //! function space type for the pressure
  typedef Dune::FunctionSpace<double, double, gridDim, 1> PressureFunctionSpaceType;

  //! discrete function space type for the pressure
  typedef GalerkinSpaceImp<PressureFunctionSpaceType, GridPartType, pressureSpaceOrder>
  DiscretePressureFunctionSpaceType;

public:

  //! discrete function space wrapper type
  typedef Dune::DiscreteOseenFunctionSpaceWrapper<Dune::DiscreteOseenFunctionSpaceWrapperTraits<
      DiscreteVelocityFunctionSpaceType, DiscretePressureFunctionSpaceType>> DiscreteOseenFunctionSpaceWrapperType;

  //! discrete function type for the velocity
  typedef Dune::AdaptiveDiscreteFunction<
      typename DiscreteOseenFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType> DiscreteVelocityFunctionType;

  //! discrete function type for the pressure
  typedef Dune::AdaptiveDiscreteFunction<
      typename DiscreteOseenFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType> DiscretePressureFunctionType;

public:

  //! discrete function wrapper type
  typedef Dune::DiscreteOseenFunctionWrapper<Dune::DiscreteOseenFunctionWrapperTraits<
      DiscreteOseenFunctionSpaceWrapperType, DiscreteVelocityFunctionType, DiscretePressureFunctionType>>
  DiscreteOseenFunctionWrapperType;

  //! function space type for sigma
  typedef Dune::MatrixFunctionSpace<double, double, gridDim, gridDim, gridDim> SigmaFunctionSpaceType;

  //! discrete function space type for sigma
  typedef GalerkinSpaceImp<SigmaFunctionSpaceType, GridPartType, sigmaSpaceOrder> DiscreteSigmaFunctionSpaceType;

public:

  //! discrete function type for sigma
  typedef Dune::AdaptiveDiscreteFunction<DiscreteSigmaFunctionSpaceType> DiscreteSigmaFunctionType;

  //! function type for the analytical force
  typedef ForceFuntionType<VelocityFunctionSpaceType, TimeProviderType> AnalyticalForceFunctionType;

  typedef AnalyticalForceFunctionType AnalyticalForceType;

  typedef AnalyticalDirichletDataImp<VelocityFunctionSpaceType, TimeProviderType> AnalyticalDirichletDataType;

  typedef DiscreteVelocityFunctionType ExtraDataFunctionType;
  /**
   *  \name   types needed for the pass
   *  \{
   **/
  //! return type of the pass
  typedef DiscreteOseenFunctionWrapperType DestinationType;
  /**
   *  \}
   **/
};

template <class CommunicatorImp, class GridImp, int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder,
          int pressureOrder = sigmaOrder>
struct Traits {
  typedef Traits<CommunicatorImp, GridImp, gridDim, sigmaOrder, velocityOrder, pressureOrder> ThisType;
  //! using DGAdaptiveLeafGridPart is mandated by DUNE-FEM, but not in any way checked...
  typedef Dune::DGAdaptiveLeafGridPart<GridImp> GridPartType;
  typedef Dune::NavierStokes::ThetaSchemeDescription<1> SchemeDescriptionType;
  typedef Dune::NavierStokes::FractionalTimeProvider<SchemeDescriptionType, CommunicatorImp> TimeProviderType;

  typedef DiscreteModelTraits<TimeProviderType, GridType, OSEEN_DATA_NAMESPACE::Force,
                              OSEEN_DATA_NAMESPACE::DirichletData, gridDim, sigmaOrder, velocityOrder,
                              pressureOrder> OseenModelTraits;
  typedef Dune::DiscreteOseenModelDefault<OseenModelTraits> OseenModelType;
  typedef Dune::OseenPass<OseenModelType> OseenPassType;

  typedef OSEEN_DATA_NAMESPACE::Pressure<typename OseenModelTraits::PressureFunctionSpaceType, TimeProviderType>
  ExactPressureType;
  typedef OSEEN_DATA_NAMESPACE::Velocity<typename OseenModelTraits::VelocityFunctionSpaceType, TimeProviderType>
  ExactVelocityType;
  typedef OSEEN_DATA_NAMESPACE::VelocityConvection<typename OseenModelTraits::VelocityFunctionSpaceType,
                                                   TimeProviderType> ConvectionType;
  typedef Dune::NavierStokes::ExactSolution<ThisType> ExactSolutionType;
  typedef typename OseenModelType::DiscreteOseenFunctionWrapperType DiscreteOseenFunctionWrapperType;
  typedef typename OseenModelType::DiscreteOseenFunctionSpaceWrapperType DiscreteOseenFunctionSpaceWrapperType;
  typedef typename DiscreteOseenFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType::FunctionSpaceType
  PressureFunctionSpaceType;
  typedef typename DiscreteOseenFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType::FunctionSpaceType
  VelocityFunctionSpaceType;
};

} // end namespace Oseen
} // end namespace Dune

#endif // OSEEN_HH

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
