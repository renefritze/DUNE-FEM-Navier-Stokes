#ifndef STOKESTRAITS_HH
#define STOKESTRAITS_HH

#include <dune/oseen/modelinterface.hh>
#include <dune/navier/rhsadapter.hh>
#include <dune/navier/weighed_function.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

namespace Dune {
namespace NavierStokes {
namespace NonlinearStep {
template <class TimeProviderType, class GridImp, template <class, class> class AnalyticalForceFunctionType,
          template <class, class, class, class> class ForceAdatperTemplateType,
          template <class, class> class AnalyticalDirichletDataImp, class ThetaValueArrayType, int gridDim,
          int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder,
          template <class, class, int, template <class> class BaseFunctionStorageImp = Dune::CachingStorage>
          class GalerkinSpaceImp = Dune::DiscontinuousGalerkinSpace>
class DiscreteOseenModelTraits {
public:
  //! using DGAdaptiveLeafGridPart is mandated by DUNE-FEM, but not in any way checked...
  typedef Dune::DGAdaptiveLeafGridPart<GridImp> GridPartType;
  //! for CRTP trick
  typedef DiscreteOseenModelDefault<DiscreteOseenModelTraits> DiscreteModelType;

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

private:
  //! function space type for sigma
  typedef Dune::MatrixFunctionSpace<double, double, gridDim, gridDim, gridDim> SigmaFunctionSpaceType;

public:
  //! discrete function space type for sigma
  typedef GalerkinSpaceImp<SigmaFunctionSpaceType, GridPartType, sigmaSpaceOrder> DiscreteSigmaFunctionSpaceType;

#if STOKES_USE_ISTL
  //! discrete function type for the velocity
  typedef Dune::BlockVectorDiscreteFunction<
      typename DiscreteOseenFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType> DiscreteVelocityFunctionType;

  //! discrete function type for the pressure
  typedef Dune::BlockVectorDiscreteFunction<
      typename DiscreteOseenFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType> DiscretePressureFunctionType;

  //! discrete function type for sigma
  typedef Dune::BlockVectorDiscreteFunction<DiscreteSigmaFunctionSpaceType> DiscreteSigmaFunctionType;
#else
  //! discrete function type for the velocity
  typedef Dune::AdaptiveDiscreteFunction<
      typename DiscreteOseenFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType> DiscreteVelocityFunctionType;

  //! discrete function type for the pressure
  typedef Dune::AdaptiveDiscreteFunction<
      typename DiscreteOseenFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType> DiscretePressureFunctionType;

  //! discrete function type for sigma
  typedef Dune::AdaptiveDiscreteFunction<DiscreteSigmaFunctionSpaceType> DiscreteSigmaFunctionType;
#endif

  //! discrete function wrapper type
  typedef Dune::DiscreteOseenFunctionWrapper<Dune::DiscreteOseenFunctionWrapperTraits<
      DiscreteOseenFunctionSpaceWrapperType, DiscreteVelocityFunctionType, DiscretePressureFunctionType>>
  DiscreteOseenFunctionWrapperType;

public:
  //! function type for the analytical force
  typedef AnalyticalForceFunctionType<VelocityFunctionSpaceType, TimeProviderType> RealAnalyticalForceType;

  typedef ForceAdatperTemplateType<TimeProviderType, RealAnalyticalForceType, DiscreteVelocityFunctionType,
                                   ThetaValueArrayType> ForceAdatperType;

  typedef ForceAdatperType AnalyticalForceType;

  //! function type for the analytical dirichlet data
  //					typedef typename Dune::NavierStokes::DirichletAdapterFunctionTraits< AnalyticalDirichletDataImp,
  //TimeProviderType >
  //										::template Implementation<VelocityFunctionSpaceType,GridPartType >
  //							AnalyticalDirichletDataTraitsImplementation;
  typedef AnalyticalDirichletDataImp<VelocityFunctionSpaceType, TimeProviderType> AnalyticalDirichletDataFunctionType;
  typedef WeighedIntersectionFunction<VelocityFunctionSpaceType, TimeProviderType, AnalyticalDirichletDataFunctionType>
  AnalyticalDirichletDataType;

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
} // namespace Nonlinear
} // end namespace NavierStokes
} // end namespace Dune

#endif // STOKESTRAITS_HH

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
