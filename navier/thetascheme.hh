#ifndef METADATA_HH
#define METADATA_HH

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/navier/nonlinear/models.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/stuff/datawriter.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>

namespace Dune {
	namespace NavierStokes {
		template <	class CommunicatorImp,
					class GridPartImp,
					template < class > class AnalyticalForceImp,
					template < class > class AnalyticalDirichletDataImp,
					template < class,class > class ExactPressureImp,
					template < class,class > class ExactVelocityImp,
					int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
		struct ThetaSchemeTraits {
			typedef GridPartImp
				GridPartType;
			typedef FractionalTimeProvider<CommunicatorImp>
				TimeProviderType;

			typedef StokesStep::DiscreteStokesModelTraits<
						TimeProviderType,
						GridPartType,
						AnalyticalForceImp,
						AnalyticalDirichletDataImp,
						gridDim,
						sigmaOrder,
						velocityOrder,
						pressureOrder >
					StokesModelTraits;
			typedef Dune::DiscreteStokesModelDefault< StokesModelTraits >
				StokesModelType;
			typedef typename StokesModelTraits::DiscreteStokesFunctionSpaceWrapperType
				DiscreteStokesFunctionSpaceWrapperType;

			typedef typename StokesModelTraits::DiscreteStokesFunctionWrapperType
				DiscreteStokesFunctionWrapperType;
			typedef typename StokesModelTraits::AnalyticalForceFunctionType
				AnalyticalForceType;
			typedef typename StokesModelTraits::AnalyticalForceAdapterType
				StokesAnalyticalForceAdapterType;
			typedef typename StokesModelTraits::AnalyticalDirichletDataType
				AnalyticalDirichletDataType;

			typedef Dune::StartPass< DiscreteStokesFunctionWrapperType, -1 >
				StokesStartPassType;
			typedef Dune::StokesPass< StokesModelType, StokesStartPassType, 0 >
				StokesPassType;

			typedef CommunicatorImp
				CommunicatorType;

			typedef ExactPressureImp< typename StokesModelTraits::PressureFunctionSpaceType,
									  TimeProviderType >
				ExactPressureType;
			typedef ExactVelocityImp< typename StokesModelTraits::VelocityFunctionSpaceType,
									  TimeProviderType >
				ExactVelocityType;
		};

		template < class T1, class T2 >
		struct TupleSerializer {
			typedef Dune::Tuple<	const typename T1::DiscreteVelocityFunctionType*,
									const typename T1::DiscretePressureFunctionType*,
									const typename T2::DiscreteVelocityFunctionType*,
									const typename T2::DiscretePressureFunctionType* >
				TupleType;

			static TupleType& getTuple( T1& t1,
										T2& t2 )
			{
				static TupleType t( &(t1.discreteVelocity()),
									&(t1.discretePressure()),
									&(t2.discreteVelocity()),
									&(t2.discretePressure()) );
				return t;
			}
		};

		template < class TraitsImp >
		class ThetaScheme {
			protected:
				typedef TraitsImp
					Traits;
				typedef typename Traits::CommunicatorType
					CommunicatorType;
				typedef ExactSolution<Traits>
					ExactSolutionType;
				typedef TupleSerializer< typename Traits::DiscreteStokesFunctionWrapperType,
										 ExactSolutionType >
					TupleSerializerType;
				typedef typename TupleSerializerType::TupleType
					OutputTupleType;
				typedef TimeAwareDataWriter<	typename Traits::TimeProviderType,
												typename Traits::GridPartType::GridType,
												OutputTupleType >
					DataWriterType;
				typedef typename Traits::DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
					DiscreteVelocityFunctionType;
				typedef typename Traits::DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
					DiscretePressureFunctionType;

				typename Traits::GridPartType gridPart_;
				const double theta_;
				const double operator_weight_alpha_;
				const double operator_weight_beta_;
				CommunicatorType& communicator_;
				typename Traits::TimeProviderType timeprovider_;
				typename Traits::DiscreteStokesFunctionSpaceWrapperType functionSpaceWrapper_;
				typename Traits::DiscreteStokesFunctionWrapperType currentFunctions_;
				typename Traits::DiscreteStokesFunctionWrapperType nextFunctions_;
				ExactSolutionType exactSolution_;
				DataWriterType dataWriter_;

			public:
				ThetaScheme( typename Traits::GridPartType gridPart,
							 const double theta = 1 - std::pow( 2.0, -1/2.0 ),
							 CommunicatorType comm = Dune::MPIManager::helper().getCommunicator()
						)
					: gridPart_( gridPart ),
					theta_(theta),
					operator_weight_alpha_( ( 1-2*theta_ ) / ( 1-theta_ ) ),
					operator_weight_beta_( 1 - operator_weight_alpha_ ),
					communicator_( comm ),
					timeprovider_( theta_,operator_weight_alpha_,operator_weight_beta_, communicator_ ),
					functionSpaceWrapper_( gridPart_ ),
					currentFunctions_(  "current_",
										functionSpaceWrapper_,
										gridPart_ ),
					nextFunctions_(  "next_",
									functionSpaceWrapper_,
									gridPart_ ),
					exactSolution_( timeprovider_,
									gridPart_,
									functionSpaceWrapper_ ),
					dataWriter_( timeprovider_,
								 gridPart_.grid(),
								 TupleSerializerType::getTuple(
										 nextFunctions_,
										 exactSolution_ )
								)
				{}

				void nextStep( const int step )
				{
					dataWriter_.write();
					currentFunctions_.assign( nextFunctions_ );
					nextFunctions_ .clear();
					timeprovider_.nextFractional();
					std::cout << "current time (substep " << step << "): " << timeprovider_.subTime() << std::endl;
				}

				void stokesStep( const typename Traits::AnalyticalForceType& force,
								 const double stokes_viscosity,
								 const double quasi_stokes_alpha,
								 const double beta_qout_re )
				{
					typename Traits::StokesAnalyticalForceAdapterType stokesForce( timeprovider_,
																				   currentFunctions_.discreteVelocity(),
																				   force,
																				   beta_qout_re );
					typename Traits::AnalyticalDirichletDataType stokesDirichletData =
							Traits::StokesModelTraits::AnalyticalDirichletDataTraitsImplementation
											::getInstance( timeprovider_,
														   functionSpaceWrapper_ );
					typename Traits::StokesModelType
							stokesModel( Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients() ,
										stokesForce,
										stokesDirichletData,
										stokes_viscosity,
										quasi_stokes_alpha );
					typename Traits::StokesStartPassType stokesStartPass;
					typename Traits::StokesPassType stokesPass( stokesStartPass,
											stokesModel,
											gridPart_,
											functionSpaceWrapper_ );
					stokesPass.apply( currentFunctions_, nextFunctions_ );
				}

				void run()
				{
					//initial flow field at t = 0
					currentFunctions_.projectInto( exactSolution_.exactVelocity(), exactSolution_.exactPressure() );

					//constants
					const double viscosity				= 1;
					const double quasi_stokes_alpha		= 1 / ( theta_ * timeprovider_.deltaT() );
					const double reynolds				= 1 / viscosity;//not really, but meh
					const double stokes_viscosity		= operator_weight_alpha_ / reynolds;
					const double beta_qout_re			= operator_weight_beta_ / reynolds;
					const int verbose					= 1;
					const typename Traits::AnalyticalForceType force ( 1.0 /*visc*/,
																 currentFunctions_.discreteVelocity().space() );

					for( timeprovider_.init( timeprovider_.deltaT() ); timeprovider_.time() < timeprovider_.endTime(); )
					{
						std::cout << "current time (substep " << 0 << "): " << timeprovider_.subTime() << std::endl;
						//stokes step A
						stokesStep( force, stokes_viscosity, quasi_stokes_alpha, beta_qout_re );

						nextStep( 1 );
						//Nonlinear step
						{

							typedef NonlinearStep::ForceAdapterFunction<	typename Traits::TimeProviderType,
																			typename Traits::AnalyticalForceType,
																			DiscreteVelocityFunctionType,
																			DiscretePressureFunctionType >
									NonlinearForceAdapterFunctionType;

							NonlinearForceAdapterFunctionType nonlinearForce( timeprovider_,
																			  currentFunctions_.discreteVelocity(),
																			  currentFunctions_.discretePressure(),
																			  force,
																			  operator_weight_alpha_ / reynolds );
							typedef NonlinearStep::Traits<	typename Traits::GridPartType,
															typename Traits::DiscreteStokesFunctionWrapperType,
															NonlinearForceAdapterFunctionType >
								NonlinearTraits;
							typename NonlinearTraits::InitialDataType problem_( currentFunctions_.discreteVelocity() );
							typename NonlinearTraits::ModelType model_( problem_,
																		currentFunctions_,
																		nonlinearForce );
							// Initial flux for advection discretization (UpwindFlux)
							typename NonlinearTraits::FluxType convectionFlux_( model_ );
							typename NonlinearTraits::DgType dg_( gridPart_.grid(),
																  convectionFlux_ );
							typename NonlinearTraits::ODEType ode( dg_,
																   timeprovider_,
																   1,
																   verbose );

							typename NonlinearTraits:: DgType :: SpaceType  sp(gridPart_);
							//this is a non-compatible df type, we'll have to project the solution back ito our space afterwards
							typedef typename NonlinearTraits:: DgType :: DestinationType
								NonLinearVelocityType;
							NonLinearVelocityType nonlinear_velocity( "de", sp );

							//set starttime to current, endtime to next
							ode.initialize( nonlinear_velocity );
							ode.solve( nonlinear_velocity );
							//recast
							Dune::BetterL2Projection::project( nonlinear_velocity, currentFunctions_.discreteVelocity() );

						}
						nextStep( 2 );
						//stokes step B
						stokesStep( force, stokes_viscosity, quasi_stokes_alpha, beta_qout_re  );

						nextStep( 3 );

						//error calc
						{
							exactSolution_.project();
							DiscretePressureFunctionType errorFunc_pressure_("",currentFunctions_.discretePressure().space());
							DiscreteVelocityFunctionType errorFunc_velocity_("",currentFunctions_.discreteVelocity().space());
							errorFunc_pressure_.assign( exactSolution_.discretePressure() );
							errorFunc_pressure_ -= currentFunctions_.discretePressure();
							errorFunc_velocity_.assign( exactSolution_.discreteVelocity() );
							errorFunc_velocity_ -= currentFunctions_.discreteVelocity();

							Dune::L2Norm< typename Traits::GridPartType > l2_Error( gridPart_ );

							const double l2_error_pressure_ = l2_Error.norm( errorFunc_pressure_ );
							const double l2_error_velocity_ = l2_Error.norm( errorFunc_velocity_ );

							Logger().Info().Resume();
							Logger().Info() << "L2-Error Pressure: " << std::setw(8) << l2_error_pressure_ << "\n"
											<< "L2-Error Velocity: " << std::setw(8) << l2_error_velocity_ << std::endl;
						}
					}
				}

		};
	}//end namespace NavierStokes
}//end namespace Dune

#endif // METADATA_HH
