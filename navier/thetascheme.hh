#ifndef METADATA_HH
#define METADATA_HH

#include <dune/navier/exactsolution.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/stuff/datawriter.hh>
#include <dune/stuff/tuple.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/stuff/error.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/profiler.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>
#include <boost/scoped_ptr.hpp>
#include <algorithm>

#include <dune/navier/thetascheme_traits.hh>

namespace Dune {
	namespace NavierStokes {
		template < class TraitsImp >
		class ThetaScheme {
			protected:
				typedef TraitsImp
					Traits;
				typedef typename Traits::CommunicatorType
					CommunicatorType;
				typedef typename Traits::ExactSolutionType
					ExactSolutionType;
				typedef Stuff::TupleSerializer<	typename Traits::DiscreteStokesFunctionWrapperType,
											typename Traits::DiscreteStokesFunctionWrapperType,
											ExactSolutionType,
											typename Traits::DiscreteStokesFunctionWrapperType>
					TupleSerializerType1;
				typedef typename TupleSerializerType1::TupleType
					OutputTupleType1;
				typedef TimeAwareDataWriter<	typename Traits::TimeProviderType,
												typename Traits::GridPartType::GridType,
												OutputTupleType1 >
					DataWriterType1;
				typedef CheckPointer< typename Traits::GridPartType::GridType,
									  OutputTupleType1 >
					CheckPointerType;
				typedef Stuff::TupleSerializer<	typename Traits::DiscreteStokesFunctionWrapperType >
					TupleSerializerType2;
				typedef typename TupleSerializerType2::TupleType
					OutputTupleType2;
				typedef TimeAwareDataWriter<	typename Traits::TimeProviderType,
												typename Traits::GridPartType::GridType,
												OutputTupleType2 >
					DataWriterType2;
				typedef typename Traits::DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
					DiscreteVelocityFunctionType;
				typedef typename Traits::DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
					DiscretePressureFunctionType;

				mutable typename Traits::GridPartType gridPart_;
				const typename Traits::ThetaSchemeDescriptionType& scheme_params_;

		protected:
				CommunicatorType& communicator_;
				typename Traits::TimeProviderType timeprovider_;
				typename Traits::DiscreteStokesFunctionSpaceWrapperType functionSpaceWrapper_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType currentFunctions_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType nextFunctions_;
				typename Traits::DiscreteStokesFunctionWrapperType errorFunctions_;
				ExactSolutionType exactSolution_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType dummyFunctions_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType updateFunctions_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType rhsFunctions_;
				OutputTupleType1& data_tuple_1;
				DataWriterType1 dataWriter1_;
				CheckPointerType check_pointer_;
				DataWriterType2 dataWriter2_;
				const typename Traits::OseenPassType::Traits::DiscreteSigmaFunctionSpaceType sigma_space_;
				mutable typename Traits::OseenPassType::RhsDatacontainer rhsDatacontainer_;

				typedef Stuff::L2Error< typename Traits::GridPartType >
					L2ErrorType;
				L2ErrorType l2Error_;


			public:
				const double viscosity_;
				const double d_t_;
				const double reynolds_;
				double current_max_gridwidth_;

			public:
				ThetaScheme( typename Traits::GridPartType gridPart,
							 const typename Traits::ThetaSchemeDescriptionType& scheme_params,
							 CommunicatorType comm			= Dune::MPIManager::helper().getCommunicator()
						)
					: gridPart_( gridPart ),
					scheme_params_( scheme_params ),
					communicator_( comm ),
					timeprovider_( scheme_params_, communicator_ ),
					functionSpaceWrapper_( gridPart_ ),
					currentFunctions_(  "current_",
										functionSpaceWrapper_,
										gridPart_ ),
					nextFunctions_(  "next_",
									functionSpaceWrapper_,
									gridPart_ ),
					errorFunctions_(  "error_",
									functionSpaceWrapper_,
									gridPart_ ),
					exactSolution_( timeprovider_,
									gridPart_,
									functionSpaceWrapper_ ),
					dummyFunctions_("dummy",
									functionSpaceWrapper_,
									gridPart_ ),
					updateFunctions_("updates",
									  functionSpaceWrapper_,
									  gridPart_ ),
					rhsFunctions_("rhs-adapter",
								 functionSpaceWrapper_,
								 gridPart_ ),
					data_tuple_1( TupleSerializerType1::getTuple(
							  currentFunctions_,
							  errorFunctions_,
							  exactSolution_,
							  dummyFunctions_) ),
					dataWriter1_( timeprovider_,
								 gridPart_.grid(),
								 data_tuple_1
								),
					check_pointer_(	gridPart_.grid(),
									"myGridName",
									data_tuple_1,
									timeprovider_.startTime(),
									timeprovider_.endTime(),
									static_cast<const LoadBalancerInterface*>( 0 )
								  ),
					dataWriter2_( timeprovider_,
								 gridPart_.grid(),
								 TupleSerializerType2::getTuple(
										 updateFunctions_,
										 rhsFunctions_)
								),
					sigma_space_( gridPart_ ),
					rhsDatacontainer_( currentFunctions_.discreteVelocity().space(), sigma_space_ ),
					l2Error_( gridPart ),
					viscosity_( Parameters().getParam( "viscosity", 1.0, Dune::ValidateNotLess<double>(0.0) ) ),
					d_t_( timeprovider_.deltaT() ),
					reynolds_( 1.0 / viscosity_ ),
					current_max_gridwidth_( Dune::GridWidth::calcGridWidth( gridPart_ ) )
				{
					Logger().Info() << scheme_params_;
				}

				void nextStep( const int step, RunInfo& info )
				{
					current_max_gridwidth_ = Dune::GridWidth::calcGridWidth( gridPart_ );
					currentFunctions_.assign( nextFunctions_ );
					exactSolution_.project();
					const bool last_substep = ( step == ( Traits::ThetaSchemeDescriptionType::numberOfSteps_ -1) );

					//error calc
					if ( Parameters().getParam( "calculate_errors", true ) ) {
						Profiler::ScopedTiming error_time("error_calc");

						errorFunctions_.discretePressure().assign( exactSolution_.discretePressure() );
						errorFunctions_.discretePressure() -= currentFunctions_.discretePressure();
						errorFunctions_.discreteVelocity().assign( exactSolution_.discreteVelocity() );
						errorFunctions_.discreteVelocity() -= currentFunctions_.discreteVelocity();

						double meanPressure_exact = Stuff::integralAndVolume( exactSolution_.exactPressure(), currentFunctions_.discretePressure().space() ).first;
						double meanPressure_discrete = Stuff::integralAndVolume( currentFunctions_.discretePressure(), currentFunctions_.discretePressure().space() ).first;

						Dune::L2Norm< typename Traits::GridPartType > l2_Error( gridPart_ );
						Dune::H1Norm< typename Traits::GridPartType > h1_Error( gridPart_ );

//						if ( Parameters().getParam( "error_scaling", false ) ) {
//								const double scale		= 1 / std::sqrt( viscosity_ );
//								errorFunctions_.discretePressure() *= scale;
//								errorFunctions_.discreteVelocity() *= scale;
//						}

						const double l2_error_pressure_				= l2_Error.norm( errorFunctions_.discretePressure() );
						const double l2_error_velocity_				= l2_Error.norm( errorFunctions_.discreteVelocity() );
						const double h1_error_pressure_				= h1_Error.norm( errorFunctions_.discretePressure() );
						const double h1_error_velocity_				= h1_Error.norm( errorFunctions_.discreteVelocity() );
						const double relative_l2_error_pressure_	= l2_error_pressure_ / l2_Error.norm( exactSolution_.discretePressure() );
						const double relative_l2_error_velocity_	= l2_error_velocity_ / l2_Error.norm( exactSolution_.discreteVelocity() );
						const double relative_h1_error_velocity_	= h1_error_velocity_ / h1_Error.norm( exactSolution_.discreteVelocity() );
						std::vector<double> error_vector;
						error_vector.push_back( l2_error_velocity_ );
						error_vector.push_back( l2_error_pressure_ );
						std::vector<double> h1_error_vector;
						h1_error_vector.push_back( h1_error_velocity_ );
						h1_error_vector.push_back( h1_error_pressure_ );

					#ifdef NDEBUG
						if ( last_substep ) //no need to be so verbose otherwise
					#endif
						{
							Logger().Info().Resume();
							Logger().Info() << boost::format ("L2-Error Pressure (abs|rel): %e | %e \t Velocity (abs|rel): %e | %e \n")
												% l2_error_pressure_ % relative_l2_error_pressure_
												% l2_error_velocity_ % relative_l2_error_velocity_
										#ifndef NDEBUG
											<< boost::format ("H1-Error Velocity (abs|rel): %e | %e") % h1_error_velocity_ % relative_h1_error_velocity_
											<< "Mean pressure (exact|discrete): " << meanPressure_exact << " | " << meanPressure_discrete
										#endif
											<< std::endl;
						}
						const double max_l2_error = 1e4;
						if ( l2_error_velocity_ > max_l2_error )
							DUNE_THROW(MathError, "Aborted, L2 error above " << max_l2_error );
						if ( std::isnan( l2_error_velocity_ ) || std::isnan( l2_error_pressure_ )  )
							throw Stuff::singlerun_abort_exception("L2 error is Nan");
						info.L2Errors		= error_vector;
						info.H1Errors		= h1_error_vector;
					}
					//end error calc

					if ( last_substep ) {
						typedef Dune::StabilizationCoefficients::ValueType
							Pair;
						Dune::StabilizationCoefficients stabil_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();

						info.codim0			= gridPart_.grid().size( 0 );
						info.grid_width		= current_max_gridwidth_;
						info.run_time		= profiler().GetTiming( "full_step" );
						info.delta_t		= timeprovider_.deltaT();
						info.current_time	= timeprovider_.subTime();
						info.viscosity		= viscosity_;
						info.reynolds		= reynolds_;

						info.c11			= Pair( stabil_coeff.Power( "C11" ), stabil_coeff.Factor( "C11" ) );
						info.c12			= Pair( stabil_coeff.Power( "C12" ), stabil_coeff.Factor( "C12" ) );
						info.d11			= Pair( stabil_coeff.Power( "D11" ), stabil_coeff.Factor( "D11" ) );
						info.d12			= Pair( stabil_coeff.Power( "D12" ), stabil_coeff.Factor( "D12" ) );
						info.bfg			= Parameters().getParam( "do-bfg", true );
						info.gridname		= gridPart_.grid().name();
						info.refine_level	= Parameters().getParam( "minref", 0, Dune::ValidateNotLess<int>(0) );

						info.polorder_pressure	= Traits::OseenModelTraits::pressureSpaceOrder;
						info.polorder_sigma		= Traits::OseenModelTraits::sigmaSpaceOrder;
						info.polorder_velocity	= Traits::OseenModelTraits::velocitySpaceOrder;

						info.solver_accuracy		= Parameters().getParam( "absLimit", 1e-4 );
						info.inner_solver_accuracy	= Parameters().getParam( "inner_absLimit", 1e-4 );
						info.bfg_tau				= Parameters().getParam( "bfg-tau", 0.1 );

						info.problemIdentifier	= TESTCASE_NAME;
						info.algo_id			= scheme_params_.algo_id;
						info.extra_info			= (boost::format("%s on %s") % COMMIT % std::getenv("HOSTNAME") ).str();

						Logger().Info() << boost::format("current time (substep %d ): %f (%f)\n")
											   % step
											   % timeprovider_.subTime()
											   % timeprovider_.previousSubTime();
					}

					if ( last_substep || !Parameters().getParam( "write_fulltimestep_only", false ) )
						writeData();
					timeprovider_.nextFractional();
				}


				void Init()
				{
					timeprovider_.init( d_t_ );
					//initial flow field at t = 0
					exactSolution_.project();
					currentFunctions_.assign( exactSolution_ );
					nextFunctions_.assign( exactSolution_ );
					writeData();
					//set current time to t_0 + dt_k
					timeprovider_.nextFractional();
				}

				RunInfoTimeMap run()
				{
					RunInfoTimeMap runInfoMap;
					Init();

					for( ;timeprovider_.time() <= timeprovider_.endTime(); )
					{
						try {
							RunInfo info = full_timestep();
							const double real_time = timeprovider_.subTime();
							nextStep( Traits::substep_count -1 , info );
							timeprovider_.printRemainderEstimate( Logger().Info() );
							runInfoMap[real_time] = info;
						}
						catch ( Stuff::singlerun_abort_exception& e ) {
							Logger().Err() << e.what() << std::endl;
							return runInfoMap;
						}
					}
					assert( runInfoMap.size() > 0 );
					return runInfoMap;
				}

				RunInfo full_timestep()
				{
					Profiler::ScopedTiming fullstep_time("full_step");
					RunInfo info;
					for ( int i=0; i < Traits::substep_count; ++i )
					{
						const double dt_k = scheme_params_.step_sizes_[i];
						substep( dt_k, scheme_params_.thetas_[i] );
						if ( i != Traits::substep_count - 1 )
							//the last step increase is done after one call level up
							nextStep( i, info );
					}
					return info;
				}

				void substep( const double dt_k, const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& theta_values )
				{
					//build rhs
					const bool first_step = timeprovider_.timeStep() <= 1;
					const typename Traits::AnalyticalForceType force ( viscosity_,
																 currentFunctions_.discreteVelocity().space() );
					if ( Parameters().getParam( "rhs_cheat", false ) )
						cheatRHS();

					boost::scoped_ptr< typename Traits::OseenForceAdapterFunctionType >
							ptr_oseenForce( first_step //in our very first step no previous computed data is avail. in rhs_container
												? new typename Traits::OseenForceAdapterFunctionType (	timeprovider_,
																										currentFunctions_.discreteVelocity(),
																										force,
																										reynolds_,
																										theta_values )
												: new typename Traits::OseenForceAdapterFunctionType (	timeprovider_,
																										currentFunctions_.discreteVelocity(),
																										force,
																										reynolds_,
																										theta_values,
																										rhsDatacontainer_ )
											);
					rhsFunctions_.discreteVelocity().assign( *ptr_oseenForce );
					typename Traits::StokesStartPassType stokesStartPass;
					typename Traits::AnalyticalDirichletDataType oseenDirichletData =
							Traits::OseenModelTraits::AnalyticalDirichletDataTraitsImplementation
											::getInstance( timeprovider_,
														   functionSpaceWrapper_ );

					unsigned int oseen_iterations = Parameters().getParam( "oseen_iterations", (unsigned int)(1), ValidateGreater<unsigned int>( 0 ) );
					const double dt_n = timeprovider_.deltaT();
					typename L2ErrorType::Errors old_error_velocity
							= l2Error_.get( currentFunctions().discreteVelocity(), exactSolution_.discreteVelocity() );
					typename L2ErrorType::Errors old_error_pressure
							= l2Error_.get( currentFunctions().discretePressure(), exactSolution_.discretePressure() );
					double velocity_error_reduction = 1.0;
					double pressure_error_reduction = 1.0;
					unsigned int i = 0;
					do
					{
						const double last_velocity_error_reduction = velocity_error_reduction;
						const double last_pressure_error_reduction = pressure_error_reduction;
						typename Traits::OseenModelType
								oseenModel( Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients(),
											*ptr_oseenForce,
											oseenDirichletData,
											theta_values[0] * dt_n / reynolds_, /*viscosity*/
											1.0f, /*alpha*/
											dt_k,/*pressure_gradient_scale_factor*/
											theta_values[0] * dt_n /*convection_scale_factor*/
						                   );
						typename Traits::OseenPassType oseenPass( stokesStartPass,
												oseenModel,
												gridPart_,
												functionSpaceWrapper_,
												currentFunctions_.discreteVelocity() /*beta*/,
												true /*do_oseen_disc*/ );
						if ( timeprovider_.timeStep() <= 1 )
							oseenPass.printInfo();
						if ( Parameters().getParam( "silent_stokes", true ) )
							Logger().Info().Suspend( Logging::LogStream::default_suspend_priority + 10 );
						oseenPass.apply( currentFunctions_, nextFunctions_, &rhsDatacontainer_ );
						Logger().Info().Resume( Logging::LogStream::default_suspend_priority + 10 );

						{
							Profiler::ScopedTiming error_time("error_calc");
							typename L2ErrorType::Errors new_error_velocity
									= l2Error_.get( nextFunctions_.discreteVelocity(), exactSolution_.discreteVelocity() );
							typename L2ErrorType::Errors new_error_pressure
									= l2Error_.get( nextFunctions_.discretePressure(), exactSolution_.discretePressure() );
							velocity_error_reduction = old_error_velocity.absolute() / new_error_velocity.absolute();
							pressure_error_reduction = old_error_pressure.absolute() / new_error_pressure.absolute() ;
						}

						currentFunctions_.assign( nextFunctions_ );

						bool abort_loop = false;
						if ( ( ( pressure_error_reduction < 1.0 )
							  && ( velocity_error_reduction < 1.0 ) ) )
						{
							Logger().Info() << "Oseen iteration increased error, aborting.. -- ";
							abort_loop = true;
						}

						if ( ( pressure_error_reduction > 10.0 )
								|| ( velocity_error_reduction > 10.0 ) )
						{
							Logger().Info() << "Oseen iteration reduced error by factor 10, aborting.. -- ";
							abort_loop = true;
						}
						if (  ( ! ( ( last_pressure_error_reduction != pressure_error_reduction )
									|| ( last_velocity_error_reduction != velocity_error_reduction ) ) )
								|| ( pressure_error_reduction < Parameters().getParam( "min_error_reduction", 1.05 ) )
								|| ( velocity_error_reduction < Parameters().getParam( "min_error_reduction", 1.05 ) ) )
						{
							Logger().Info() << "Oseen iteration reduced no error, aborting.. -- ";
							abort_loop = true;
						}
						if ( abort_loop )
						{
							Logger().Info() << boost::format(" iteration %d, error reduction: pressure  %e | velocity %e")
																		   % i % pressure_error_reduction % velocity_error_reduction
																		<< std::endl;
							break;
						}
					} while ( i++ < oseen_iterations ) ;
				}

				void setUpdateFunctions() const
				{
					updateFunctions_.assign( nextFunctions_);
					updateFunctions_ -= currentFunctions_ ;
				}

				void writeData()
				{
					Profiler::ScopedTiming io_time("IO");
					dataWriter1_.write();
					dataWriter2_.write();
//					check_pointer_.write( timeprovider_.time(), timeprovider_.timeStep() );
				}

				typename Traits::OseenPassType::RhsDatacontainer& rhsDatacontainer()
				{
					return rhsDatacontainer_;
				}

				const ExactSolutionType& exactSolution() const
				{
					return exactSolution_;
				}

				const typename Traits::DiscreteStokesFunctionWrapperType& currentFunctions() const
				{
					return currentFunctions_;
				}

				const typename Traits::TimeProviderType& timeprovider() const
				{
					return timeprovider_;
				}

				void cheatRHS()
				{
					typedef typename DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType
						VelocityFunctionSpaceType;
					VelocityFunctionSpaceType continousVelocitySpace_;

					// ----
					typedef TESTING_NS::VelocityLaplace<	VelocityFunctionSpaceType,
																				typename Traits::TimeProviderType >
							VelocityLaplace;
					VelocityLaplace velocity_laplace( timeprovider_, continousVelocitySpace_ );
					Dune::BetterL2Projection
						::project( timeprovider_.previousSubTime(), velocity_laplace, rhsDatacontainer_.velocity_laplace );
					currentFunctions_.discreteVelocity().assign( exactSolution_.discreteVelocity() );
				}
		};
	}//end namespace NavierStokes
}//end namespace Dune

#endif // METADATA_HH
