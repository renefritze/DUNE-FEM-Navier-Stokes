#ifndef CONV_DIFF_HH
#define CONV_DIFF_HH

#include <dune/navier/global_defines.hh>

#ifndef CONVDIFF_DATA_NAMESPACE
//	#define CONVDIFF_DATA_NAMESPACE Dune::ConvDiff::AdapterFunctionsVectorial
#define CONVDIFF_DATA_NAMESPACE Dune::ConvDiff::TimeDisc
#endif
#define TESTING_NS CONVDIFF_DATA_NAMESPACE

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/misc/l2error.hh>

#include <dune/oseen/functionspacewrapper.hh>
#include <dune/oseen/modelinterface.hh>
#include <dune/oseen/pass.hh>
#include <dune/oseen/boundarydata.hh>

#include <dune/stuff/printing.hh>
#include <dune/stuff/femeoc.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/oseen/postprocessing.hh>
#include <dune/stuff/profiler.hh>
#include <dune/stuff/timeseries.hh>
#include <dune/stuff/signals.hh>
#include <dune/stuff/tuple.hh>
#include <dune/stuff/error.hh>
#include <dune/stuff/functionadapter.hh>

namespace Dune {
namespace ConvDiff {
namespace TestCase2D {
template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
   **/
  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha)
    , lambda_(Parameters().getParam("lambda", 0.0))
    , gamma_(Parameters().getParam("alpha", 0.0)) {}

  /**
   *  \brief  destructor
   *  doing nothing
   **/
  ~Force() {}

  /**
   *  \brief  evaluates the force
   *  \param  arg
   *          point to evaluate at
   *  \param  ret
   *          value of force at given point
   **/
  inline void evaluate(const double /*time*/, const DomainType& arg, RangeType& ret) const {
    const double lambda = lambda_;
    const double x = arg[0];
    const double y = arg[1];
    const double e_lambda_x = std::exp(lambda * x);
    const double cos = std::cos(2 * M_PI * y);
    const double sin = std::sin(2 * M_PI * y);
    const double lambda_square = lambda * lambda;
    RangeType u;
    VelocityEvaluate(lambda_, 0, arg, u);
    // convection
    //					  assert( false ); //M_2_PI == 2 / PI
    ret[0] = u[0] * (-lambda * e_lambda_x * cos) + u[1] * (-M_2_PI * e_lambda_x * sin);
    ret[1] = u[0] * ((lambda_square / M_2_PI) * e_lambda_x * sin) + u[1] * lambda * e_lambda_x * cos;
    // laplace
    ret[0] -= viscosity_ * (-lambda_square * e_lambda_x * cos);
    ret[1] -= viscosity_ * (-lambda * e_lambda_x * sin);
    // pressure grad
    ret[0] += lambda * std::exp(lambda * 2 * x);
    ret[1] += 0;
    // u
    ret[0] += gamma_ * u[0];
    ret[1] += gamma_ * u[1];
  }
  inline void evaluate(const DomainType& arg, RangeType& ret) const { evaluate(0, arg, ret); }

private:
  const double viscosity_;
  const double alpha_;
  const double lambda_;
  const double gamma_;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double lambda, const double time, const DomainType& arg, RangeType& ret) {
  const double x = arg[0];
  const double y = arg[1];
  const double e_lambda_x = std::exp(lambda * x);

  ret[0] = 1 - e_lambda_x * std::cos(2 * M_PI * y);
  ret[1] = (lambda / (2 * M_PI)) * e_lambda_x * std::sin(2 * M_PI * y);
}
template <class FunctionSpaceImp>
class Beta : public Fem::Function<FunctionSpaceImp, Beta<FunctionSpaceImp>> {
public:
  typedef Beta<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  Beta(double dummy, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
       const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~Beta() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

  /**
  * \brief  evaluates the dirichlet data
  * \param  arg
  *         point to evaluate at
  * \param  ret
  *         value of dirichlet boundary data at given point
  **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const { VelocityEvaluate(lambda_, 0, arg, ret); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double lambda_;
};
template <class FunctionSpaceImp>
class ExactConvection : public Fem::Function<FunctionSpaceImp, ExactConvection<FunctionSpaceImp>> {
public:
  typedef ExactConvection<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  ExactConvection(double dummy, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                  const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~ExactConvection() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    const double lambda = lambda_;
    const double x = arg[0];
    const double y = arg[1];
    const double e_lambda_x = std::exp(lambda * x);
    const double cos = std::cos(2 * M_PI * y);
    const double sin = std::sin(2 * M_PI * y);
    const double lambda_square = lambda * lambda;
    RangeType u;
    VelocityEvaluate(lambda_, time, arg, u);
    // convection
    //					  assert( false ); //M_2_PI == 2 / PI
    ret[0] = u[0] * (-lambda * e_lambda_x * cos) + u[1] * (-M_2_PI * e_lambda_x * sin);
    ret[1] = u[0] * ((lambda_square / M_2_PI) * e_lambda_x * sin) + u[1] * lambda * e_lambda_x * cos;
  }

  /**
  * \brief  evaluates the dirichlet data
  * \param  arg
  *         point to evaluate at
  * \param  ret
  *         value of dirichlet boundary data at given point
  **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const { VelocityEvaluate(lambda_, 0, arg, ret); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double lambda_;
};
/**
*  \brief  describes the dirichlet boundary data
*
*  \tparam DirichletTraitsImp
*          types like functionspace, range type, etc
*
*  \todo   extensive docu with latex
**/
template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

  /**
  * \brief  evaluates the dirichlet data
  * \param  arg
  *         point to evaluate at
  * \param  ret
  *         value of dirichlet boundary data at given point
  **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const { assert(false); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double lambda_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class Velocity : public TimeFunction<FunctionSpaceImp, Velocity<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Velocity<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  Velocity(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~Velocity() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

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

template <class FunctionSpaceImp, class TimeProviderImp>
class Pressure : public TimeFunction<FunctionSpaceImp, Pressure<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Pressure<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *
   *  doing nothing besides Base init
   **/
  Pressure(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0))
    , shift_(0.0) {}

  /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
  ~Pressure() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    Dune::CompileTimeChecker<(dim_ == 2)> Pressure_Unsuitable_WorldDim;
    const double x = arg[0];
    const double y = arg[1];
    const double e_2lambda_x = std::exp(2 * lambda_ * x);

    ret[0] = 0.5 * e_2lambda_x + shift_;
  }

  void setShift(const double shift) {
    shift_ = shift;
    Logger().Info() << "Set pressure shift to: " << shift_ << std::endl;
  }

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
  double shift_;
};

} // end namespace Oseem::TestCase2D

namespace TestCaseTaylor2D {

template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
   **/
  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  /**
   *  \brief  destructor
   *  doing nothing
   **/
  ~Force() {}

  /**
   *  \brief  evaluates the force
   *  \param  arg
   *          point to evaluate at
   *  \param  ret
   *          value of force at given point
   **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const {
    const double x = arg[0];
    const double y = arg[1];
    const double v = viscosity_;
    const double a = Parameters().getParam("alpha", 1.0);
    ;
    const double P = 2 * M_PI; // pi_factor;
    const double E = 1; // std::exp( -2 * std::pow( P, 2 ) * viscosity_ * time );
    const double F = 1; // std::exp( -4 * std::pow( P, 2 ) * viscosity_ * time );
    const double S_x = std::sin(P * x);
    const double S_y = std::sin(P * y);
    const double S_2x = std::sin(2 * P * x);
    const double S_2y = std::sin(2 * P * y);
    const double C_x = std::cos(P * x);
    const double C_y = std::cos(P * y);
    RangeType u, u_a;
    VelocityEvaluate(0, 0, arg, u);
    u_a = (u * a);
    ret = u_a;
    ret[0] = -C_x * E * P * (S_x * E + v * S_y * P);
    ret[1] = -C_y * E * P * (S_y * E - v * S_x * P);
    //					  ret = RangeType(0);
  }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double lambda, const double time, const DomainType& arg, RangeType& ret) {
  const double x = arg[0];
  const double y = arg[1];
  const double v = Parameters().getParam("viscosity", 1.0);
  const double F = 1; // std::exp( -8 * std::pow( M_PI, 2 ) * time );
  const double P = 2 * M_PI; // pi_factor;
  const double S_x = std::sin(P * x);
  const double S_y = std::sin(P * y);
  const double S_2x = std::sin(2 * P * x);
  const double S_2y = std::sin(2 * P * y);
  const double C_x = std::cos(P * x);
  const double C_y = std::cos(P * y);

  ret[0] = (-1 / v) * C_x * S_y * F;
  ret[1] = (1 / v) * S_x * C_y * F;
}

/**
*  \brief  describes the dirichlet boundary data
*
*  \tparam DirichletTraitsImp
*          types like functionspace, range type, etc
*
*  \todo   extensive docu with latex
**/
template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

  /**
  * \brief  evaluates the dirichlet data
  * \param  arg
  *         point to evaluate at
  * \param  ret
  *         value of dirichlet boundary data at given point
  **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const { assert(false); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double lambda_;
};
template <class FunctionSpaceImp>
class Beta : public Fem::Function<FunctionSpaceImp, Beta<FunctionSpaceImp>> {
public:
  typedef Beta<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  Beta(double dummy, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
       const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~Beta() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

  /**
  * \brief  evaluates the dirichlet data
  * \param  arg
  *         point to evaluate at
  * \param  ret
  *         value of dirichlet boundary data at given point
  **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const { VelocityEvaluate(lambda_, 0, arg, ret); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double lambda_;
};
template <class FunctionSpaceImp>
class ExactConvection : public Fem::Function<FunctionSpaceImp, ExactConvection<FunctionSpaceImp>> {
public:
  typedef ExactConvection<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  ExactConvection(double dummy, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                  const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~ExactConvection() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

  /**
  * \brief  evaluates the dirichlet data
  * \param  arg
  *         point to evaluate at
  * \param  ret
  *         value of dirichlet boundary data at given point
  **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const { VelocityEvaluate(lambda_, 0, arg, ret); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double lambda_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class Velocity : public TimeFunction<FunctionSpaceImp, Velocity<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Velocity<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  Velocity(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~Velocity() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

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

template <class FunctionSpaceImp, class TimeProviderImp>
class Pressure : public TimeFunction<FunctionSpaceImp, Pressure<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Pressure<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *
   *  doing nothing besides Base init
   **/
  Pressure(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0))
    , shift_(0.0) {}

  /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
  ~Pressure() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    Dune::CompileTimeChecker<(dim_ == 2)> Pressure_Unsuitable_WorldDim;
    const double x = arg[0];
    const double y = arg[1];
    const double v = Parameters().getParam("viscosity", 1.0);
    const double F = std::exp(-16 * std::pow(M_PI, 2) * time);
    const double C1 = std::cos(4 * M_PI * (x + 0.25));
    const double C2 = std::cos(4 * M_PI * (y + 0.5));

    ret = (-1 / (4 * v)) * (C1 + C2) * F;
  }

  template <class DiscreteFunctionSpace>
  void setShift(const DiscreteFunctionSpace& space) {
    //					shift_ = -1 * Stuff::meanValue( *this, space );
  }

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
  double shift_;
};
} // end namespace TestCaseTaylor2D

namespace TrivialTestCase {
template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
   **/
  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  ~Force() {}

  inline void evaluate(const DomainType& arg, RangeType& ret) const {
    const double x = arg[0];
    const double y = arg[1];
    const double v = viscosity_;
    const double a = Parameters().getParam("alpha", 1.0);
    const double P = 2 * M_PI; // pi_factor;
    const double E = 1; // std::exp( -2 * std::pow( P, 2 ) * viscosity_ * time );
    const double F = 1; // std::exp( -4 * std::pow( P, 2 ) * viscosity_ * time );
    const double S_x = std::sin(P * x);
    const double S_y = std::sin(P * y);
    const double S_2x = std::sin(2 * P * x);
    const double S_2y = std::sin(2 * P * y);
    const double C_x = std::cos(P * x);
    const double C_y = std::cos(P * y);

    // alpha * u part
    RangeType u;
    VelocityEvaluate(0, 0, arg, u);
    RangeType u_a = u;
    u_a *= a;
    ret = u_a;

    // laplace
    //					  ret[0] -= viscosity_ * 2;
    //					  ret[1] += viscosity_ * 2;
    // beta = u
    //					  ret[0] = x * ( 1 + a );
    //					  ret[1] = y * ( 1 - a );
    // beta = 0
    //					  ret[0] = x * ( a );
    //					  ret[1] = y * ( - a );

    //					  ret[0] = ( a + v * P * P ) * S_y * C_y ;
    ////					  ret[0] = a ;
    //					  ret[1] = 0;

    // beta = (1,0)
    //					  ret[0] += u[0];
    //					  ret[1] += -u[1];
    ret[0] += Parameters().getParam("conv", 1.0);
  }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double lambda, const double time, const DomainType& arg, RangeType& ret) {
  const double x = arg[0];
  const double y = arg[1];
  const double v = Parameters().getParam("viscosity", 1.0);
  const double F = std::exp(-8 * std::pow(M_PI, 2) * time);
  const double P = 2 * M_PI;
  const double S_x = std::sin(P * x);
  const double S_y = std::sin(P * y);
  const double S_2x = std::sin(2 * P * x);
  const double S_2y = std::sin(2 * P * y);
  const double C_x = std::cos(P * x);
  const double C_y = std::cos(P * y);

  //			ret[0] = 4;
  //			ret[1] = 0;
  ret[0] = x;
  ret[1] = -y;
}
template <class FunctionSpaceImp>
class Beta : public Fem::Function<FunctionSpaceImp, Beta<FunctionSpaceImp>> {
public:
  typedef Beta<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
   **/
  Beta(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  ~Beta() {}

  inline void evaluate(const DomainType& arg, RangeType& ret) const {
    VelocityEvaluate(0, 0, arg, ret);
    ret = RangeType(0);
    ret[0] = Parameters().getParam("conv", 1.0);
    //					  ret[1] = Parameters().getParam( "conv", 1.0 );
  }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class FunctionSpaceImp>
class ExactConvection : public Fem::Function<FunctionSpaceImp, ExactConvection<FunctionSpaceImp>> {
public:
  typedef ExactConvection<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
   **/
  ExactConvection(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  ~ExactConvection() {}

  inline void evaluate(const DomainType& arg, RangeType& ret) const {
    ret = RangeType(0);
    ret[0] = Parameters().getParam("conv", 1.0);
    //					 ret[1] = -arg[1];
  }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

/**
*  \brief  describes the dirichlet boundary data
*
*  \tparam DirichletTraitsImp
*          types like functionspace, range type, etc
*
*  \todo   extensive docu with latex
**/
template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

  /**
  * \brief  evaluates the dirichlet data
  * \param  arg
  *         point to evaluate at
  * \param  ret
  *         value of dirichlet boundary data at given point
  **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const { assert(false); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double lambda_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class Velocity : public TimeFunction<FunctionSpaceImp, Velocity<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Velocity<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  Velocity(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0)) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~Velocity() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(lambda_, time, arg, ret);
  }

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
template <class FunctionSpaceImp, class TimeProviderImp>
class Pressure : public TimeFunction<FunctionSpaceImp, Pressure<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Pressure<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *
   *  doing nothing besides Base init
   **/
  Pressure(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , lambda_(Parameters().getParam("lambda", 0.0))
    , shift_(0.0) {}

  /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
  ~Pressure() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    Dune::CompileTimeChecker<(dim_ == 2)> Pressure_Unsuitable_WorldDim;
    const double x = arg[0];
    const double y = arg[1];
    const double v = Parameters().getParam("viscosity", 1.0);
    const double F = std::exp(-16 * std::pow(M_PI, 2) * time);
    const double C1 = std::cos(4 * M_PI * (x + 0.25));
    const double C2 = std::cos(4 * M_PI * (y + 0.5));

    ret = 0;
  }

  template <class DiscreteFunctionSpace>
  void setShift(const DiscreteFunctionSpace& space) {
    //					shift_ = -1 * Stuff::meanValue( *this, space );
  }

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
  double shift_;
};

} // end namespace TrivialTestCase

namespace AdapterFunctionsVectorial {

static const double pi_factor = M_PI / 2.; // controls number of vortices
struct Evals {
  template <class DomainType>
  Evals(const DomainType& arg, const double time)
    : x(arg[0])
    , y(arg[1])
    ,
    //			time_(time),
    time_(0)
    , v(Parameters().getParam("viscosity", 1.0))
    , alpha(Parameters().getParam("alpha", 1.0))
    , P(pi_factor)
    , E(std::exp(-2. * std::pow(P, 2.) * v * time_))
    , F(std::exp(-4. * std::pow(P, 2.) * v * time_))
    , S_x(std::sin(P * x))
    , S_y(std::sin(P * y))
    , S_2x(std::sin(2. * P * x))
    , S_2y(std::sin(2. * P * y))
    , C_2x(std::cos(2. * P * x))
    , C_2y(std::cos(2. * P * y))
    , C_x(std::cos(P * x))
    , C_y(std::cos(P * y)) {}
  double x;
  double y;
  double time_;
  double v;
  double alpha;
  double P;
  double E;
  double F;
  double S_x;
  double S_y;
  double C_2x;
  double C_2y;
  double S_2x;
  double S_2y;
  double C_x;
  double C_y;
};

template <class FunctionSpaceImp>
class Force : public Dune::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Dune::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
   **/
  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  /**
   *  \brief  destructor
   *  doing nothing
   **/
  ~Force() {}

  /**
   *  \brief  evaluates the force
   *  \param  arg
   *          point to evaluate at
   *  \param  ret
   *          value of force at given point
   **/
  inline void evaluate(const double time, const DomainType& arg, RangeType& ret) const {
    Evals evals(arg, time);
    ret = RangeType(0);
    // diff
    RangeType laplace;
    VelocityLaplaceEvaluateTime(time, arg, laplace);
    laplace *= evals.v;
    ret -= laplace;

    // druck
    //					  ret[0] += 0.5 * evals.P * evals.F * evals.S_2x;
    //					  ret[1] += 0.5 * evals.P * evals.F * evals.S_2y;

    // conv
    RangeType conv;
    VelocityConvectionEvaluateTime(time, arg, conv);
    ret += conv;

    // zeitableitung
    RangeType u;
    VelocityEvaluate(0, 0, time, arg, u);
    //				  ret[0] += ( -2 * std::pow( evals.P, 2 ) * evals.v ) * u[0];
    //				  ret[1] += ( -2 * std::pow( evals.P, 2 ) * evals.v ) * u[1];

    u *= evals.alpha;
    ret += u;

    //					  ret *=  Parameters().getParam( "rhs_factor", 1.0 );

    ret[0] = -2;
    ret[1] = -2;
    ret[0] += 2 * std::pow(arg[0], 2.0) * arg[1];
    ret[1] += 2 * std::pow(arg[1], 2.0) * arg[0];
  }
  inline void evaluate(const DomainType& arg, RangeType& ret) const { evaluate(0, arg, ret); }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double /*parameter_a*/, const double /*parameter_d*/, const double time,
                      const DomainType& arg, RangeType& ret) {
  //			Evals evals( arg, time );
  //			ret[0] = -1 *	evals.C_x * evals.S_y * evals.E;
  //			ret[1] =		evals.S_x * evals.C_y * evals.E;
  ret[0] = arg[1] * arg[1];
  ret[1] = arg[0] * arg[0];
}

/**
*  \brief  describes the dirichlet boundary data
*
*  \tparam DirichletTraitsImp
*          types like functionspace, range type, etc
*
*  \todo   extensive docu with latex
**/
template <class FunctionSpaceImp>
class DirichletData : public Dune::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Dune::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
  }

  /**
  * \brief  evaluates the dirichlet data
  * \param  arg
  *         point to evaluate at
  * \param  ret
  *         value of dirichlet boundary data at given point
  **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const {
    VelocityEvaluate(parameter_a_, parameter_d_, 0, arg, ret);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

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
  Velocity(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~Velocity() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
  }

  /**
 * \brief  evaluates the dirichlet data
 * \param  arg
 *         point to evaluate at
 * \param  ret
 *         value of dirichlet boundary data at given point
 **/
  inline void evaluate(const DomainType& arg, RangeType& ret) const {
    VelocityEvaluate(parameter_a_, parameter_d_, 0, arg, ret);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class Pressure
    : public Dune::TimeFunction<FunctionSpaceImp, Pressure<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Pressure<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *
   *  doing nothing besides Base init
   **/
  Pressure(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
  ~Pressure() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    Dune::CompileTimeChecker<(dim_ == 2)> Pressure_Unsuitable_WorldDim;
    Evals evals(arg, time);
    ret[0] = -0.25 * (evals.C_2x + evals.C_2y) * evals.F;
  }

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
  const double parameter_a_;
  const double parameter_d_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class PressureGradient : public Dune::TimeFunction<
    FunctionSpaceImp, PressureGradient<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef PressureGradient<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *
   *  doing nothing besides Base init
   **/
  PressureGradient(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                   const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
  ~PressureGradient() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    Evals evals(arg, time);
    ret[0] = 0.5 * evals.P * evals.F * evals.S_2x;
    ret[1] = 0.5 * evals.P * evals.F * evals.S_2y;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

template <class DomainType, class RangeType>
void VelocityLaplaceEvaluateTime(const double time, const DomainType& arg, RangeType& ret) {
  Evals evals(arg, time);
  ret[0] = 2 * evals.C_x * evals.E * evals.P * evals.S_y * evals.P;
  ret[1] = -2 * evals.C_y * evals.E * evals.P * evals.S_x * evals.P;
}
template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityLaplace
    : public Dune::TimeFunction<FunctionSpaceImp, VelocityLaplace<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityLaplace<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  VelocityLaplace(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                  const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~VelocityLaplace() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    //				dune_static_assert( dim_ == 2 , "DirichletData_Unsuitable_WorldDim");
    VelocityLaplaceEvaluateTime(time, arg, ret);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
template <class DomainType, class RangeType>
void VelocityConvectionEvaluateTime(const double time, const DomainType& arg, RangeType& ret) {
  Evals evals(arg, time);
  RangeType u;
  VelocityEvaluate(0, 0, time, arg, u);
  ret[0] = u[0] * evals.P * evals.S_x * evals.S_y;
  ret[0] -= u[1] * evals.P * evals.C_x * evals.C_y;

  ret[1] = u[0] * evals.P * evals.C_x * evals.C_y;
  ret[1] -= u[1] * evals.P * evals.S_x * evals.S_y;

  //		ret[0] =  E * E *P * C_x * S_x * ( C_y * C_y - S_y * S_y );//eigentlich falsch
  //		ret[1] = - E * E *P * S_y * C_y * ( S_x * S_x - C_x * C_x );
}
template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityConvection : public Dune::TimeFunction<
    FunctionSpaceImp, VelocityConvection<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityConvection<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  VelocityConvection(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                     const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~VelocityConvection() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    VelocityConvectionEvaluateTime(time, arg, ret);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityGradientX : public Dune::TimeFunction<
    FunctionSpaceImp, VelocityGradientX<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityGradientX<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  VelocityGradientX(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                    const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~VelocityGradientX() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    Evals evals(arg, time);
    ret[0] = evals.P * evals.S_x * evals.S_y * evals.E;
    ret[1] = -evals.P * evals.C_x * evals.C_y * evals.E;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityGradientY : public Dune::TimeFunction<
    FunctionSpaceImp, VelocityGradientY<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityGradientY<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  VelocityGradientY(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                    const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~VelocityGradientY() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    Evals evals(arg, time);
    ret[0] = evals.P * evals.C_x * evals.C_y * evals.E;
    ret[1] = -evals.P * evals.S_x * evals.S_y * evals.E;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityGradient : public Dune::TimeFunction<
    FunctionSpaceImp, VelocityGradient<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityGradient<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
  *  \brief  constructor
  *
  *  doing nothing besides Base init
  **/
  VelocityGradient(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                   const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  /**
  *  \brief  destructor
  *
  *  doing nothing
  **/
  ~VelocityGradient() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    Evals evals(arg, time);
    // grad_x
    ret(0, 0) = evals.P * evals.S_x * evals.S_y * evals.E;
    ret(0, 1) = -evals.P * evals.C_x * evals.C_y * evals.E;
    // grad_y
    ret(1, 0) = evals.P * evals.C_x * evals.C_y * evals.E;
    ret(1, 1) = -evals.P * evals.S_x * evals.S_y * evals.E;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
} // end namespace AdapterFunctionsVectorial
namespace TimeDisc {

NULLFUNCTION_TP(VelocityGradientX)
NULLFUNCTION_TP(VelocityGradientY)
NULLFUNCTION_TP(VelocityGradient)

template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  /**
   *  \brief  constructor
   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
   **/
  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  /**
   *  \brief  destructor
   *  doing nothing
   **/
  ~Force() {}

  /**
   *  \brief  evaluates the force
   *  \param  arg
   *          point to evaluate at
   *  \param  ret
   *          value of force at given point
   **/
  inline void evaluate(const double time, const DomainType& arg, RangeType& ret) const {
    const double x = arg[0];
    const double y = arg[1];
    const double v = viscosity_;
    const double alpha = Parameters().getParam("alpha", 1.0);
    ret[0] = std::pow(time, 3.0) * arg[1] * arg[1];
    ret[1] = std::pow(time, 2.0) * arg[0];
    ret *= alpha;
    // laplce
    ret[0] += -2 * std::pow(time, 3.0) * v;
    ret[1] += 0;
    // grad p
    //					  ret[0] += time;
    //					  ret[1] += 1;
    // conv
    ret[0] += std::pow(time, 5.0) * 2 * x * y;
    ret[1] += std::pow(time, 5.0) * y * y;
    // dt u
    //					  ret[0] += std::pow(time,2.0)*3*y*y;
    //					  ret[1] += 2*time*x;

    //					  ret *=0;
  }
  inline void evaluate(const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret, const IntersectionType& intersection) const {
    ret[0] = std::pow(time, 3.0) * arg[1] * arg[1];
    ret[1] = std::pow(time, 2.0) * arg[0];
  }

  inline void evaluate(const DomainType& arg, RangeType& ret) const { assert(false); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class Velocity : public TimeFunction<FunctionSpaceImp, Velocity<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Velocity<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  Velocity(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~Velocity() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(dim_ == 2, "Wrong world dim");
    ret[0] = std::pow(time, 3.0) * arg[1] * arg[1];
    ret[1] = std::pow(time, 2.0) * arg[0];
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class PressureGradient
    : public TimeFunction<FunctionSpaceImp, PressureGradient<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef PressureGradient<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  PressureGradient(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                   const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~PressureGradient() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    ret[0] = time;
    ret[1] = 1;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class Pressure : public TimeFunction<FunctionSpaceImp, Pressure<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef Pressure<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  Pressure(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
           const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~Pressure() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    ret = time * arg[0] + arg[1] - ((time + 1) / 2.0);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityLaplace
    : public Dune::TimeFunction<FunctionSpaceImp, VelocityLaplace<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityLaplace<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  VelocityLaplace(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                  const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~VelocityLaplace() {}

  void evaluateTime(const double time, const DomainType& /*arg*/, RangeType& ret) const {
    ret[0] = 2 * std::pow(time, 3.0);
    ret[1] = 0;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityConvection : public Dune::TimeFunction<
    FunctionSpaceImp, VelocityConvection<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityConvection<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef Dune::TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  VelocityConvection(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                     const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~VelocityConvection() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    const double x = arg[0];
    const double y = arg[1];
    ret[0] = std::pow(time, 5.0) * 2 * x * y;
    ret[1] = std::pow(time, 5.0) * y * y;
    ;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
} // end namespace TimeDisc
}
}

#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/oseen/modelinterface.hh>
#include <dune/oseen/pass.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/navier/fractionaldatawriter.hh>
#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/navier/thetascheme_traits.hh>
#include <cmath>

namespace Dune {
namespace ConvDiff {

template <class TimeProviderType, class GridPartImp, template <class> class ForceFuntionType,
          template <class> class AnalyticalDirichletDataImp, int gridDim, int sigmaOrder,
          int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder,
          template <class, class, int, template <class> class BaseFunctionStorageImp = Dune::CachingStorage>
          class GalerkinSpaceImp = Dune::DiscontinuousGalerkinSpace>
class DiscreteModelTraits {
public:

  //! for CRTP trick
  typedef Dune::DiscreteOseenModelDefault<DiscreteModelTraits> DiscreteModelType;

  //! we use caching quadratures for the entities
  typedef Dune::CachingQuadrature<GridPartImp, 0> VolumeQuadratureType;

  //! we use caching quadratures for the faces
  typedef Dune::CachingQuadrature<GridPartImp, 1> FaceQuadratureType;

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
  typedef GalerkinSpaceImp<VelocityFunctionSpaceType, GridPartImp, velocitySpaceOrder>
  DiscreteVelocityFunctionSpaceType;

  //! function space type for the pressure
  typedef Dune::FunctionSpace<double, double, gridDim, 1> PressureFunctionSpaceType;

  //! discrete function space type for the pressure
  typedef GalerkinSpaceImp<PressureFunctionSpaceType, GridPartImp, pressureSpaceOrder>
  DiscretePressureFunctionSpaceType;

public:

  //! discrete function space wrapper type
  typedef Dune::DiscreteOseenFunctionSpaceWrapper<Dune::DiscreteOseenFunctionSpaceWrapperTraits<
      DiscreteVelocityFunctionSpaceType, DiscretePressureFunctionSpaceType>> DiscreteOseenFunctionSpaceWrapperType;

private:

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

public:

  //! function space type for sigma
  typedef Dune::MatrixFunctionSpace<double, double, gridDim, gridDim, gridDim> SigmaFunctionSpaceType;

  //! discrete function space type for sigma
  typedef GalerkinSpaceImp<SigmaFunctionSpaceType, GridPartImp, sigmaSpaceOrder> DiscreteSigmaFunctionSpaceType;

public:

  //! discrete function type for sigma
  typedef Dune::AdaptiveDiscreteFunction<DiscreteSigmaFunctionSpaceType> DiscreteSigmaFunctionType;

  //! function type for the analytical force
  typedef ForceFuntionType<VelocityFunctionSpaceType> AnalyticalForceFunctionType;

  typedef AnalyticalForceFunctionType AnalyticalForceType;

  //! function type for the analytical dirichlet data
  typedef typename Dune::NavierStokes::StokesStep::DirichletAdapterFunctionTraits<
      AnalyticalDirichletDataImp, TimeProviderType>::template Implementation<VelocityFunctionSpaceType, GridPartImp>
  AnalyticalDirichletDataTraitsImplementation;
  typedef typename AnalyticalDirichletDataTraitsImplementation::AnalyticalDirichletDataType AnalyticalDirichletDataType;

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

template <class CommunicatorImp, class GridPartImp, int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder,
          int pressureOrder = sigmaOrder>
struct Traits {
  typedef Traits<CommunicatorImp, GridPartImp, gridDim, sigmaOrder, velocityOrder, pressureOrder> ThisType;
  typedef GridPartImp GridPartType;
  typedef Dune::NavierStokes::ThetaSchemeDescription<1> SchemeDescriptionType;
  typedef Dune::NavierStokes::FractionalTimeProvider<SchemeDescriptionType, CommunicatorImp> TimeProviderType;

  typedef DiscreteModelTraits<TimeProviderType, GridPartType, CONVDIFF_DATA_NAMESPACE::Force,
                              CONVDIFF_DATA_NAMESPACE::DirichletData, gridDim, sigmaOrder, velocityOrder,
                              pressureOrder> OseenModelTraits;
  typedef Dune::DiscreteOseenModelDefault<OseenModelTraits> OseenModelType;
  typedef Dune::StartPass<typename OseenModelTraits::DiscreteOseenFunctionWrapperType, -1> StartPassType;
  typedef Dune::OseenPass<OseenModelType, StartPassType> OseenPassType;

  typedef CONVDIFF_DATA_NAMESPACE::Pressure<typename OseenModelTraits::PressureFunctionSpaceType, TimeProviderType>
  ExactPressureType;
  typedef CONVDIFF_DATA_NAMESPACE::Velocity<typename OseenModelTraits::VelocityFunctionSpaceType, TimeProviderType>
  ExactVelocityType;
  typedef CONVDIFF_DATA_NAMESPACE::Velocity<typename OseenModelTraits::VelocityFunctionSpaceType, TimeProviderType>
  ConvectionType;
  typedef CONVDIFF_DATA_NAMESPACE::VelocityConvection<typename OseenModelTraits::VelocityFunctionSpaceType,
                                                      TimeProviderType> ExactConvectionType;
  typedef CONVDIFF_DATA_NAMESPACE::VelocityGradientX<typename OseenModelTraits::VelocityFunctionSpaceType,
                                                     TimeProviderType> VelocityGradientXType;
  typedef CONVDIFF_DATA_NAMESPACE::VelocityGradientY<typename OseenModelTraits::VelocityFunctionSpaceType,
                                                     TimeProviderType> VelocityGradientYType;
  typedef CONVDIFF_DATA_NAMESPACE::VelocityGradient<typename OseenModelTraits::SigmaFunctionSpaceType, TimeProviderType>
  VelocityGradientType;
  typedef CONVDIFF_DATA_NAMESPACE::VelocityLaplace<typename OseenModelTraits::VelocityFunctionSpaceType,
                                                   TimeProviderType> VelocityLaplaceType;

  typedef Dune::NavierStokes::ExactSolution<ThisType> ExactSolutionType;
  typedef typename OseenModelType::DiscreteOseenFunctionWrapperType DiscreteOseenFunctionWrapperType;
  typedef typename OseenModelType::DiscreteOseenFunctionSpaceWrapperType DiscreteOseenFunctionSpaceWrapperType;
  typedef typename DiscreteOseenFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType::FunctionSpaceType
  PressureFunctionSpaceType;
  typedef typename DiscreteOseenFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType::FunctionSpaceType
  VelocityFunctionSpaceType;
};

} // end namespace ConvDiff
} // end namespace DUne
#endif // CONV_DIFF_HH

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
