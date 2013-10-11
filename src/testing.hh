#ifndef NAVIER_TESTING_HH
#define NAVIER_TESTING_HH

#include <dune/stuff/misc.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/common/tuples.hh>

#include <boost/format.hpp>

namespace Testing {

namespace AdapterFunctions {

template <class FunctionSpaceImp>
class Force : public Dune::Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Dune::Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
  inline void evaluate(const double /*time*/, const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }
  inline void evaluate(const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }

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
class DirichletData : public Dune::Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Dune::Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];

    ret[0] = x * x + y;
    ret[1] = y * y + x;
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
    //				dune_static_assert( ( dim_ == 2 ) > Pressure_Unsuitable_WorldDim;
    const double x = arg[0];
    const double y = arg[1];

    ret[0] = std::sin(x);
  }

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
    //				dune_static_assert( ( dim_ == 2 ) > Pressure_Unsuitable_WorldDim;
    const double x = arg[0];
    const double y = arg[1];

    ret[0] = std::cos(x);
    ret[1] = 0;
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];

    ret[0] = 2;
    ret[1] = 2;
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
    dune_static_assert((dim_ == 2), "VelocityConvection_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];
    const double u_1 = x * x + y;
    const double u_2 = y * y + x;
    ret[0] = u_1 * 2 * x + u_2;
    ret[1] = u_1 + u_2 * 2 * y;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
} // end namespace TestCase2D

namespace AdapterFunctionsScalar {

template <class FunctionSpaceImp>
class Force : public Dune::Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Dune::Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    const double viscosity = Parameters().getParam("viscosity", 1.0);
    const double x = arg[0];
    const double y = arg[1];
    // conv
    ret[0] = 2 * x * y;
    ret[1] = y * y;
    ret *= std::pow(time, 5.0f);

    // laplace
    ret[0] -= viscosity * 2 * std::pow(time, 3.0f);
    //				  ret[1] -= 0;

    // press grad
    ret[0] += time;
    ret[1] += 1;

    // time diff
    ret[0] += std::pow(arg[0], 3.0) * 2 * time * time;
    ret[1] += time * time;
  }
  inline void evaluate(const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};
template <class DomainType, class RangeType>
void VelocityEvaluate(const double /*parameter_a*/, const double /*parameter_d*/, const double time,
                      const DomainType& arg, RangeType& ret) {
  ret[0] = std::pow(time, 3.0f) * arg[1];
  ret[1] = std::pow(time, 2.0f) * arg[0];
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
class DirichletData : public Dune::Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Dune::Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert((dim_ == 2), "Pressure_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];

    ret[0] = (time * arg[0] + arg[1]) - ((time + 1) / 2.0f);
  }

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
    //				dune_static_assert( ( dim_ == 2 ) > Pressure_Unsuitable_WorldDim;
    //				const double x			= arg[0];
    //				const double y			= arg[1];

    ret[0] = time;
    ret[1] = 1;
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    ret[0] = 2 * std::pow(time, 3.0f);
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];
    ret[0] = 2 * x * y;
    ret[1] = y * y;
    ret *= std::pow(time, 5.0f);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
} // end namespace AdapterFunctionsScalar

namespace AdapterFunctionsVectorial {

static const double pi_factor = M_PI; // controls number of vortices
struct Evals {
  template <class DomainType>
  Evals(const DomainType& arg, const double time)
    : x(arg[0])
    , y(arg[1])
    , time_(time)
    ,
    //			time_(0),
    v(Parameters().getParam("viscosity", 1.0))
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
class Force : public Dune::Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Dune::Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    ret[0] += 0.5 * evals.P * evals.F * evals.S_2x;
    ret[1] += 0.5 * evals.P * evals.F * evals.S_2y;

    // conv
    RangeType conv;
    VelocityConvectionEvaluateTime(time, arg, conv);
    ret += conv;

    // zeitableitung
    RangeType u;
    VelocityEvaluate(0, 0, time, arg, u);
    ret[0] += (-2 * std::pow(evals.P, 2) * evals.v) * u[0];
    ret[1] += (-2 * std::pow(evals.P, 2) * evals.v) * u[1];

    ret *= Parameters().getParam("rhs_factor", 1.0);
  }
  inline void evaluate(const DomainType& /*arg*/, RangeType& ret) const { assert(false); }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double /*parameter_a*/, const double /*parameter_d*/, const double time,
                      const DomainType& arg, RangeType& ret) {
  Evals evals(arg, time);
  ret[0] = -1 * evals.C_x * evals.S_y * evals.E;
  ret[1] = evals.S_x * evals.C_y * evals.E;
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
class DirichletData : public Dune::Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Dune::Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert((dim_ == 2), "Pressure_Unsuitable_WorldDim");
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
    //				dune_static_assert( ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
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
  ret[0] = -evals.E * evals.E * evals.P * evals.C_x * evals.S_x; // eigentlich richtig
  ret[1] = -evals.E * evals.E * evals.P * evals.S_y * evals.C_y;

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
    ret = RangeType(0);
    Evals evals(arg, time);
    // velo eval
    //				ret[0] = -1 *	evals.C_x * evals.S_y * evals.E;
    //				ret[1] =		evals.S_x * evals.C_y * evals.E;
    ret(0, 0) = evals.P * evals.S_x * evals.S_y * evals.E;
    ret(0, 1) = -evals.P * evals.C_x * evals.C_y * evals.E;
    ret(1, 0) = evals.P * evals.C_x * evals.C_y * evals.E;
    ret(1, 1) = -evals.P * evals.S_x * evals.S_y * evals.E;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
} // end namespace AdapterFunctionsVectorial

namespace AdapterFunctionsVisco {

static const double pi_factor = M_PI; // controls number of vortices
struct Evals {
  template <class DomainType>
  Evals(const DomainType& arg, const double time)
    : x(arg[0])
    , y(arg[1])
    , time_(time)
    ,
    //			time_(1),
    v(Parameters().getParam("viscosity", 1.0))
    , P(M_PI)
    , E(std::exp(-8 * std::pow(M_PI, 2) * time_))
    , F(std::exp(-16 * std::pow(M_PI, 2) * time_))
    , S_x(std::sin(2 * M_PI * (x + 0.25)))
    , S_y(std::sin(2 * M_PI * (y + 0.5)))
    , S_2x(std::sin(4 * M_PI * (x + 0.25)))
    , S_2y(std::sin(4 * M_PI * (x + 0.5)))
    , C_2x(std::cos(4 * M_PI * (x + 0.25)))
    , C_2y(std::cos(4 * M_PI * (x + 0.5)))
    , C_x(std::cos(2 * M_PI * (y + 0.25)))
    , C_y(std::cos(2 * M_PI * (y + 0.5))) {}
  const double x;
  const double y;
  const double time_;
  const double v;
  const double P;
  const double E;
  const double F;
  const double S_x;
  const double S_y;
  const double C_2x;
  const double C_2y;
  const double S_2x;
  const double S_2y;
  const double C_x;
  const double C_y;
};

template <class DomainType, class RangeType>
void PressureGradientEvaluateTime(const double time, const DomainType& arg, RangeType& ret) {
  // ret[0] = 1 - std::pow( evals.x + evals.y);
  Evals evals(arg, time);
  ret = RangeType(0);
  //		ret[0] += evals.y;
  //		ret[1] += evals.x;
}

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
    PressureGradientEvaluateTime(time, arg, ret);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

template <class FunctionSpaceImp>
class Force : public Dune::Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Dune::Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    //				  RangeType laplace;
    //				  VelocityLaplaceEvaluateTime( time, arg, laplace );
    //				  laplace *= evals.v;
    //				  ret -= laplace;

    //				  //u
    //				  RangeType u;
    //				  VelocityEvaluate( 0,0,time, arg, u);
    //				  ret += u;

    // druck
    RangeType pressure_gradient;
    PressureGradientEvaluateTime(time, arg, pressure_gradient);
    ret += pressure_gradient;

    // conv
    RangeType conv;
    VelocityConvectionEvaluateTime(time, arg, conv);
    ret += conv;

    // zeitableitung
    //				  RangeType u;
    //				  VelocityEvaluate( 0, 0, time, arg, u);
    //				  ret[0] += ( -8 * std::pow( M_PI, 2 ) ) * u[0];
    //				  ret[1] += ( -8 * std::pow( M_PI, 2 ) ) * u[1];

    //				  ret *=  Parameters().getParam( "rhs_factor", 1.0 );
  }
  inline void evaluate(const DomainType& /*arg*/, RangeType& ret) const { assert(false); }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double /*parameter_a*/, const double /*parameter_d*/, const double time,
                      const DomainType& arg, RangeType& ret) {
  Evals evals(arg, time);
  ret[0] = evals.x;
  ret[1] = -evals.y;
  ret[0] = 1;
  ret[1] = 0;
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
class DirichletData : public Dune::Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Dune::Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert((dim_ == 2), "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert((dim_ == 2), "Pressure_Unsuitable_WorldDim");
    Evals evals(arg, time);
    ret[0] = evals.x * evals.y;
    ret[0] = 0;
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

template <class DomainType, class RangeType>
void VelocityLaplaceEvaluateTime(const double time, const DomainType& arg, RangeType& ret) {
  Evals evals(arg, time);
  ret[0] = 0;
  ret[1] = 0;
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
    //				dune_static_assert( ( dim_ == 2 ), "DirichletData_Unsuitable_WorldDim");
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
  // beta = (1,0)
  ret[0] = 0;
  ret[1] = 0;
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

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const { ret = RangeType(0); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
} // end namespace AdapterFunctionsVisco

} // end namespace AdapterFunctions

template <class T1, class T2, class T3>
struct TupleSerializer {
  typedef Dune::tuple<
      const typename T1::DiscreteVelocityFunctionType*, const typename T1::DiscretePressureFunctionType*,
      const typename T2::DiscreteVelocityFunctionType*, const typename T2::DiscretePressureFunctionType*,
      const typename T3::DiscreteVelocityFunctionType*, const typename T3::DiscretePressureFunctionType*> TupleType;

  static TupleType& getTuple(T1& t1, T2& t2, T3& t3) {
    static TupleType t(&(t1.discreteVelocity()), &(t1.discretePressure()), &(t2.discreteVelocity()),
                       &(t2.discretePressure()), &(t3.discreteVelocity()), &(t3.discretePressure()));
    return t;
  }
};

#endif // TESTING_HH

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
