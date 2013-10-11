#ifndef TESTDATA_HH
#define TESTDATA_HH

#include <dune/stuff/timefunction.hh>
#include <dune/stuff/functions.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/fem/misc/validator.hh>

namespace Dune {
namespace NavierStokes {

namespace TestCase3D {
template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  ~Force() {}

  inline void evaluate(const double /*time*/, const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }
  inline void evaluate(const DomainType& arg, RangeType& ret) const { assert(false); }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double parameter_a, const double parameter_d, const double time, const DomainType& arg,
                      RangeType& ret) {
  const double x = arg[0];
  const double y = arg[1];
  const double z = arg[2];
  const double e_d_d_t = std::exp(-1 * parameter_d * parameter_d * time);
#define PM +
  ret[0] = -parameter_a * ((std::exp(parameter_a * x) * std::sin(parameter_a * y PM /*+-*/ parameter_d * z)) +
                           (std::exp(parameter_a * z) * std::cos(parameter_a * x PM /*+-*/ parameter_d * y))) *
           e_d_d_t;
  ret[1] = -parameter_a * ((std::exp(parameter_a * y) * std::sin(parameter_a * z PM /*+-*/ parameter_d * x)) +
                           (std::exp(parameter_a * x) * std::cos(parameter_a * y PM /*+-*/ parameter_d * z))) *
           e_d_d_t;
  ret[2] = -parameter_a * ((std::exp(parameter_a * z) * std::sin(parameter_a * x PM /*+-*/ parameter_d * y)) +
                           (std::exp(parameter_a * y) * std::cos(parameter_a * z PM /*+-*/ parameter_d * x))) *
           e_d_d_t;
#undef PM
}

template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 3, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert(dim_ == 3, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert(dim_ == 3, "Pressure_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];
    const double z = arg[2];
    const double e_d_d_t = std::exp(-2 * parameter_d_ * parameter_d_ * time);
#define PM +
    ret = (-parameter_a_ * parameter_a_ / 2.0) *
          (std::exp(2.0 * parameter_a_ * x) + std::exp(2.0 * parameter_a_ * y) + std::exp(2.0 * parameter_a_ * z) +
           (2 * std::sin(parameter_a_ * x PM parameter_d_ * y) * std::cos(parameter_a_ * z PM parameter_d_ * x) *
            std::exp(parameter_a_ * (y + z))) +
           (2 * std::sin(parameter_a_ * y PM parameter_d_ * z) * std::cos(parameter_a_ * x PM parameter_d_ * y) *
            std::exp(parameter_a_ * (z + x))) +
           (2 * std::sin(parameter_a_ * z PM parameter_d_ * x) * std::cos(parameter_a_ * y PM parameter_d_ * z) *
            std::exp(parameter_a_ * (x + y)))) *
          e_d_d_t;
#undef PM
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

  void evaluateTime(const double t, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(dim_ == 3, "Pressure_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];
    const double z = arg[2];
    const double a = parameter_a_;
    const double b = parameter_d_;

    const double sxy = std::sin(a * x + b * y);
    const double sxz = std::sin(a * x + b * z);
    const double syz = std::sin(a * y + b * z);
    const double syx = std::sin(a * y + b * x);
    const double szy = std::sin(a * z + b * y);
    const double szx = std::sin(a * z + b * x);

    const double cxy = std::cos(a * x + b * y);
    const double cxz = std::cos(a * x + b * z);
    const double cyz = std::cos(a * y + b * z);
    const double cyx = std::cos(a * y + b * x);
    const double czy = std::cos(a * z + b * y);
    const double czx = std::cos(a * z + b * x);

    const double exy = std::exp(x + y);
    const double exz = std::exp(x + z);
    const double eyz = std::exp(z + y);
    const double fr = -0.5 * a * a * std::exp(-2 * b * b * t);

    ret[0] = -2 * a * exz * sxy * syz - 2 * b * eyz * sxy * szx + 2 * a * eyz * cxy * czx + 2 * b * exy * czx * cyz +
             2 * a * exy * szx * cyz + 2 * a * exz * cxy * cyz + (2 * a * std::exp(2 * a * x));

    ret[1] = -2 * a * exy * szx * syz - 2 * b * exz * sxy * syz + 2 * a * exz * cxy * cyz + 2 * b * eyz * cxy * czx +
             2 * a * eyz * sxy * czx + 2 * a * exy * czx * cyz + (2 * a * std::exp(2 * a * y));

    ret[2] = -2 * a * eyz * sxy * syz - 2 * b * exy * szx * syz + 2 * a * exy * czx * cyz + 2 * b * eyz * cxy * czx +
             2 * a * eyz * sxy * czx + 2 * a * exz * cxy * syz + (2 * a * std::exp(2 * a * y));

    ret *= fr;
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
    const double v = Parameters().getParam("viscosity", 1.0, Dune::ValidateNotLess<double>(0.0));
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

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    //				dune_static_assert( dim_ == 2  , "DirichletData_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

} // end namespace TestCase3D

namespace TestCase2D {
static const double pi_factor = M_PI; // controls number of vortices
template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  ~Force() {}

  inline void evaluate(const double time, const DomainType& arg, RangeType& ret) const {
    const double x = arg[0];
    const double y = arg[1];
    const double v = viscosity_;
    const double P = pi_factor;
    const double E = std::exp(-2 * std::pow(P, 2) * viscosity_ * time);
    const double F = std::exp(-4 * std::pow(P, 2) * viscosity_ * time);
    const double S_x = std::sin(P * x);
    const double S_y = std::sin(P * y);
    const double S_2x = std::sin(2 * P * x);
    const double S_2y = std::sin(2 * P * y);
    const double C_x = std::cos(P * x);
    const double C_y = std::cos(P * y);
    //						  ret[0] = - C_x * E * P * ( S_x * E + v * S_y * P )	+ 0.5 * P * F * S_2x;
    //						  ret[1] = - C_y * E * P * ( S_y * E - v * S_x * P )	+ 0.5 * P * F * S_2y;

    // diff
    ret[0] = -2 * C_x * E * P * (v * S_y * P);
    ret[1] = -2 * C_y * E * P * (-v * S_x * P);

    // druck
    ret[0] += 0.5 * P * F * S_2x;
    ret[1] += 0.5 * P * F * S_2y;

    // conv
    ret[0] += -E * E * P * C_x * S_x;
    ret[1] += -E * E * P * S_y * C_y;

    // zeitableitung
    RangeType u;
    VelocityEvaluate(0, 0, time, arg, u);
    ret[0] += (-2 * M_PI * M_PI * v) * u[0];
    ret[1] += (-2 * M_PI * M_PI * v) * u[1];

    ret *= 0;
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
  const double x = arg[0];
  const double y = arg[1];
  const double v = Parameters().getParam("viscosity", 1.0, Dune::ValidateNotLess<double>(0.0));
  const double e_minus_2_t = std::exp(-2 * std::pow(pi_factor, 2) * v * time);

  ret[0] = -1 * std::cos(pi_factor * x) * std::sin(pi_factor * y) * e_minus_2_t;
  ret[1] = std::sin(pi_factor * x) * std::cos(pi_factor * y) * e_minus_2_t;
}

template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert(dim_ == 2, "Pressure_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];
    const double v = Parameters().getParam("viscosity", 1.0, Dune::ValidateNotLess<double>(0.0));
    const double e_minus_4_t = std::exp(-4 * std::pow(pi_factor, 2) * time * v);

    ret[0] = -0.25 * (std::cos(2 * pi_factor * x) + std::cos(2 * pi_factor * y)) * e_minus_4_t;
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
    dune_static_assert(dim_ == 2, "Pressure_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];
    const double v = Parameters().getParam("viscosity", 1.0, Dune::ValidateNotLess<double>(0.0));
    const double e_minus_4_t = std::exp(-4 * std::pow(pi_factor, 2) * time * v);

    ret[0] = 2 * pi_factor * 0.25 * (std::sin(2 * pi_factor * x)) * e_minus_4_t;
    ret[1] = 2 * pi_factor * 0.25 * (std::sin(2 * pi_factor * y)) * e_minus_4_t;
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
    //				dune_static_assert( dim_ == 2 , "DirichletData_Unsuitable_WorldDim");

    const double x = arg[0];
    const double y = arg[1];
    const double v = Parameters().getParam("viscosity", 1.0, Dune::ValidateNotLess<double>(0.0));
    ;
    const double P = pi_factor;
    const double E = std::exp(-2 * std::pow(P, 2) * v * time);
    const double F = std::exp(-4 * std::pow(P, 2) * v * time);
    const double S_x = std::sin(P * x);
    const double S_y = std::sin(P * y);
    const double S_2x = std::sin(2 * P * x);
    const double S_2y = std::sin(2 * P * y);
    const double C_x = std::cos(P * x);
    const double C_y = std::cos(P * y);
    ret[0] = -E * E * P * C_x * S_x;
    ret[1] = -E * E * P * S_y * C_y;
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

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    //				dune_static_assert( dim_ == 2  , "DirichletData_Unsuitable_WorldDim");
    const double x = arg[0];
    const double y = arg[1];
    const double v = Parameters().getParam("viscosity", 1.0, Dune::ValidateNotLess<double>(0.0));
    const double P = pi_factor;
    const double E = std::exp(-2 * std::pow(P, 2) * v * time);
    const double F = std::exp(-4 * std::pow(P, 2) * v * time);
    const double S_x = std::sin(P * x);
    const double S_y = std::sin(P * y);
    const double S_2x = std::sin(2 * P * x);
    const double S_2y = std::sin(2 * P * y);
    const double C_x = std::cos(P * x);
    const double C_y = std::cos(P * y);
    ret[0] = -2 * C_x * E * P * (v * S_y * P);
    ret[1] = -2 * C_y * E * P * (-v * S_x * P);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

} // end namespace TestCase2D

namespace TestCase2DAnimation {
static const double pi_factor = M_PI; // controls number of vortices
template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  ~Force() {}

  inline void evaluate(const double time, const DomainType& arg, RangeType& ret) const { ret = RangeType(0); }
  inline void evaluate(const DomainType& /*arg*/, RangeType& ret) const { assert(false); }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    dune_static_assert(dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    const int id = intersection.boundaryId();
    const double y = arg[1];
    const double f = -4 * (y - 0.5) * (y + 0.5);
    ret = RangeType(0.0);
    if (id == 2) { // bottom
      ret[0] = 0.0;
      ret[1] = 0.0;
    } else if (id == 3) { // right
      ret[0] = std::abs(std::sin(2.0 * M_PI * time)) * f;
      ret[1] = 0.0;
    } else if (id == 4) { // top
      ret[0] = 0.0;
      ret[1] = 0.0;
    } else if (id == 5) { // left
      ret[0] = std::abs(std::sin(2.0 * M_PI * time)) * f;
      ret[1] = 0.0;
    }
  }

  inline void evaluate(const DomainType& arg, RangeType& ret) const { assert(false); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

NULLFUNCTION_TP(Velocity)
NULLFUNCTION_TP(Pressure)
NULLFUNCTION_TP(PressureGradient)
NULLFUNCTION_TP(VelocityLaplace)
NULLFUNCTION_TP(VelocityConvection)

} // end namespace TestCaseAnimation

namespace TrivialTestCase {
template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
  inline void evaluate(const double /*time*/, const DomainType& arg, RangeType& ret) const {
    ret[0] = 0;
    ret[1] = 0;
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
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    ret[0] = 1;
    ret[1] = 0;
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
    ret[0] = 1;
    ret[1] = 0;
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

  void evaluateTime(const double /*time*/, const DomainType& arg, RangeType& ret) const { ret = 0; }

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

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    //				dune_static_assert( dim_ == 2, "DirichletData_Unsuitable_WorldDim");
    ret = RangeType(0);
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

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const { ret = RangeType(0); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
} // end namespace TrivialTestCase

namespace GreenTaylor {
template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  ~Force() {}

  inline void evaluate(const double /*time*/, const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }
  inline void evaluate(const DomainType& /*arg*/, RangeType& ret) const { ret = RangeType(0); }

private:
  const double viscosity_;
  const double alpha_;
  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template <class DomainType, class RangeType>
void VelocityEvaluate(const double /*parameter_a*/, const double /*parameter_d*/, const double time,
                      const DomainType& arg, RangeType& ret) {
  const double x = arg[0];
  const double y = arg[1];
  const double v = Parameters().getParam("viscosity", 1.0, Dune::ValidateNotLess<double>(0.0));
  const double F = std::exp(-8 * std::pow(M_PI, 2) * time);
  const double C1 = std::cos(2 * M_PI * (x + 0.25));
  const double S1 = std::sin(2 * M_PI * (x + 0.25));
  const double S2 = std::sin(2 * M_PI * (y + 0.5));
  const double C2 = std::cos(2 * M_PI * (y + 0.5));

  ret[0] = (-1 / v) * C1 * S2 * F;
  ret[1] = (1 / v) * S1 * C2 * F;
}

template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    const double x = arg[0];
    const double y = arg[1];
    const double v = Parameters().getParam("viscosity", 1.0, Dune::ValidateNotLess<double>(0.0));
    const double F = std::exp(-16 * std::pow(M_PI, 2) * time);
    const double C1 = std::cos(4 * M_PI * (x + 0.25));
    const double C2 = std::cos(4 * M_PI * (y + 0.5));

    ret = (-1 / (4 * v)) * (C1 + C2) * F;

    //							  - 0.920735694 ;
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

} // end namespace TestCase2D_KOKO
namespace DrivenCavity {
NULLFUNCTION(Force)
NULLFUNCTION_TP(VelocityLaplace)
NULLFUNCTION_TP(VelocityConvection)
NULLFUNCTION_TP(Velocity)
NULLFUNCTION_TP(Pressure)

template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    ret = RangeType(0);
    if (arg[1] == 1)
      ret[0] = 1;
  }

  inline void evaluate(const DomainType& arg, RangeType& ret) const { assert(false); }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};
} // end namespace DrivenCavity

namespace NullTest {
NULLFUNCTION(Force)
NULLFUNCTION(DirichletData)
NULLFUNCTION_TP(VelocityLaplace)
NULLFUNCTION_TP(VelocityConvection)
NULLFUNCTION_TP(Velocity)
NULLFUNCTION_TP(Pressure)
NULLFUNCTION_TP(PressureGradient)
}

namespace TimeDisc {
template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
    // laplce
    ret[0] = -2 * std::pow(time, 3.0) * v;
    ret[1] = 0;
    // grad p
    ret[0] += time;
    ret[1] += 1;
    // conv
    ret[0] += std::pow(time, 5.0) * 2 * x * y;
    ret[1] += std::pow(time, 5.0) * y * y;
    // dt u
    ret[0] += std::pow(time, 2.0) * 3 * y * y;
    ret[1] += 2 * time * x;
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
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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

namespace TimeDiscConst {
template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
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
  inline void evaluate(const double /*t*/, const DomainType& arg, RangeType& ret) const {
    const double time = 1.0;
    const double x = arg[0];
    const double y = arg[1];
    const double v = viscosity_;
    // laplce
    ret[0] = -2 * std::pow(time, 3.0) * v;
    ret[1] = 0;
    // grad p
    ret[0] += time;
    ret[1] += 1;
    // conv
    ret[0] += std::pow(time, 5.0) * 2 * x * y;
    ret[1] += std::pow(time, 5.0) * y * y;
    // dt u
    //						  ret[0] += std::pow(time,2.0)*3*y*y;
    //						  ret[1] += 2*time*x;
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
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double /*time*/, const DomainType& arg, RangeType& ret,
                const IntersectionType& intersection) const {
    const double time = 1.0;
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

  void evaluateTime(const double /*time*/, const DomainType& arg, RangeType& ret) const {
    const double time = 1.0;
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

  void evaluateTime(const double /*time*/, const DomainType& arg, RangeType& ret) const {
    const double time = 1.0;
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

  void evaluateTime(const double /*time*/, const DomainType& arg, RangeType& ret) const {
    const double time = 1.0;
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

  void evaluateTime(const double /*time*/, const DomainType& /*arg*/, RangeType& ret) const {
    const double time = 1.0;
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

  void evaluateTime(const double /*time*/, const DomainType& arg, RangeType& ret) const {
    const double time = 1.0;
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
} // end namespace TimeDiscConst

namespace TestCase1D {

template <class R>
double H_eval(const double t, const double sigma, const double mu, const double alpha, const R& x) {
  //(d/dx(d/dx((1/(sqrt(2*pi)*s))*exp(-(x-m)^2/(2*s*s))*exp(-a*t))))
  return (1 / std::sqrt(2 * M_PI * sigma * sigma)) * std::exp(-std::pow(x - mu, 2) / (2 * sigma * sigma)) *
         std::exp(-alpha * t);
}

template <class FunctionSpaceImp>
class Force : public Fem::Function<FunctionSpaceImp, Force<FunctionSpaceImp>> {
public:
  typedef Force<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  Force(const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0)
    : BaseType(space)
    , viscosity_(viscosity)
    , alpha_(alpha) {}

  ~Force() {}

  inline void evaluate(const double time, const DomainType& arg, RangeType& ret) const {
    static const double sigma = 0.1;
    static const double alpha = 3.0;
    static const double mu = 0.5;
    ret = -std::pow(sigma, -4.0) * (arg * arg - 2 * mu * arg + alpha * std::pow(sigma, 4.0) + mu * mu - sigma * sigma) *
          H_eval(time, sigma, mu, alpha, arg);
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
  static const double sigma = 0.1;
  static const double alpha = 3.0;
  static const double mu = 0.5;

  ret = H_eval(time, sigma, mu, alpha, arg);
}

template <class FunctionSpaceImp>
class DirichletData : public Fem::Function<FunctionSpaceImp, DirichletData<FunctionSpaceImp>> {
public:
  typedef DirichletData<FunctionSpaceImp> ThisType;
  typedef Fem::Function<FunctionSpaceImp, ThisType> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  DirichletData(const FunctionSpaceImp& space, const double parameter_a = M_PI / 2.0,
                const double parameter_d = M_PI / 4.0)
    : BaseType(space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~DirichletData() {}

  template <class IntersectionType>
  void evaluate(const double time, const DomainType& arg, RangeType& ret,
                const IntersectionType& /*intersection */) const {
    dune_static_assert(dim_ == 1, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
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
    dune_static_assert(dim_ == 1, "DirichletData_Unsuitable_WorldDim");
    VelocityEvaluate(parameter_a_, parameter_d_, time, arg, ret);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

template <class FunctionSpaceImp, class TimeProviderImp>
class VelocityLaplace
    : public TimeFunction<FunctionSpaceImp, VelocityLaplace<FunctionSpaceImp, TimeProviderImp>, TimeProviderImp> {
public:
  typedef VelocityLaplace<FunctionSpaceImp, TimeProviderImp> ThisType;
  typedef TimeFunction<FunctionSpaceImp, ThisType, TimeProviderImp> BaseType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;

  VelocityLaplace(const TimeProviderImp& timeprovider, const FunctionSpaceImp& space,
                  const double parameter_a = M_PI / 2.0, const double parameter_d = M_PI / 4.0)
    : BaseType(timeprovider, space)
    , parameter_a_(parameter_a)
    , parameter_d_(parameter_d) {}

  ~VelocityLaplace() {}

  void evaluateTime(const double time, const DomainType& arg, RangeType& ret) const {
    dune_static_assert(dim_ == 1, "DirichletData_Unsuitable_WorldDim");
    static const double sigma = 0.1;
    static const double alpha = 3.0;
    static const double mu = 0.5;

    ret = H_eval(time, sigma, mu, alpha, arg);
    ret *= std::pow(sigma, -4.0) * (mu * mu - 2 * mu * arg - sigma * sigma + arg * arg);
  }

private:
  static const int dim_ = FunctionSpaceImp::dimDomain;
  const double parameter_a_;
  const double parameter_d_;
};

NULLFUNCTION_TP(PressureGradient)
NULLFUNCTION_TP(VelocityConvection)
NULLFUNCTION_TP(Pressure)
} // end namespace TestCase1D

} // end namespace NavierStokes
} // end namespace Dune
#endif // TESTDATA_HH

/** Copyright (c) 2012, Felix Tobias Albrecht, Rene Milk
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
