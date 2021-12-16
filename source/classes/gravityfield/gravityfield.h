/***********************************************/
/**
* @file gravityfield.h
*
* @brief Function values of (time variable) gravity fields.
*
* @author Torsten Mayer-Guerr
* @date 2001-08-15
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELD__
#define __GROOPS_GRAVITYFIELD__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfield = R"(
\section{Gravityfield}\label{gravityfieldType}
This class computes functionals of the time depending gravity field,
e.g potential, gravity anomalies or gravity gradients.

If several instances of the class are given the results are summed up.
Before summation every single result is multiplicated by a \config{factor}.
To subtract a normal field like GRS80 from a potential
to get the disturbance potential you must choose one factor by 1
and the other by -1. To get the mean of two fields just set each factor to 0.5.

Some of the instances gives also information about the accuracy.
The variance of the result (sum) is computed by means of variance propagation.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"

/**
* @defgroup gravityfieldGroup GravityField
* @brief Function values of (time variable) gravity fields.
* @ingroup classesGroup
* The interface is given by @ref Gravityfield.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Gravityfield;
class GravityfieldBase;
typedef std::shared_ptr<Gravityfield> GravityfieldPtr;

/***** CLASS ***********************************/

/** @brief Function values of (time variable) gravity fields.
* Calculate at the computation point functionals of (time variable) gravity fields.
* An Instance of this class can be created by @ref readConfig. */
class Gravityfield // : public GravityfieldBase
{
public:
  /// Constructor.
  Gravityfield(Config &config, const std::string &name);

  /// Destructor.
  ~Gravityfield();

  /** @brief Function values of the gravity field.
  * Convolution of an inverse kernel with the potential.
  * Example: Computes gravity anomalies with the stokes kernel.
  * @param time If time==Time(), only the static part will be computed.
  * @param point computation point in an Earth fixed reference system [m].
  * @param kernel inverse kernel as differential operator. */
  Double field(const Time &time, const Vector3d &point, const Kernel &kernel) const;

  /** @brief Potential.
  * @param time If time==Time(), only the static part will be computed.
  * @param point computation point in an Earth fixed reference system [m].
  * @return potential [m^2/s^2]. */
  Double potential(const Time &time, const Vector3d &point) const;

  /** @brief Radial gradient of the potential.
  * @f$ \frac{\partial V}{\partial r} @f$.
  * @param time If time==Time(), only the static part will be computed.
  * @param point computation point in an Earth fixed reference system [m].
  * @return radial gradient [m/s^2]. */
  Double radialGradient(const Time &time, const Vector3d &point) const;

  /** @brief Gravity vector.
  * @f$ \nabla V @f$.
  * @param time If time==Time(), only the static part will be computed.
  * @param point computation point in an Earth fixed reference system [m].
  * @return Result is given in an Earth fixed system [m/s^2]. */
  Vector3d gravity(const Time &time, const Vector3d &point) const;

  /** @brief Gravity gradient.
  * @f$ \nabla\nabla V @f$
  * @param time If time==Time(), only the static part will be computed.
  * @param point computation point in an Earth fixed reference system [m].
  * @return Result is given in an Earth fixed system [1/s^2]. */
  Tensor3d gravityGradient(const Time &time, const Vector3d &point) const;

  /** @brief Deformation due to loading.
  * @param  time point in time
  * @param  point station position in TRF [m]
  * @param  gravity local gravity at station [m/s**2]
  * @param  hn vertical load love numbers
  * @param  ln horizontal load love numbers
  * @return deformation in TRF [m] */
  Vector3d deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const;

  /** @brief Deformation due to loading.
  * @param  time points in time
  * @param  point station positions in TRF [m]
  * @param  gravity local gravity at stations [m/s**2]
  * @param  hn vertical load love numbers
  * @param  ln horizontal load love numbers
  * @param[out] disp series (stationSize x TimeSize) in TRF [m] */
  void deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity, const Vector &hn, const Vector &ln, std::vector< std::vector<Vector3d> > &disp) const;

  /** @brief Conversion into a series of spherical harmonics.
  * If @a time==Time(), only the static part will be computed. */
  SphericalHarmonics sphericalHarmonics(const Time &time, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const;

  /** @brief Variance-Covariance-Matrix.
  * The result is a full covariance matrix or a vector containing only the variances
  * of a spherical harmonics expansion.
  * The coefficients are given in a degree wise sequence with alternating sin, cos:
  * ..., c20,c21,s21,c22,s22,..., beginning with c00. */
  Matrix sphericalHarmonicsCovariance(const Time &time, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const;

  /** @brief Variance-Covariance-Matrix.
  * Variance-Covariance-Matrix of gravity field functionals at a list of grid points. */
  Matrix variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel) const;

  /** @brief Variance of gravity field functional.
  * Variance of a gravity field functional at one point. */
  Double variance(const Time &time, const Vector3d &point, const Kernel &kernel) const;

  /** @brief Covariance between gravity field functional at two different points. */
  Double covariance(const Time &time, const Vector3d &point1, const Vector3d &point2, const Kernel &kernel) const;

  /** @brief creates an derived instance of this class. */
  static GravityfieldPtr create(Config &config, const std::string &name) {return GravityfieldPtr(new Gravityfield(config, name));}

private:
  std::vector<GravityfieldBase*> gravityfield;

  void variance(const Time &/*time*/, const std::vector<Vector3d> &/*point*/, const Kernel &/*kernel*/, Matrix &/*D*/) const {}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Gravityfield.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and a class with zero gravity is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] gravityfield Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Gravityfield */
template<> Bool readConfig(Config &config, const std::string &name, GravityfieldPtr &gravityfield, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class GravityfieldBase
{
public:
virtual ~GravityfieldBase() {}

virtual Double   field          (const Time &time, const Vector3d &point, const Kernel &kernel) const;
virtual Double   potential      (const Time &time, const Vector3d &point) const = 0;
virtual Double   radialGradient (const Time &time, const Vector3d &point) const = 0;
virtual Vector3d gravity        (const Time &time, const Vector3d &point) const = 0;
virtual Tensor3d gravityGradient(const Time &time, const Vector3d &point) const = 0;
virtual Vector3d deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const = 0;
virtual void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                 const Vector &hn, const Vector &ln, std::vector< std::vector<Vector3d> > &disp) const = 0;

static Matrix deformationMatrix(const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                const Vector &hn, const Vector &ln, Double GM, Double R, UInt maxDegree);

virtual SphericalHarmonics sphericalHarmonics(const Time &time, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const = 0;
virtual Matrix   sphericalHarmonicsCovariance(const Time &time, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const;

virtual void   variance  (const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const=0;
virtual Double variance  (const Time &time, const Vector3d &point, const Kernel &kernel) const;
virtual Double covariance(const Time &time, const Vector3d &point1, const Vector3d &point2, const Kernel &kernel) const;
};

/***********************************************/

#endif /* __GROOPS_GRAVITYFIELD__ */

