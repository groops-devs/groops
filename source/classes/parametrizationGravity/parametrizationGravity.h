/***********************************************/
/**
* @file parametrizationGravity.h
*
* @brief Parametrization of the gravity field.
*
* @author Torsten Mayer-Guerr
* @date 2001-08-25
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGRAVITY__
#define __GROOPS_PARAMETRIZATIONGRAVITY__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGravity
static const char *docstringParametrizationGravity = R"(
\section{ParametrizationGravity}\label{parametrizationGravityType}
This class gives a parametrization of the time depending gravity field.
Together with the class \configClass{oberservation}{observationType} it will be used
to set up the design matrix in a least squares adjustment.
If multiple parametrizations are given the coefficents in the parameter vector
are sequently appended.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/parameterName.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"

/**
* @defgroup parametrizationGravityGroup ParametrizationGravity.
* @brief Parametrization of the gravity field.
* @ingroup classesGroup
* The interface is given by @ref ParametrizationGravity.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class ParametrizationGravity;
class ParametrizationGravityBase;
typedef std::shared_ptr<ParametrizationGravity> ParametrizationGravityPtr;

/***** CLASS ***********************************/

/** @brief Parametrization of the gravity field.
* Parametrization of the (time variable) gravity field.
* Creates the observation equtions for different functions of the field (design matrix).
* An Instance of this class can be created by @ref readConfig. */
class ParametrizationGravity
{
  UInt _parameterCount;
  std::vector<UInt> index;
  std::vector<ParametrizationGravityBase*> parametrizations;

  void computeIndices();

public:
  /// Constructor.
  ParametrizationGravity(Config &config, const std::string &name);

  /// Destructor.
 ~ParametrizationGravity();

  /** @brief Estimate parameter in the given interval only.
  * Change result of @a parameterCount(), @a parameterName().
  * @return TRUE if parameters are changed */
  Bool setInterval(const Time &timeStart, const Time &timeEnd);

  /** @brief Number of parameters.
  * This is the column count of the design matrix. */
  UInt parameterCount() const {return _parameterCount;}

  /** @brief Name of parameters.
  * The names are appended to @a name. */
  void parameterName(std::vector<ParameterName> &name) const;

  /** @brief Values of the gravity field.
  * Observation equations for values of the gravity field.
  * The values are computed by the convolution of the potential with the inverse @a kernel.
  * Example: If @a kernel is the Stokes function, @a field gives the observation equations for gravity anomalies.
  * @param time Time of observation.
  * @param point Computational point in TRF [m].
  * @param kernel Stellt den Typ der Beobachtung dar (z.B. Stokes-Kern fuer Schwereanomalien).
  * @param A Must be a (sub)matrix with the dimension (1 x parameterCount()). It is filled with the partial derivatives with respect to the parameters. */
  void field(const Time &time, const Vector3d &point, const Kernel &kernel, MatrixSliceRef A) const;

  /** @brief Gravitational potential.
  * Observation equation for the gravitational potential at @a point [m^2/s^2].
  * @param time Time of observation.
  * @param point Computational point in TRF [m].
  * @param A Must be a (sub)matrix with the dimension (1 x parameterCount()). It is filled with the partial derivatives with respect to the parameters. */
  void potential(const Time &time, const Vector3d &point, MatrixSliceRef A) const;

  /** @brief Radial derivative.
  * Observation equations for the radial derivative of the potential in TRF at @a point [m/s^2].
  * @param time Time of observation.
  * @param point Computational point in TRF [m].
  * @param A Must be a (sub)matrix with the dimension (1 x parameterCount()). It is filled with the partial derivatives with respect to the parameters. */
  void radialGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const;

  /** @brief Gravity vector.
  * Three observation equations for the gravity (x,y,z) in TRF at @a point [m/s^2].
  * @param time Time of observation.
  * @param point Computational point in TRF [m].
  * @param A Must be a (sub)matrix with the dimension (3 x parameterCount()). It is filled with the partial derivatives with respect to the parameters. */
  void gravity(const Time &time, const Vector3d &point, MatrixSliceRef A) const;

  /** @brief Gravity Gradient.
  * Six observations equations for gravity gradients (xx,xy,xz,yy,yz,zz) in TRF at @a point [1/s].
  * @param time Time of observation.
  * @param point Computational point in TRF [m].
  * @param A Must be a (sub)matrix with the dimension (6 x parameterCount()). It is filled with the partial derivatives with respect to the parameters. */
  void gravityGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const;

  /** @brief Loading deformation.
  * Three observation equations for loading deformation at a station (x,y,z) in TRF [m].
  * @param time Time of observation.
  * @param point station position in TRF [m]
  * @param gravity local gravity at station [m/s**2]
  * @param hn vertical load love numbers
  * @param ln horizontal load love numbers
  * @param A Must be a (sub)matrix with the dimension (3 x parameterCount()). It is filled with the partial derivatives with respect to the parameters. */
  void deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, MatrixSliceRef A) const;

  /** @brief Conversion of parameter Vector into SphericalHarmonics. */
  SphericalHarmonics sphericalHarmonics(const Time &time, const Vector &x, UInt maxDegree=INFINITYDEGREE) const;

  /** @brief Conversion of parameter Vector into SphericalHarmonics.
  * Additionaly the variances of the parameters can be provided. */
  SphericalHarmonics sphericalHarmonics(const Time &time, const Vector &x, const Vector &sigma2x, UInt maxDegree=INFINITYDEGREE) const;

  /** @brief creates an derived instance of this class. */
  static ParametrizationGravityPtr create(Config &config, const std::string &name) {return ParametrizationGravityPtr(new ParametrizationGravity(config, name));}
};


/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class ParametrizationGravity.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var without parameters is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates ParametrizationGravity */
template<> Bool readConfig(Config &config, const std::string &name, ParametrizationGravityPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Interne Klasse
class ParametrizationGravityBase
{
public:
  virtual ~ParametrizationGravityBase() {}

  virtual Bool setInterval(const Time &/*timeStart*/, const Time &/*timeEnd*/) {return FALSE;}
  virtual UInt parameterCount() const = 0;
  virtual void parameterName(std::vector<ParameterName> &name) const = 0;
  virtual void field          (const Time &time, const Vector3d &point, const Kernel &kernel, MatrixSliceRef A) const = 0;
  virtual void potential      (const Time &time, const Vector3d &point, MatrixSliceRef A) const = 0;
  virtual void radialGradient (const Time &time, const Vector3d &point, MatrixSliceRef A) const = 0;
  virtual void gravity        (const Time &time, const Vector3d &point, MatrixSliceRef A) const = 0;
  virtual void gravityGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const = 0;
  virtual void deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, MatrixSliceRef A) const = 0;

  virtual SphericalHarmonics sphericalHarmonics(const Time &time, const Vector &x, UInt maxDegree) const = 0;
  virtual SphericalHarmonics sphericalHarmonics(const Time &time, const Vector &x, const Vector &sigma2x, UInt maxDegree)  const = 0;
};

/***********************************************/

#endif
