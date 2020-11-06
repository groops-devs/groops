/***********************************************/
/**
* @file tides.h
*
* @brief Tidal forces.
* Astronomic tides, earth tides, ocean tides, pole tides and so on.
*
* @author Torsten Mayer-Guerr
* @date 2002-12-13.
*
*/
/***********************************************/

#ifndef __GROOPS_TIDES__
#define __GROOPS_TIDES__

// Latex documentation
#ifdef DOCSTRING_Tides
static const char *docstringTides = R"(
\section{Tides}\label{tidesType}
This class computes functionals of the time depending tide potential,
e.g potential, acceleration or gravity gradients.

If several instances of the class are given the results are summed up.
Before summation every single result is multiplicated by a \config{factor}.
To get the difference between two ocean tide models you must choose one factor by 1
and the other by -1. To get the mean of two models just set each factor to 0.5.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "base/sphericalHarmonics.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"

/**
* @defgroup tidesGroup Tides
* @brief Tidal forces.
* @ingroup classesGroup
* The interface is given by @ref Tides.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Tides;
class TidesBase;
typedef std::shared_ptr<Tides> TidesPtr;

/***** CLASS ***********************************/

/** @brief Tidal forces.
* Astronomic tides, earth tides, ocean tides, pole tides and so on.
* An Instance of this class can be created by @ref readConfig. */
class Tides
{
  std::vector<TidesBase*> tides;

public:
  /// Constructor.
  Tides(Config &config, const std::string &name);

  /// Destructor.
  ~Tides();

  /** @brief Tidal potential.
  * @param timeGPS point in time (GPS)
  * @param point Computation point in TRF [m].
  * @param rotEarth CRF -> TRF
  * @param rotation need for computation of polar motion.
  * @param ephemerides ephemerides of sun and moon.
  * @return Potential [@f$m^2/s^2@f$] */
  Double potential(const Time &timeGPS, const Vector3d &point,
                   const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const;

  /** @brief Radial derivative of tidal potential.
  * @param timeGPS point in time (GPS)
  * @param point Computation point in TRF [m].
  * @param rotEarth CRF -> TRF
  * @param rotation need for computation of polar motion.
  * @param ephemerides ephemerides of sun and moon.
  * @return Radial derivative of tidal potential [@f$m/s^2@f$]
  */
  Double radialGradient(const Time &timeGPS, const Vector3d &point,
                        const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const;

  /** @brief Tidal acceleration.
  * @param timeGPS point in time (GPS)
  * @param point Computation point in TRF [m].
  * @param rotEarth CRF -> TRF
  * @param rotation need for computation of polar motion.
  * @param ephemerides ephemerides of sun and moon.
  * @return acceleration in TRF [@f$m/s^2@f$] */
  Vector3d acceleration(const Time &timeGPS, const Vector3d &point,
                        const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const;

  /** @brief Tidal acceleration gradient.
  * @param timeGPS point in time (GPS)
  * @param point Computation point in TRF [m].
  * @param rotEarth CRF -> TRF
  * @param rotation need for computation of polar motion.
  * @param ephemerides ephemerides of sun and moon.
  * @return gradient in TRF [@f$1/s^2@f$] */
  Tensor3d gradient(const Time &timeGPS, const Vector3d &point,
                    const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const;

  /** @brief Deformation.
  * Direct tides creates no deformation.
  * Earth tides and polar tides do not use @a gravity, @a hn, @a ln.
  * @param timeGPS point in time (GPS)
  * @param point station position in TRF [m]
  * @param rotEarth CRF -> TRF
  * @param rotation need for computation of polar motion.
  * @param ephemerides ephemerides of sun and moon.
  * @param gravity local gravity at station [m/s**2]
  * @param hn vertical load love numbers
  * @param ln horizontal load love numbers
  * @return deformation in TRF [m] */
  Vector3d deformation(const Time &timeGPS, const Vector3d &point,
                       const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                       Double gravity, const Vector &hn, const Vector &ln) const;

  /** @brief Deformation.
  * Direct tides creates no deformation.
  * Earth tides and polar tides do not use @a gravity, @a hn, @a ln.
  * @param timeGPS point in time (GPS)
  * @param point station position in TRF [m]
  * @param rotEarth CRF -> TRF
  * @param rotation need for computation of polar motion.
  * @param ephemerides ephemerides of sun and moon.
  * @param gravity local gravity at station [m/s**2]
  * @param hn vertical load love numbers
  * @param ln horizontal load love numbers
  * @param[out] disp series (stationSize x TimeSize) in TRF [m] */
  void deformation(const std::vector<Time> &timeGPS, const std::vector<Vector3d> &point,
                   const std::vector<Rotary3d> &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                   const std::vector<Double> &gravity, const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const;

  /** @brief Tides in terms of SphericalHarmonics.
  * Only gravitational potential of mass distributions of the Earth is considered
  * (Earth tides, ocean tides, without planetary tides). */
  SphericalHarmonics sphericalHarmonics(const Time &timeGPS, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                        UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const;

  /** @brief creates an derived instance of this class. */
  static TidesPtr create(Config &config, const std::string &name) {return TidesPtr(new Tides(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Tides.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and a class without effect is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] tides Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Tides */
template<> Bool readConfig(Config &config, const std::string &name, TidesPtr &tides, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class TidesBase
{
public:
  virtual ~TidesBase() {}

  virtual Double   potential      (const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const;
  virtual Double   radialGradient (const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const;
  virtual Vector3d gravity        (const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const;
  virtual Tensor3d gravityGradient(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const;
  virtual Vector3d deformation    (const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                   Double gravity, const Vector &hn, const Vector &ln) const;
  virtual void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Rotary3d> &rotEarth,
                                   EarthRotationPtr rotation, EphemeridesPtr ephemerides, const std::vector<Double> &gravity, const Vector &hn, const Vector &ln,
                                   std::vector<std::vector<Vector3d>> &disp) const;

  static Matrix deformationMatrix (const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                   const Vector &hn, const Vector &ln, Double GM, Double R, UInt maxDegree);

  virtual SphericalHarmonics sphericalHarmonics(const Time &time, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                                UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const = 0;
};

/***********************************************/

#endif /* __GROOPS_TIDES__ */
