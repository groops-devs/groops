/***********************************************/
/**
* @file troposphere.h
*
* @brief Signal delay in the atmosphere.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-04-03
*/
/***********************************************/

#ifndef __GROOPS_TROPOSPHERE__
#define __GROOPS_TROPOSPHERE__

// Latex documentation
#ifdef DOCSTRING_Troposphere
static const char *docstringTroposphere = R"(
\section{Troposphere}\label{troposphereType}
This class provides functions for calculating and estimating
the signal delay in the dry and wet atmosphere.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup troposphereGroup Troposphere
* @brief Signal delay in the atmosphere.
* @ingroup classesGroup
* Functions for calculating and estimating
* the signal delay in the dry and wet atmosphere.
* The interface is given by @ref Troposphere. */
/// @{

/***** TYPES ***********************************/

class Troposphere;
typedef std::shared_ptr<Troposphere> TropospherePtr;

/***** CLASS ***********************************/

/** @brief Signal delay in the atmosphere.
* Functions for calculating and estimating
* the signal delay in the dry and wet atmosphere.
* An Instance of this class can be created by @ref readConfig. */
class Troposphere
{
public:
  /// Destructor.
  virtual ~Troposphere() {}

  /** @brief Init the station list.
  * To speed up the computation and to reduce memory consumption the coordinates
  * of the station (given in TRF [m]) must be set before calling the other functions. */
  virtual void init(const std::vector<std::string> &stationNames, const std::vector<Vector3d> &stationPositions) = 0;

  /** @brief Approx value of the slant delay.
  * @param stationId Station number from the list given at init.
  * @param time Time of the measurement.
  * @param frequency of the electromagnetic signal [Hz].
  * @param azimuth  Azimuth.
  * @param elevation  Elevation.
  * @return delay [m] */
  virtual Double slantDelay(UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const = 0;

  /** @brief Mapping Function of the hydrostatic atmosphere.
  * @param stationId Station number from the list given at init.
  * @param time Time of the measurement.
  * @param frequency of the electromagnetic signal [Hz].
  * @param azimuth  Azimuth.
  * @param elevation  Elevation.
  * @return function value [] */
  virtual Double mappingFunctionHydrostatic(UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const = 0;

  /** @brief Mapping Function of the wet atmosphere.
  * @param stationId Station number from the list given at init.
  * @param time Time of the measurement.
  * @param frequency of the electromagnetic signal [Hz].
  * @param azimuth  Azimuth.
  * @param elevation  Elevation.
  * @return function value [] */
  virtual Double mappingFunctionWet(UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const = 0;

  /** @brief Gradient of the mapping function.
  * @param stationId Station number from the list given at init.
  * @param time Time of the measurement.
  * @param frequency of the electromagnetic signal [Hz].
  * @param azimuth  Azimuth.
  * @param elevation  Elevation.
  * @param[out] dx  Gradient function value in North direction [].
  * @param[out] dy  Gradient function value in East direction []. */
  virtual void mappingFunctionGradient(UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation, Double &dx, Double &dy) const = 0;


  /** @brief Get tropospheric zenith dry/wet delay and dry/wet gradients in North and East directions at a specific time stamp.
  * @param stationId Station number from the list given at init.
  * @param time Time of the measurement.
  * @param frequency of the electromagnetic signal [Hz].
  * @param[out] zenithDryDelay Zenith dry delay [m].
  * @param[out] zenithWetDelay Zenith wet delay [m].
  * @param[out] gradientDryNorth Dry gradient in North direction [m].
  * @param[out] gradientWetNorth Wet gradient in North direction [m].
  * @param[out] gradientDryEast Dry gradient in East direction [m].
  * @param[out] gradientWetEast Wet gradient in East direction [m].
  * @param[out] aDry Dry mapping function coefficient a [].
  * @param[out] aWet Wet mapping function coefficient a []. */
  virtual void getAprioriValues(UInt stationId, const Time &time, Double frequency, Double &zenithDryDelay, Double &zenithWetDelay, Double &gradientDryNorth,
                                Double &gradientWetNorth, Double &gradientDryEast, Double &gradientWetEast, Double &aDry, Double &aWet) const = 0;

  static Double mappingFunction(Double sinE, Double a, Double b, Double c) {return (1+(a/(1+(b/(1+c)))))/(sinE+(a/(sinE+(b/(sinE+c)))));}

  /** @brief creates an derived instance of this class. */
  static TropospherePtr create(Config &config, const std::string &name);

private:
  mutable Time timeRef;
  Matrix coeff;

protected:
  Vector longitude, latitude, height;

  mutable Vector topo;
  mutable Vector ah, aw, bh, bw, ch, cw;
  mutable Vector p, T, Q, dT, la, Tm;
  mutable Vector zhd, zwd;
  mutable Vector gnh, geh, gnw, gew;

  void initEmpiricalCoefficients(const FileName &fileNameGpt, const std::vector<Vector3d> &stationPositions);
  void computeEmpiricalCoefficients(const Time &time) const;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Troposphere.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a troposphere is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] troposphere Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Troposphere */
template<> Bool readConfig(Config &config, const std::string &name, TropospherePtr &troposphere, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
