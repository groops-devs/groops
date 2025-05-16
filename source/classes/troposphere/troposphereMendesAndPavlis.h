/***********************************************/
/**
* @file troposphereMendesAndPavlis.h
*
* @brief Atmospheric Delay Correction by Mendes and Pavlis.
* @see Troposphere
*
* @author Torsten Mayer-Guerr
* @author Barbara Suesser-Rechberger
* @date 2022-09-24
*/
/***********************************************/

#ifndef __GROOPS_TROPOSPHEREMENDESANDPAVLIS__
#define __GROOPS_TROPOSPHEREMENDESANDPAVLIS__

// Latex documentation
#ifdef DOCSTRING_Troposphere
static const char *docstringTroposphereMendesAndPavlis = R"(
\subsection{MendesAndPavlis}\label{troposphereType:mendesAndPavlis}

Tropospheric delays based on the Mendes-Pavlis model that employs meteorological data.
(Mendes et al. (2002), \href{https://doi.org/10.1029/2001GL014394}{10.1029/2001GL014394} and
Mendes and Pavlis (2004), \href{https://doi.org/10.1029/2004GL020308}{110.1029/2004GL020308})

The meteorological data have to be provided via \configFile{inputfileStationMeteorology}{instrument}.
This file contains the temperature, air pressure and humidity and must be first generated using the
programs \program{Crd2NormalPoints}, \program{Cstg2NormalPoints}, \program{Merit2NormalPoints} or \program{Merit2FullRate}.
)";
#endif

/***********************************************/

#include "files/fileInstrument.h"
#include "troposphere.h"
#include "base/polynomial.h"

/***** CLASS ***********************************/

/** @brief Vienna Mapping Function.
* @ingroup troposphereGroup
* @see Troposphere */
class TroposphereMendesAndPavlis : public Troposphere
{
  FileName                       fileNameMeteo, fileNameGeoid;
  Vector                         latitude, height, f;
  std::vector<std::vector<Time>> times;       // for each station
  std::vector<Matrix>            meteorology; // for each station: times x (temperature, pressure, waterVapor)

  Matrix interpolate(UInt stationId, const Time &time, const_MatrixSliceRef A) const;

public:
  TroposphereMendesAndPavlis(Config &config);
 ~TroposphereMendesAndPavlis() {}

  void   init(const std::vector<std::string> &stationNames, const std::vector<Vector3d> &stationPositions) override;
  Double slantDelay                (UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const override;
  Double mappingFunctionHydrostatic(UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const override;
  Double mappingFunctionWet        (UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const override;
  void   mappingFunctionGradient   (UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation, Double &dx, Double &dy) const override;
  void   getAprioriValues          (UInt stationId, const Time &time, Double frequency, Double &zenithDryDelay, Double &zenithWetDelay, Double &gradientDryNorth,
                                    Double &gradientWetNorth, Double &gradientDryEast, Double &gradientWetEast, Double &aDry, Double &aWet) const override;
};

/***********************************************/

inline TroposphereMendesAndPavlis::TroposphereMendesAndPavlis(Config &config)
{
  try
  {
    readConfig(config, "inputfileStationMeteorology", fileNameMeteo, Config::MUSTSET, "meteorology_{loopTime:%y}.{station}.dat",  "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void TroposphereMendesAndPavlis::init(const std::vector<std::string> &stationNames, const std::vector<Vector3d> &stationPositions)
{
  try
  {
    latitude = Vector(stationPositions.size());
    height   = Vector(stationPositions.size());
    f        = Vector(stationPositions.size());
    Ellipsoid ellipsoid;
    for(UInt stationId=0; stationId<stationPositions.size(); stationId++)
    {
      Angle L, B;
      ellipsoid(stationPositions.at(stationId), L, B, height(stationId));
      latitude(stationId) = B;
      f(stationId)        = (1. - 0.00266*std::cos(2.*B) - 0.00028e-3 * height(stationId));
    }

    // read station wise meteorology data
    meteorology.resize(stationNames.size());
    times.resize(stationNames.size());
    for(UInt stationId=0; stationId<stationPositions.size(); stationId++)
    {
      VariableList varList;
      varList.setVariable("station", stationNames.at(stationId));
      MeteorologicalArc arc = InstrumentFile::read(fileNameMeteo(varList));
      times.at(stationId)       = arc.times();
      meteorology.at(stationId) = arc.matrix().column(1, 3); // temperature [K], pressure [Pa], humidity [%]

      // Calculation of the water vapour pressure [Pa] from relative humidity [%]:
      // Equation source: Marini, J. W. and Murray Jr, C. (1973).
      // Correction of laser range tracking data for atmospheric refraction at elevations above 10 degrees.
      for(UInt i=0; i<meteorology.at(stationId).rows(); i++)
      {
        const Double T = meteorology.at(stationId)(i, 0)-273.15; // temperature [C°]
        meteorology.at(stationId)(i, 2) *= 6.11*std::pow(10, 7.5*T/(237.3+T));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Matrix TroposphereMendesAndPavlis::interpolate(UInt stationId, const Time &time, const_MatrixSliceRef A) const
{
  try
  {
    auto closestTimeIter = std::min_element(times.at(stationId).begin(), times.at(stationId).end(), [&](const Time &t1, const Time &t2)
                          {return std::fabs((t1-time).seconds()) < std::fabs((t2-time).seconds());});

    Time closestTime = *closestTimeIter;

    if(times.at(stationId).size() >= 2)
    {
      // Linear interpolation should be possible but shall be only done if the distance between the points is less than one hour
      Double distanceToClosestPoint = closestTime.seconds() - time.seconds();
      auto closestTimeNeighbourIter = closestTimeIter;

      if(distanceToClosestPoint < 0)
      {
        if(closestTimeNeighbourIter < (times.at(stationId).end()-1))
          closestTimeNeighbourIter++;
      }
      else
      {
        if(closestTimeNeighbourIter > times.at(stationId).begin())
          closestTimeNeighbourIter--;
      }

      if((std::fabs(closestTime.seconds() - (*closestTimeNeighbourIter).seconds()) < 60*60) && (closestTimeNeighbourIter != closestTimeIter))
      {
        Polynomial polynomial(times.at(stationId), 1, FALSE);
        Matrix interpolatedValues = polynomial.interpolate({time}, A);

        // Check if interpolation was done
        if(!std::isnan(sum(interpolatedValues)))
        {
          // Interpolation was successful
          return interpolatedValues;
        }
        else
        {
          // Interpolation was not done, element which is closest to the given time is returned
          return A.row(std::distance(times.at(stationId).begin(), closestTimeIter));
        }
      }
      else
      {
        // Interpolation was not done, element which is closest to the given time is returned
        return A.row(std::distance(times.at(stationId).begin(), closestTimeIter));
      }
    }
    else
    {
      // Interpolation was not done, element which is closest to the given time is returned
      return A.row(std::distance(times.at(stationId).begin(), closestTimeIter));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Double TroposphereMendesAndPavlis::slantDelay(UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const
{
  try
  {
    Double zenithDryDelay, zenithWetDelay, gradientDryNorth, gradientWetNorth, gradientDryEast, gradientWetEast, aDry, aWet;
    getAprioriValues(stationId, time, frequency, zenithDryDelay, zenithWetDelay, gradientDryNorth, gradientWetNorth, gradientDryEast, gradientWetEast, aDry, aWet);
    return (zenithDryDelay + zenithWetDelay) * mappingFunctionHydrostatic(stationId, time, frequency, azimuth, elevation);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Double TroposphereMendesAndPavlis::mappingFunctionHydrostatic(UInt stationId, const Time &time, Double /*frequency*/, Angle /*azimuth*/, Angle elevation) const
{
  try
  {
    // FCULa source: IERS 2010 conventions, Table 9.1
    constexpr Double a10 = 12100.8e-7;
    constexpr Double a11 =  1729.5e-9;
    constexpr Double a12 =   319.1e-7;
    constexpr Double a13 = -1847.8e-11;
    constexpr Double a20 = 30496.5e-7;
    constexpr Double a21 =   234.6e-8;
    constexpr Double a22 =  -103.5e-6;
    constexpr Double a23 =  -185.6e-10;
    constexpr Double a30 =  6877.7e-5;
    constexpr Double a31 =   197.2e-7;
    constexpr Double a32 =  -345.8e-5;
    constexpr Double a33 =   106.0e-9;

    const Double T    = interpolate(stationId, time, meteorology.at(stationId).column(0))(0, 0)-273.15; // convert temperature in [°C]
    const Double cosL = std::cos(latitude(stationId));
    const Double a1   = a10 + a11*T + a12*cosL + a13*height(stationId);
    const Double a2   = a20 + a21*T + a22*cosL + a23*height(stationId);
    const Double a3   = a30 + a31*T + a32*cosL + a33*height(stationId);

    return mappingFunction(std::sin(elevation), a1, a2, a3);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Double TroposphereMendesAndPavlis::mappingFunctionWet(UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const
{
  return mappingFunctionHydrostatic(stationId, time, frequency, azimuth, elevation);
}

/***********************************************/

inline void TroposphereMendesAndPavlis::mappingFunctionGradient(UInt /*stationId*/, const Time &/*time*/, Double /*frequency*/, Angle /*azimuth*/, Angle /*elevation*/, Double &dx, Double &dy) const
{
  try
  {
    throw(Exception("not implemented yet"));
    dx = dy = 0.;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void TroposphereMendesAndPavlis::getAprioriValues(UInt stationId, const Time &time, Double frequency,
                                                         Double &zenithDryDelay, Double &zenithWetDelay,
                                                         Double &gradientDryNorth, Double &gradientWetNorth, Double &gradientDryEast, Double &gradientWetEast,
                                                         Double &aDry, Double &aWet) const
{
  try
  {
    // Values source: Mendes, V. B., and Pavlis, E. C. (2004), High-accuracy zenith delay prediction at optical wavelengths, Geophys. Res. Lett., 31, L14602
    constexpr Double k0     = 238.0185;  // [um^-2]
    constexpr Double k1     = 19990.975; // [um^-2]
    constexpr Double k2     = 57.362;    // [um^-2]
    constexpr Double k3     = 579.55174; // [um^-2]
    constexpr Double xc     = 375;       // [ppm] carbon dioxide content according to IAG recommendations
    constexpr Double cco2   = 1 + 0.534e-6 * (xc - 450); // carbon dioxide coefficient

    constexpr Double omega0 =  295.235;  // [-]
    constexpr Double omega1 =  2.6422;   // [um^2]
    constexpr Double omega2 = -0.032380; // [um^4]
    constexpr Double omega3 =  0.004028; // [um^6]

    const Double sigma2 = std::pow(1e-6*frequency/LIGHT_VELOCITY, 2); // [1/um^2]
    const Double fh   = 1e-2 * (k1*(k0+sigma2)/std::pow((k0-sigma2), 2) + k3*(k2+sigma2)/std::pow(k2-sigma2, 2)) * cco2;
    const Double fnh  = 0.003101 * (omega0 + 3*omega1*sigma2 + 5*omega2*std::pow(sigma2, 2) + 7*omega3*std::pow(sigma2, 3));

    const Matrix M = interpolate(stationId, time, meteorology.at(stationId)).column(1, 2); // pressure [Pa], waterVaporPressure [Pa]
    zenithDryDelay = 0.00002416579*fh/f(stationId) * M(0, 0)/*pressure*/;
    zenithWetDelay = (5.316e-6*fnh - 3.759e-6*fh)/f(stationId) * M(0, 1)/*waterVaporPressure*/;

    gradientDryNorth = gradientWetNorth = gradientDryEast = gradientWetEast = 0.;
    aDry = aWet = 0.;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
