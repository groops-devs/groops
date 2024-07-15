/***********************************************/
/**
* @file troposphereHopfield.h
*
* @brief Hopfield troposphere model.
* @see Troposphere
*
* @author André Hauschild
* @date 2023-04-27
*/
/***********************************************/

#ifndef __GROOPS_TROPOSPHEREHOPFIELD__
#define __GROOPS_TROPOSPHEREHOPFIELD__

// Latex documentation
#ifdef DOCSTRING_Troposphere
static const char *docstringTroposphereHopfield = R"(
\subsection{Hopfield}\label{troposphereType:Hopfield}

Tropospheric delays based on a simplified version of the Hopfield model (Hopfield, 1969) assuming identical and constant meteorological parameters at all stations.

Hopfield, H.S. (1969) Two-Quartic Tropospheric Refractivity Profile for Correcting Satellite Data. Journal of Geophysical Research, 74, 4487-4499
http://dx.doi.org/10.1029/JC074i018p04487

)";
#endif

/***********************************************/

#include "troposphere.h"

/***** CLASS ***********************************/

/** @brief Hopfield troposphere model.
* @ingroup troposphereGroup
* @see Troposphere */
class TroposphereHopfield : public Troposphere
{
  FileName fileNameGpt;
  Double   a_ht, b_ht, c_ht;

public:
  TroposphereHopfield(Config &config);
 ~TroposphereHopfield() {}

  void   init(const std::vector<std::string> &stationNames, const std::vector<Vector3d> &stationPositions) override;
  Double slantDelay                (UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const override;
  Double mappingFunctionHydrostatic(UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const override;
  Double mappingFunctionWet        (UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const override;
  void   mappingFunctionGradient   (UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation, Double &dx, Double &dy) const override;
  void   getAprioriValues          (UInt stationId, const Time &time, Double frequency, Double &zenithDryDelay, Double &zenithWetDelay, Double &gradientDryNorth,
                                    Double &gradientWetNorth, Double &gradientDryEast, Double &gradientWetEast, Double &aDry, Double &aWet) const override;
};

/***********************************************/

#endif
