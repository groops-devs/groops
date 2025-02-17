/***********************************************/
/**
* @file troposphereSaastamoinen.h
*
* @brief Saastamoinen troposphere model.
* @see Troposphere
*
* @author André Hauschild
* @date 2023-04-27
*/
/***********************************************/

#ifndef __GROOPS_TROPOSPHERESAASTAMOINEN__
#define __GROOPS_TROPOSPHERESAASTAMOINEN__

// Latex documentation
#ifdef DOCSTRING_Troposphere
static const char *docstringTroposphereSaastamoinen = R"(
\subsection{Saastamoinen}\label{troposphereType:Saastamoinen}

Tropospheric delays based on the Saastamoinen model for the zenith hydrostatic delay (Saastamoinen, 1972) and the model by Askne and Nordius for the zenith wet delay (Askne and Nordius, 1987).

Saastamoinen, J. (1972) Atmospheric Correction for the Troposphere and Stratosphere in Radioranging of Satellites, in the Use of Artificial Satellites for Geodesy. Geophys. Monogr. Ser. (edited by Henriksen, S.W., et al.), 15, 247-251.
Askne, J., and Nordius, H. (1987), Estimation of tropospheric delay for microwaves from surface weather data, Radio Sci., 22( 3), 379– 386, doi:10.1029/RS022i003p00379.

)";
#endif

/***********************************************/

#include "troposphere.h"

/***** CLASS ***********************************/

/** @brief Saastamoinen troposphere model.
* @ingroup troposphereGroup
* @see Troposphere */
class TroposphereSaastamoinen : public Troposphere
{
  FileName fileNameGpt;
  Double   a_ht, b_ht, c_ht;

public:
  TroposphereSaastamoinen(Config &config);
 ~TroposphereSaastamoinen() {}

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
