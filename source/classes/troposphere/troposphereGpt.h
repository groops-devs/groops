/***********************************************/
/**
* @file troposphereGpt.h
*
* @brief GPT empirical troposphere model.
* @see Troposphere
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2017-02-23
*/
/***********************************************/

#ifndef __GROOPS_TROPOSPHEREGPT__
#define __GROOPS_TROPOSPHEREGPT__

// Latex documentation
#ifdef DOCSTRING_Troposphere
static const char *docstringTroposphereGpt = R"(
\subsection{GPT}\label{troposphereType:gpt}

Tropospheric delays based on the Global Pressure and Temperature 3 (GPT3) model
(Landskron and Boehm 2017, DOI: \href{https://doi.org/10.1007/s00190-017-1066-2}{10.1007/s00190-017-1066-2}).

It is an empirical model derived from the Vienna Mapping Functions 3
(VMF3, see \configClass{viennaMapping}{troposphereType:viennaMapping}) and thus does not require
additional mapping coefficients and zenith delay values.
)";
#endif

/***********************************************/

#include "troposphere.h"

/***** CLASS ***********************************/

/** @brief GPT empirical troposphere model.
* @ingroup troposphereGroup
* @see Troposphere */
class TroposphereGpt : public Troposphere
{
  FileName fileNameGpt;
  Double   a_ht, b_ht, c_ht;

public:
  TroposphereGpt(Config &config);
 ~TroposphereGpt() {}

  void   init(const std::vector<Vector3d> &stationPositions) override;
  Double slantDelay(const Time &time, UInt stationId, Angle azimuth, Angle elevation) const override;
  Double mappingFunctionHydrostatic(const Time &time, UInt stationId, Angle azimuth, Angle elevation) const override;
  Double mappingFunctionWet(const Time &time, UInt stationId, Angle azimuth, Angle elevation) const override;
  void   mappingFunctionGradient(const Time &time, UInt stationId, Angle azimuth, Angle elevation, Double &dx, Double &dy) const override;
  void   getAprioriValues(const Time&time, UInt stationId, Double &zenithDryDelay, Double &zenithWetDelay, Double &gradientDryNorth,
                          Double &gradientWetNorth, Double &gradientDryEast, Double &gradientWetEast, Double &aDry, Double &aWet) const override;
};

/***********************************************/

#endif
