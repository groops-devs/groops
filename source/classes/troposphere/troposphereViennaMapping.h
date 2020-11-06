/***********************************************/
/**
* @file troposphereViennaMapping.h
*
* @brief Vienna Mapping Functions.
* @see Troposphere
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-04-03
*/
/***********************************************/

#ifndef __GROOPS_TROPOSPHEREVIENNAMAPPING__
#define __GROOPS_TROPOSPHEREVIENNAMAPPING__

// Latex documentation
#ifdef DOCSTRING_Troposphere
static const char *docstringTroposphereViennaMapping = R"(
\subsection{ViennaMapping}\label{troposphereType:viennaMapping}

Tropospheric delays based on the Vienna Mapping Functions 3 (VMF3) model
(Landskron and Boehm 2017, DOI: \href{https://doi.org/10.1007/s00190-017-1066-2}{10.1007/s00190-017-1066-2}).

Hydrostatic and wet mapping function coefficients ($a_h$, $a_w$) and zenith delays (ZHD, ZWD) have to be provided
via \configFile{inputfileVmfCoefficients}{griddedDataTimeSeries}. This file can contain either station-specific data
(see \program{ViennaMappingFunctionStation2File}) or data on a regular global grid
(see \program{ViennaMappingFunctionGrid2File}). In the second case mapping coefficients and zenith delays are
interpolated to the requested coordinates. This includes a height correction that requires approximate meteorological
data provided via \configFile{inputfileGpt}{griddedData}.
)";
#endif

/***********************************************/

#include "troposphere.h"

/***** CLASS ***********************************/

/** @brief Vienna Mapping Function.
* @ingroup troposphereGroup
* @see Troposphere */
class TroposphereViennaMapping : public Troposphere
{
  std::vector<FileName> fileNamesCoefficients;
  FileName              fileNameGpt;
  Double                a_ht, b_ht, c_ht;
  std::vector<Time>     times;
  Double                sampling;
  Matrix                ah, aw, zhd, zwd, gnh, geh, gnw, gew;

  void findIndex(const Time &time, UInt &idx, Double &tau) const;

public:
  TroposphereViennaMapping(Config &config);
 ~TroposphereViennaMapping() {}

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
