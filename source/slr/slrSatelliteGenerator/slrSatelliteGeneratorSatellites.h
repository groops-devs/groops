/***********************************************/
/**
* @file slrSatelliteGeneratorSatellites.h
*
* @brief Provides a list of SLR satellites.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRSLRSATELLITEGENERATORSATELLITES__
#define __GROOPS_SLRSLRSATELLITEGENERATORSATELLITES__

// Latex documentation
#ifdef DOCSTRING_SlrSatelliteGenerator
static const char *docstringSlrSatelliteGeneratorSatellites = R"(
\subsection{Satellites}\label{slrSatelliteGeneratorType:satellites}
A list of satellite names must be provided via \configFile{inputfileSatelliteList}{stringList}.
The other input files are read for each satellite, where the file name is interpreted as a template
with the variable \verb|{satellite}| being replaced by the satellite name from list.
The \configFile{inputfileSatelliteInfo}{platform} contains information about laser retro-reflector,
optical reference point, retro-reflector orientation, range corrections and center of mass.
It can be created via \program{PlatformCreate}.
If \configFile{inputfileAttitude}{instrument} ist not provided an orbit reference frame
(along, cross, nearly nadir) is assumed.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "slr/slrSatelliteGenerator/slrSatelliteGenerator.h"

/***** CLASS ***********************************/

/** @brief Provides a list of SLR satellites.
* @ingroup slrSatelliteGeneratorGroup
* @see SlrSatelliteGenerator */
class SlrSatelliteGeneratorSatellites : public SlrSatelliteGeneratorBase
{
  std::vector<FileName> fileNamesSatelliteList;
  FileName              fileNameSatelliteInfo;
  FileName              fileNameOrbit, fileNameAttitude;
  UInt                  interpolationDegree;

public:
  SlrSatelliteGeneratorSatellites(Config &config);
  void init(const std::vector<Time> &times, std::vector<SlrSatellitePtr> &satellites) override;
};

/***********************************************/

#endif
