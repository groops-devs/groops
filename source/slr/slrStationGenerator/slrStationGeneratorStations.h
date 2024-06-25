/***********************************************/
/**
* @file slrStationGeneratorStations.h
*
* @brief SLR ground station network.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRRECEIVERGENERATORSTATIONNETWORK__
#define __GROOPS_SLRRECEIVERGENERATORSTATIONNETWORK__

// Latex documentation
#ifdef DOCSTRING_SlrStationGenerator
static const char *docstringSlrStationGeneratorStations = R"(
\subsection{Stations}\label{slrStationGeneratorType:stations}
A list of station names must be provided via \configFile{inputfileStationList}{stringList}.
It defines the variable \verb|{station}| for the station specific input files.
The \configFile{inputfileStationInfo}{platform} contains metadata information like station number,
station name and approximate station postion in terrestrial reference frame (TRF)
considering the station eccentricities. They can be created via \program{SinexEccentricties2SlrPlatform}
or \program{PlatformCreate}. The \configFile{inputfileObservations}{instrument} are separate files
for each \verb|{station}|-\verb|{satellite}| pair. They can be converted from CRD
format via \program{Crd2NormalPoints}, CSTG format via \program{Cstg2NormalPoints}
and MERIT II format via \program{Merit2NormalPoints} and \program{Merit2FullRate}.

The apriori observation weighting is defined by the expression \config{accuracy} in $[m]$.
The following variables are defined for each observation from the
\configFile{inputfileObservations}{instrument}: \verb|{residual}|, \verb|{accuracy}|,
\verb|{redundancy}|, \verb|{laserWavelength}|, \verb|{azimut}|, \verb|{elevation}|.
Observations with non-positive accuracies are removed.
This can be used for a rough outlier removal by an expression such as
\config{accuracy} = \verb|if(abs(residual)>30, NAN, accuracy)|.

The effects of loading and tidal deformation on station positions can be corrected for
via \configClass{loadingDisplacement}{gravityfieldType} and
\configClass{tidalDisplacement}{tidesType}, respectively.
Tidal deformations typically include:
\begin{itemize}
  \item \configClass{earthTide}{tidesType:earthTide}: Earth tidal deformations (IERS conventions)
  \item \configClass{doodsonHarmonicTide}{tidesType:doodsonHarmonicTide}: ocean tidal deformations
        (e.g. fes2014b\_n720, \config{minDegree}=\verb|1|)
  \item \configClass{doodsonHarmonicTide}{tidesType:doodsonHarmonicTide}: atmospheric tidal deformation
        (e.g. AOD1B RL06, \config{minDegree}=\verb|1|)
  \item \configClass{poleTide}{tidesType:poleTide}: pole tidal deformations (IERS conventions)
  \item \configClass{poleOceanTide}{tidesType:oceanPoleTide}: ocean pole tidal deformations (IERS conventions)
\end{itemize}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/tides/tides.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "slr/slrStationGenerator/slrStationGenerator.h"

/***** CLASS ***********************************/

/** @brief SLR ground station network.
* @ingroup slrStationGeneratorGroup
* @see SlrStationGenerator */
class SlrStationGeneratorStations : public SlrStationGeneratorBase
{
  FileName              fileNameStationList, fileNameStationInfo;
  FileName              fileNameStationPosition;
  FileName              fileNameObs;
  ExpressionVariablePtr accuracyExpr;
  GravityfieldPtr       gravityfield;
  TidesPtr              tides;
  EphemeridesPtr        ephemerides;
  FileName              deformationName, potentialName;
  Angle                 elevationCutOff;
  UInt                  interpolationDegree;
  std::vector<SlrStationPtr> stations;

public:
  SlrStationGeneratorStations(Config &config);

  void init(const std::vector<Time> &times, const std::vector<SlrSatellitePtr> &satellites,
            EarthRotationPtr earthRotation, std::vector<SlrStationPtr> &stations) override;
};

/***********************************************/

#endif
