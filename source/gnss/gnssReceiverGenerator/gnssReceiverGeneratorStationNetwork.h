/***********************************************/
/**
* @file gnssReceiverGeneratorStationNetwork.h
*
* @brief GNSS ground station network.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-02-25
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSRECEIVERGENERATORSTATIONNETWORK__
#define __GROOPS_GNSSRECEIVERGENERATORSTATIONNETWORK__

// Latex documentation
#ifdef DOCSTRING_GnssReceiverGenerator
static const char *docstringGnssReceiverGeneratorStationNetwork = R"(
\subsection{StationNetwork}\label{gnssReceiverGeneratorType:stationNetwork}
A network of GNSS ground stations is defined via \configFile{inputfileStationList}{stringTable}.
Each line can contain more than one station. The first station in each line for which \configFile{inputfileObservations}{instrument}
exists and contains enough observations is used for the processing.
All input files except \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition},
\configFile{inputfileReceiverDefinition}{gnssReceiverDefinition}, and
\configFile{inputfileAccuracyDefinition}{gnssAntennaDefinition} are read for each station.
The file name is interpreted as a template with the variable \verb|{station}| being replaced by the station name.

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
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"

/***** CLASS ***********************************/

/** @brief GNSS ground station network.
* @ingroup gnssReceiverGeneratorGroup
* @see GnssReceiverGenerator */
class GnssReceiverGeneratorStationNetwork : public GnssReceiverGeneratorBase
{
  FileName              fileNameStationList, fileNameStationInfo;
  FileName              fileNameAntennaDef, fileNameReceiverDef, fileNameAccuracyDef;
  FileName              fileNameStationPosition, fileNameClock, fileNameObs;
  UInt                  maxStationCount;
  GravityfieldPtr       gravityfield;
  TidesPtr              tides;
  EphemeridesPtr        ephemerides;
  FileName              deformationName, potentialName;
  std::vector<GnssType> useType, ignoreType;
  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;
  Angle                 elevationCutOff;
  Bool                  printInfo;
  UInt                  minObsCountPerTrack;
  Angle                 elevationTrackMinimum;
  Double                minEstimableEpochsRatio;
  Double                huber, huberPower;
  Double                codeMaxPosDiff;
  Double                denoisingLambda;
  UInt                  tecWindowSize;
  Double                tecSigmaFactor;
  FileName              fileNameTrackBefore, fileNameTrackAfter;
  std::vector<GnssReceiverPtr> receivers;

public:
  GnssReceiverGeneratorStationNetwork(Config &config);

  void init(const std::vector<Time> &times, const Time &timeMargin, const std::vector<GnssTransmitterPtr> &transmitters,
            EarthRotationPtr earthRotation, Parallel::CommunicatorPtr comm, std::vector<GnssReceiverPtr> &receivers) override;

  void preprocessing(Gnss *gnss, Parallel::CommunicatorPtr comm) override;

  void simulation(const std::vector<GnssType> &types, NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                  Gnss *gnss, Parallel::CommunicatorPtr comm) override;
};

/***********************************************/

#endif
