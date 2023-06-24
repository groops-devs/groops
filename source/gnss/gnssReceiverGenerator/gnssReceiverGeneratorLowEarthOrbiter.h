/***********************************************/
/**
* @file gnssReceiverGeneratorLowEarthOrbiter.h
*
* @brief GNSS for Low Earth Orbiter (LEO).
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-02-25
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSRECEIVERGENERATORLOWEARTHORBITER__
#define __GROOPS_GNSSRECEIVERGENERATORLOWEARTHORBITER__

// Latex documentation
#ifdef DOCSTRING_GnssReceiverGenerator
static const char *docstringGnssReceiverGeneratorLowEarthOrbiter = R"(
\subsection{LowEarthOrbiter}\label{gnssReceiverGeneratorType:lowEarthOrbiter}
A single low-Earth orbiting (LEO) satellite with an onboard GNSS receiver.
An apriori orbit is needed as \configFile{inputfileOrbit}{instrument}.
Attitude data must be provided via \configFile{inputfileStarCamera}{instrument}.
If no attitude data is available from the satellite operator,
the star camera data can be simulated by using \program{SimulateStarCamera}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"

/***** CLASS ***********************************/

/** @brief GNSS for Low Earth Orbiter (LEO).
* @ingroup gnssReceiverGeneratorGroup
* @see GnssReceiverGenerator */
class GnssReceiverGeneratorLowEarthOrbiter : public GnssReceiverGeneratorBase
{
  FileName              fileNameStationInfo;
  FileName              fileNameAntennaDef, fileNameReceiverDef, fileNameAccuracyDef;
  FileName              fileNameObs, fileNameOrbit, fileNameStarCamera;
  ExpressionVariablePtr exprSigmaPhase, exprSigmaCode;
  Bool                  integerAmbiguities;
  Double                wavelengthFactor;
  Angle                 elevationCutOff;
  std::vector<GnssType> useType, ignoreType;
  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;
  Bool                  printInfo;
  Double                huber, huberPower;
  Double                codeMaxPosDiff;
  UInt                  minObsCountPerTrack;
  Double                denoisingLambda;
  UInt                  tecWindowSize;
  Double                tecSigmaFactor;
  FileName              fileNameTrackBefore, fileNameTrackAfter;
  GnssReceiverPtr       recv;

public:
  GnssReceiverGeneratorLowEarthOrbiter(Config &config);

  void init(const std::vector<Time> &times, const Time &timeMargin, const std::vector<GnssTransmitterPtr> &transmitters,
            EarthRotationPtr earthRotation, Parallel::CommunicatorPtr comm, std::vector<GnssReceiverPtr> &receivers) override;

  void preprocessing(Gnss *gnss, Parallel::CommunicatorPtr comm) override;

  void simulation(const std::vector<GnssType> &types, NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                  Gnss *gnss, Parallel::CommunicatorPtr comm) override;
};

/***********************************************/

#endif
