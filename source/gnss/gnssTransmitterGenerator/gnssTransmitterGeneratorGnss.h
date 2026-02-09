/***********************************************/
/**
* @file gnssTransmitterGeneratorGnss.h
*
* @brief Provides a list of GNSS transmitters.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-02-25
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSGNSSTRANSMITTERGENERATORGNSS__
#define __GROOPS_GNSSGNSSTRANSMITTERGENERATORGNSS__

// Latex documentation
#ifdef DOCSTRING_GnssTransmitterGenerator
static const char *docstringGnssTransmitterGeneratorGnss = R"(
\subsection{GNSS}\label{gnssTransmitterGeneratorType:gnss}
A list of satellite PRNs (i.e for GPS: G01, G02, G03, ...) must be provided via
\configFile{inputfileTransmitterList}{stringList}. Satellite system codes follow the
\href{https://files.igs.org/pub/data/format/rinex305.pdf}{RINEX 3 definition}, see \reference{GnssType}{gnssType}.
All input files except \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition},
and \configFile{inputfileReceiverDefinition}{gnssReceiverDefinition} are read for each satellite.
The file name is interpreted as a template with the variable \verb|{prn}| being replaced by the satellite PRN.

Metadata input files (marked with \textbf{*} below) are provided in GROOPS file formats at
\url{https://ftp.tugraz.at/pub/ITSG/groops}. These files are regularly updated.
\begin{itemize}
  \item \configFile{inputfileTransmitterInfo}{platform}\textbf{*}:
        PRN-SVN mapping, antenna offsets and orientations.
        Created via \program{GnssAntex2AntennaDefinition} or \program{PlatformCreate}.
  \item \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition}\textbf{*}:
        Antenna center variations.
        Created via \program{GnssAntex2AntennaDefinition} or \program{GnssAntennaDefinitionCreate}.
  \item \configFile{inputfileSignalDefintion}{gnssReceiverDefinition}\textbf{*}:
        Transmitted signal types.
        Created via \program{GnssReceiverDefinitionCreate} in case you want to define which signal
        types a satellite transmits.
  \item \configFile{inputfileClockFrequencyScale}{instrument}\textbf{*}:
        Scale factor of transmitted signals due to frequency offset/clock drift.
        Can be dreived from broadcast clocks drifts.
  \item \configFile{inputfileOrbit}{instrument}: Converted via \program{Sp3Format2Orbit} or
        output of \program{GnssProcessing}.
  \item \configFile{inputfileAttitude}{instrument}:
        Rotation from body frame to CRF. Created via \program{SimulateStarCameraGnss} or converted via \program{GnssOrbex2StarCamera}.
  \item \configFile{inputfileClock}{instrument}:
        Converted via \program{GnssClockRinex2InstrumentClock} or \program{GnssRinexNavigation2OrbitClock} or
        output of \program{GnssProcessing}.
\end{itemize}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssTransmitterGenerator/gnssTransmitterGenerator.h"

/***** CLASS ***********************************/

/** @brief Provides a list of GNSS transmitters.
* @ingroup gnssTransmitterGeneratorGroup
* @see GnssTransmitterGenerator */
class GnssTransmitterGeneratorGnss : public GnssTransmitterGeneratorBase
{
  std::vector<FileName> fileNamesTransmitterList;
  FileName              fileNameTransmitterInfo, fileNameAntennaDef, fileNameSignalDef;
  FileName              fileNameOrbit, fileNameAttitude, fileNameClock, fileNameScale;
  Bool                  interpolateClock;
  UInt                  interpolationDegree;
  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;

public:
  GnssTransmitterGeneratorGnss(Config &config);
  void init(const std::vector<Time> &times, std::vector<GnssTransmitterPtr> &transmitters) override;
};

/***********************************************/

#endif
