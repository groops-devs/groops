/***********************************************/
/**
* @file gnssParametrizationTransmitterAntennas.h
*
* @brief Antenna center variations.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERANTENNAS__
#define __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERANTENNAS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationTransmitterAntennas = R"(
\subsection{TransmitterAntennas}\label{gnssParametrizationType:transmitterAntennas}
Same as \configClass{receiverAntennas}{gnssParametrizationType:receiverAntennas} but
for transmitting antennas (GNSS satellites).

The \file{parameter names}{parameterName} are
\verb|<antennaName>:<antennaCenterVariations>.<gnssType>::|.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Antenna center variations.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationTransmitterAntennas : public GnssParametrizationBase
{
  Gnss                                        *gnss;
  std::string                                  name;
  PlatformSelectorPtr                          selectTransmitters;
  ParametrizationGnssAntennaPtr                parametrization;
  std::vector<GnssType>                        typesPattern;
  Bool                                         addNonMatchingTypes;
  Bool                                         ignoreSerial;
  std::vector<UInt>                            transmitter2antenna;
  std::vector<std::vector<GnssParameterIndex>> index; // for each antenna and pattern
  std::vector<std::vector<GnssType>>           types; // for each antenna and pattern

public:
  GnssParametrizationTransmitterAntennas(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
};

/***********************************************/

#endif
