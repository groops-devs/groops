/***********************************************/
/**
* @file gnssParametrizationTransmitterGps.h
*
* @brief GPS satellites (transmitter).
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-04-27
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERGPS__
#define __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERGPS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationTransmitter
static const char *docstringGnssParametrizationTransmitterGps = R"(
\subsection{GPS}\label{gnssParametrizationTransmitterType:gps}

Since GPS block IIF satellites are affected by time-variable signal biases
(see Montenbruck et al. 2011, DOI: \href{https://doi.org/10.1007/s10291-011-0232-x}{10.1007/s10291-011-0232-x}),
a \config{biasModel} can be set up for phase observations on the L5 frequency
by defining \configClass{type}{gnssType}=\verb|L5*G| and an appropriate
\configClass{parametrizationAcceleration}{parametrizationAccelerationType}.
The resulting \configFile{outputfileBias}{instrument} contains the time-variable part of the
respective bias evaluated at each epoch, while the constant part is included in
\configFile{outputfileSignalBias}{gnssSignalBias}. In case of PPP the time-variable
part can then be provided via \config{timeVariableBias}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssParametrizationTransmitter.h"

/***** CLASS ***********************************/

/** @brief GPS satellites (transmitter).
* @ingroup gnssParametrizationTransmitterGroup
* @see GnssParametrizationTransmitter */
class GnssParametrizationTransmitterGps : public GnssParametrizationTransmitter
{
public:
  GnssParametrizationTransmitterGps(Config &config);
  virtual ~GnssParametrizationTransmitterGps() {}
  std::string name()   const override {return "gps";}
  GnssType    system() const override {return GnssType::GPS;}
};

/***********************************************/

#endif /* __GROOPS___ */
