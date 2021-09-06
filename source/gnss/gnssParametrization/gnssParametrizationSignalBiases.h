/***********************************************/
/**
* @file gnssParametrizationSignalBiases.h
*
* @brief Signal biases.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONSIGNALBIASES__
#define __GROOPS_GNSSPARAMETRIZATIONSIGNALBIASES__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationSignalBiases = R"(
\subsection{SignalBiases}\label{gnssParametrizationType:signalBiases}
Each code and phase observation (e.g \verb|C1C| or \verb|L2W|) contains a bias at transmitter/receiver level
\begin{equation}
  [\tau\nu a]_r^s(t) = \dots + \text{bias}[\tau\nu a]^s + \text{bias}[\tau\nu a]_r + \dots
\end{equation}
This class provides the apriori model $\M f(\M x_0)$ of eq. \eqref{gnssParametrizationType:model} only.

The \configFile{inputfileSignalBiasTransmitter/Receiver}{gnssSignalBias} are read
for each receiver and transmitter. The file name is interpreted as a template with
the variables \verb|{prn}| and \verb|{station}| being replaced by the name.
(Infos regarding the variables \verb|{prn}| and \verb|{station}| can be found in
\configClass{gnssTransmitterGeneratorType}{gnssTransmitterGeneratorType} and
\configClass{gnssReceiverGeneratorType}{gnssReceiverGeneratorType} respectively). The files can
be converted with \program{GnssSinexBias2SignalBias}.

The estimation of the biases is complex due to different linear dependencies, which
result in rank deficiencies in the system of normal equations.
For simplification the parametrization for $\Delta\M x$ has been split into:
\configClass{parametrization:codeBiases}{gnssParametrizationType:codeBiases},
\configClass{parametrization:tecBiases}{gnssParametrizationType:tecBiases}, and
\configClass{parametrization:ambiguities}{gnssParametrizationType:ambiguities} (including phase biases).
The file handling on the other hand still remains within this class. Any prior
values for the receiver/transmitter biases are read with the respective \config{inputFileSignalBias}.
All biases for a receiver/transmitter are accumulated and written to the respective \config{outputFileSignalBias}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Signal biases.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationSignalBiases : public GnssParametrizationBase
{
  Gnss                      *gnss;
  std::string                name;
  GnssTransceiverSelectorPtr selectTransmitters, selectReceivers;
  FileName                   fileNameOutTransmitter, fileNameOutReceiver;
  FileName                   fileNameInTransmitter, fileNameInReceiver;

public:
  GnssParametrizationSignalBiases(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
