/***********************************************/
/**
* @file gnssParametrizationAmbiguities.h
*
* @brief Integer and float ambiguities.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONAMBIGUITIES__
#define __GROOPS_GNSSPARAMETRIZATIONAMBIGUITIES__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationAmbiguities = R"(
\subsection{Ambiguities}\label{gnssParametrizationType:ambiguities}
Sets up an ambiguity parameter for each track and phase observation type.
\begin{equation}
  [L\nu a]_r^s(t) = \dots + \text{bias}[L\nu a]^s + \text{bias}[L\nu a]_r + \lambda[L\nu] N[L\nu a]_r^s
\end{equation}
As the phase observations contain a float bias at transmitter/receiver level, not all ambiguities
are resolvable to integer values. The number of resolvable ambiguities can be increased with
known phase biases read from file via \configClass{parametrization:signalBiases}{gnssParametrizationType:signalBiases}.
In this case, \configClass{estimateTransmitter/ReceiverPhaseBiasTransmitter}{platformSelectorType} should
not be used for the corresponding transmitters and receivers.

In case of GLONASS, the phase biases at receiver level differ between different frequency channels
(frequency division multiple access, FDMA) and for each channel an extra float phase bias is estimated.
With \config{linearGlonassBias} a linear relationship between bias and frequency channel is assumed,
which reduces the number of float bias parameters and increases the number of resolvable integer ambiguities.

The integer ambiguities can be resolved and fixed in
\configClass{GnssProcessing:processingStep:resolveAmbiguities}{gnssProcessingStepType:resolveAmbiguities}.
Resolved integer ambiguities are not estimated as unknown parameters in
\configClass{gnssProcessingStepType:estimate}{gnssProcessingStepType} anymore
and are removed from the system of normal equations.

The estimated phase biases can be written to files in
\configClass{parametrization:signalBiases}{gnssParametrizationType:signalBiases}.

The \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|<station>:phaseBias(<gnssType>)::|,
\item \verb|<prn>:phaseBias(<gnssType>)::|,
\item \verb|<station>.<prn>:ambiguity<index>of<count>(<GnssTypes>)::<track interval>|.
\end{itemize}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssLambda.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Integer and float ambiguities.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationAmbiguities : public GnssParametrizationBase
{
  // float phase biases
  class ParameterRecv
  {
  public:
    GnssReceiverPtr       recv;
    GnssParameterIndex    index;
    std::vector<GnssType> types;
    Matrix                Bias;
  };

  class ParameterTrans
  {
  public:
    GnssTransmitterPtr    trans;
    GnssParameterIndex    index;
    std::vector<GnssType> types;
  };

  class Ambiguity : public GnssAmbiguity
  {
  public:
    std::vector<GnssType> types;
    Matrix                T;          // Matrix to transform ambiguities to observations [cycles -> m]
    Vector                value;      // ambiguities in cycles
    GnssParameterIndex    index;      // index in the parameter vector
    Vector                resolved;
    Bool                  isInteger;

    explicit Ambiguity(GnssTrack *track) : GnssAmbiguity(track) {}
    Vector ambiguities(const std::vector<GnssType> &types) const override;
  };

  class AmbiguityInfo; // to share information between processes

  Gnss                           *gnss;
  std::string                     name;
  PlatformSelectorPtr             selectTransmitters, selectReceivers;
  std::vector<ParameterTrans*>    paraTrans;
  std::vector<ParameterRecv*>     paraRecv;
  Bool                            isLinearBias;

  std::vector<Ambiguity*> getAmbiguities() const;

public:
  GnssParametrizationAmbiguities(Config &config);
 ~GnssParametrizationAmbiguities();

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  Double ambiguityResolve(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                          const std::vector<Byte> &selectedTransmitters, const std::vector<Byte> &selectedReceivers,
                          const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger) override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
};

/***********************************************/

#endif
