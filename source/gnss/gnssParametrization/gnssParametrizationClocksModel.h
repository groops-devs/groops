/***********************************************/
/**
* @file gnssParametrizationClocksModel.h
*
* @brief Clock errors.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONCLOCKSMODEL__
#define __GROOPS_GNSSPARAMETRIZATIONCLOCKSMODEL__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationClocksModel = R"(
\subsection{ClocksModel}\label{gnssParametrizationType:clocksModel}
This parametrization is an alternative to \configClass{parametrization:clocks}{gnssParametrizationType:clocks}.
Clock errors are estimated epoch-wise for each \configClass{selectTransmitter/Receiver}{platformSelectorType}
and, opposed to \configClass{parametrization:clocks}{gnssParametrizationType:clocks}, are also estimated for epochs
that have no valid observations available (e.g. data gaps).

The clock error of the an epoch can be predicted by the clock error
of the preceding epoch and an unknown clock drift
\begin{equation}
  \Delta t_{i+1} = \Delta t_{i} + t_{drift} dt + \epsilon_i.
\end{equation}
This equation is applied as an additional constraint equation in each epoch
\begin{equation}
  0 = \Delta t_{i+1} - \Delta t_{i} - t_{drift} dt + \epsilon_i.
\end{equation}
The variance $\sigma^2(\epsilon)$ is estimated iteratively by variance component estimation (VCE).
Clock jumps are treated as outliers and are automatically downweighted as described in
\configClass{GnssProcessing:processingStep:estimate}{gnssProcessingStepType:estimate}.

The absolute initial clock error and clock drift cannot be determined if all receiver
and transmitter clocks are estimated together due to their linear dependency.
This linear dependency would lead to a rank deficiency in the normal equation matrix in the same
manner as described in \configClass{parametrization:clocks}{gnssParametrizationType:clocks}.
To circumvent the rank deficiency additional zero-mean constraints are required for the first and last epoch.
The realization of the constraint is done as an additional observation equation in the form
\begin{equation}
  0 = \frac{1}{n_i + n_k} (\sum_i \Delta t^{s_i} + \sum_k \Delta t_{r_k})
\end{equation}
summed over all \configClass{selectTransmitters/ReceiversZeroMean}{platformSelectorType}
with a standard deviation of \config{sigmaZeroMeanConstraint}.

The \file{parameter names}{parameterName} are \verb|<station or prn>:clock::<time>|
and \verb|<station or prn>:clockDrift::|.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Clock errors.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationClocksModel : public GnssParametrizationBase
{
  Gnss                                        *gnss;
  std::string                                  name, nameConstraint;
  PlatformSelectorPtr                          selectTransmitters, selectReceivers;
  PlatformSelectorPtr                          selectTransmittersZeroMean, selectReceiversZeroMean;
  std::vector<Byte>                            selectedTransmitters, selectedReceivers;
  std::vector<Byte>                            selectedTransmittersZeroMean, selectedReceiversZeroMean;
  Double                                       huber, huberPower;
  FileName                                     fileNameReceiver, fileNameTransmitter;
  Bool                                         applyConstraint;
  Double                                       sigmaZeroMean;
  std::vector<std::vector<GnssParameterIndex>> indexTrans, indexRecv;                 // for each trans/recv and epoch
  std::vector<GnssParameterIndex>              indexDriftTrans, indexDriftRecv;       // for each trans/recv
  Vector                                       isMyRank;                              // for each trans
  Vector                                       driftTrans, drift0Trans, clock0Trans;  // for each trans [m/s], [m/s], [m]
  Vector                                       driftRecv,  drift0Recv,  clock0Recv;   // for each recv  [m/s], [m/s], [m]
  Vector                                       sigma0Trans,     sigma0Recv;           // for each trans/recv  [m]
  std::vector<std::vector<Double>>             sigmaEpochTrans, sigmaEpochRecv;       // for each trans/recv and epoch  [m]
  std::vector<Bool>                            isFirstTrans,    isFirstRecv;          // for each trans/recv

public:
  GnssParametrizationClocksModel(Config &config);

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   requirements(GnssNormalEquationInfo &normalEquationInfo, std::vector<UInt> &transCount, std::vector<UInt> &transCountEpoch,
                      std::vector<UInt> &recvCount, std::vector<UInt> &recvCountEpoch) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  void   constraintsEpoch(const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

//***********************************************/

#endif
