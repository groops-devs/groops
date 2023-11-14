/***********************************************/
/**
* @file gnssParametrizationClocks.h
*
* @brief Clocks errors.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONCLOCKS__
#define __GROOPS_GNSSPARAMETRIZATIONCLOCKS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationClocks = R"(
\subsection{Clocks}\label{gnssParametrizationType:clocks}
Clock errors are estimated epoch-wise for each \configClass{selectTransmitter/Receiver}{platformSelectorType}.
No clock errors are estimated if no valid observations are available (e.g. data gaps in the observations).

These parameters are lineary dependent and would lead to a rank deficiency in the normal equation
matrix. To circumvent this issue, the estimation requires an additional zero-mean constraint added in each epoch.
This is realized with an additional observation equation
\begin{equation}
 0 = \sum_i \delta t^{s_i} + \sum_k \delta t_{r_k}
\end{equation}
summed over all \configClass{selectTransmitters/ReceiversZeroMean}{platformSelectorType}
with a standard deviation of \config{sigmaZeroMeanConstraint}.

The \file{parameter names}{parameterName} are \verb|<station or prn>:clock::<time>|.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Clocks errors.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationClocks : public GnssParametrizationBase
{
  Gnss                                        *gnss;
  std::string                                  name, nameConstraint;
  PlatformSelectorPtr                          selectTransmitters, selectReceivers;
  PlatformSelectorPtr                          selectTransmittersZeroMean, selectReceiversZeroMean;
  std::vector<Byte>                            selectedTransmitters, selectedReceivers;
  std::vector<Byte>                            selectedTransmittersZeroMean, selectedReceiversZeroMean;
  FileName                                     fileNameReceiver, fileNameTransmitter;
  Bool                                         applyConstraint;
  Double                                       sigmaZeroMean;
  std::vector<std::vector<GnssParameterIndex>> indexTrans, indexRecv; // for each trans/recv and epoch
  std::vector<std::vector<Double>>             x0Trans, x0Recv;       // for each trans/recv and epoch

public:
  GnssParametrizationClocks(Config &config);

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
