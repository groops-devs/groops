/***********************************************/
/**
* @file gnssParametrizationCodeBiases.h
*
* @brief Code biases.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONCODEBIASES__
#define __GROOPS_GNSSPARAMETRIZATIONCODEBIASES__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationCodeBiases = R"(
\subsection{CodeBiases}\label{gnssParametrizationType:codeBiases}
Each code observation (e.g \verb|C1C| or \verb|C2W|) contains a bias at transmitter/receiver level
\begin{equation}
  [C\nu a]_r^s(t) = \dots + \text{bias}[C\nu a]^s + \text{bias}[C\nu a]_r + \dots
\end{equation}
The code biases cannot be estimated together with clock errors and ionospheric delays in an absolute sense
as rank deficiencies will occur in the system of normal equations. Therefore, the biases are not initialized and set up
as parameters directly but only estimable linear combinations are parametrized.

The basic idea is to set up simplified normal equations with the biases,
clock and STEC parameters of one single receiver or transmitter,
eliminate clock and STEC parameters and perform an eigen value decomposition
of the normal equation matrix
\begin{equation}
  \M N = \M Q \M\Lambda \M Q^T.
\end{equation}
Instead of estimating the original bias parameter $\M x$ a transformed set $\bar{\M x}$
is introduced:
\begin{equation}
  \bar{\M x} = \M Q^T \M x.
\end{equation}
The new parameters corresponding to eigen values $\lambda>0$ are estimable,
the others are left out (set to zero). The missing linear combinations,
which depend on the STEC parameters, can be added with
\configClass{parametrization:tecBiases}{gnssParametrizationType:tecBiases}.

Additional rank deficiencies may also occur when biases of transmitters and receivers are estimated together.
The minimum norm nullspace (also via eigen value decomposition)
is formulated as zero constraint equations and added with a standard deviation of \config{sigmaZeroMeanConstraint}.

In case of GLONASS the code biases at receiver level can differ between different frequency channels
(frequency division multiple access, FDMA) and for each channel an extra code bias is estimated.
With \config{linearGlonassBias} a linear relationship between bias and frequency channel is assumed,
which reduces the number of bias parameters.

The estimated biases can be written to files in
\configClass{parametrization:signalBiases}{gnssParametrizationType:signalBiases}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Code biases.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationCodeBiases : public GnssParametrizationBase
{
  class Parameter
  {
  public:
    GnssTransmitterPtr trans;
    GnssReceiverPtr    recv;
    GnssParameterIndex index;
    Matrix             Bias;
  };

  Gnss                      *gnss;
  std::string                name, nameConstraint;
  GnssTransceiverSelectorPtr selectTransmitters, selectReceivers;
  Bool                       isLinearBias;
  Bool                       applyConstraint;
  Double                     sigmaZeroMean;
  std::vector<Parameter*>    paraTrans, paraRecv;
  std::vector<UInt>          idxBiasTrans, idxBiasRecv; // indices in zeroMean matrix
  Matrix                     zeroMeanDesign;            // zero mean observation equations

public:
  GnssParametrizationCodeBiases(Config &config);
 ~GnssParametrizationCodeBiases();

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  void   constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
};

/***********************************************/

#endif
