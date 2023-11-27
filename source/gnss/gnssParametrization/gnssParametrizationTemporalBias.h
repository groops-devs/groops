/***********************************************/
/**
* @file gnssParametrizationTemporalBias.h
*
* @brief Temporal changing signal bias.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTEMPORALBIAS__
#define __GROOPS_GNSSPARAMETRIZATIONTEMPORALBIAS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationTemporalBias = R"(
\subsection{TemporalBias}\label{gnssParametrizationType:temporalBias}
This parametrization resolves the issue of some phase observations suffering from time-variable biases.
Such a phenomenon has been found to affect GPS block IIF satellites on the L5 phase measurements
(see Montenbruck et al. 2011, DOI: \href{https://doi.org/10.1007/s10291-011-0232-x}{10.1007/s10291-011-0232-x}).

For these time-variable biases an appropriate temporal representation has to be defined in
\configClass{parametrizationTemporal}{parametrizationTemporalType}.
For example, time-variable biases for GPS block IIF L5 phase observations (\configClass{type}{gnssType}=\verb|L5*G|)
can be represented by a cubic spline with a nodal distance of one hour.

The result is written as a \file{times series file}{instrument} at the processing sampling
or the sampling set by \configClass{GnssProcessing:processingStep:selectEpochs}{gnssProcessingStepType:selectEpochs}).

This parametrization should be set up in addition to the constant
\configClass{parametrization:signalBiases}{gnssParametrizationType:signalBiases}.
Depending on the temporal representation a temporal zero-mean constraint is needed
to separate this parametrization from the constant component. The constraint equations are added with
a standard deviation of \config{sigmaZeroMeanConstraint}.

The \file{parameter names}{parameterName} are
\verb|<prn>:signalBias.<gnssType>:<temporal>:<interval>|.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Temporal changing signal bias.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationTemporalBias : public GnssParametrizationBase
{
  class Parameter
  {
  public:
    GnssTransmitterPtr trans;
    GnssParameterIndex index;
    Vector             x;
    Vector             bias;
  };

  Gnss                      *gnss;
  std::string                name, nameConstraint;
  PlatformSelectorPtr        selectTransmitters;
  FileName                   fileNameOut, fileNameIn;
  GnssType                   type;
  ParametrizationTemporalPtr temporal;
  Bool                       applyConstraint;
  Double                     sigmaZeroMean;
  std::vector<Parameter*>    parameters;

public:
  GnssParametrizationTemporalBias(Config &config);
 ~GnssParametrizationTemporalBias();

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   observationCorrections(GnssObservationEquation &eqn) const override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  void   constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
