/***********************************************/
/**
* @file gnssParametrizationTroposphere.h
*
* @brief Troposphere.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTROPOSPHERE__
#define __GROOPS_GNSSPARAMETRIZATIONTROPOSPHERE__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationTroposphere = R"(
\subsection{Troposphere}\label{gnssParametrizationType:troposphere}
A priori tropospheric correction is handled by a \configClass{troposphere}{troposphereType} model (e.g. Vienna Mapping Functions 3).
Additional parameters for zenith wet delay and gradients can be set up via
\configClass{troposphereWetEstimation}{parametrizationTemporalType} (usually 2-hourly linear splines)
and \configClass{troposphereGradientEstimation}{parametrizationTemporalType} (usually a daily trend).
These parameters can be soft-constrained using
\configClass{parametrization:constraints}{gnssParametrizationType:constraints}
to avoid an unsolvable system of normal equations in case of data gaps.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "classes/troposphere/troposphere.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Troposphere.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationTroposphere : public GnssParametrizationBase
{
  class Parameter
  {
  public:
    UInt                idRecv, idTropo;
    GnssParameterIndex  indexWet, indexGradient;
    Vector              xWet, xGradient;
    std::vector<Double> zenitDelayWet, gradientX, gradientY;
  };

  Gnss                      *gnss;
  std::string                name;
  PlatformSelectorPtr        selectReceivers;
  FileName                   fileNameTropo;
  TropospherePtr             troposphere;
  ParametrizationTemporalPtr parametrizationWet;
  ParametrizationTemporalPtr parametrizationGradient;
  std::vector<Parameter*>    parameters;


public:
  GnssParametrizationTroposphere(Config &config);
 ~GnssParametrizationTroposphere();

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   observationCorrections(GnssObservationEquation &eqn) const override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
