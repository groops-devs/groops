/***********************************************/
/**
* @file slrParametrizationTroposphere.h
*
* @brief Troposphere.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONTROPOSPHERE__
#define __GROOPS_SLRPARAMETRIZATIONTROPOSPHERE__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationTroposphere = R"(
\subsection{Troposphere}\label{slrParametrizationType:troposphere}
A priori tropospheric correction is handled by a \configClass{troposphere}{troposphereType} model (e.g. Mendes and Pavlis).
Additional parameters in $[m]$ for zenith delay can be set up via
\configClass{troposphereEstimation}{parametrizationTemporalType}.
These parameters can be soft-constrained using
\configClass{parametrization:constraints}{slrParametrizationType:constraints}
to avoid an unsolvable system of normal equations in case of data gaps.

The \file{parameter names}{parameterName} are \verb|<station>:troposphere:<temporal>:<interval>|.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/troposphere/troposphere.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slr.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Troposphere.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationTroposphere : public SlrParametrizationBase
{
  class Parameter
  {
  public:
    UInt                idStat, idTropo;
    SlrParameterIndex   index;
    Vector              x;
    std::vector<Double> zenitDelay;
  };

  Slr                       *slr;
  std::string                name;
  PlatformSelectorPtr        selectorStations;
  FileName                   fileNameTropo;
  TropospherePtr             troposphere;
  ParametrizationTemporalPtr parametrization;
  std::vector<Parameter*>    parameters;


public:
  SlrParametrizationTroposphere(Config &config);
 ~SlrParametrizationTroposphere();

  void   init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void   observationCorrections(SlrObservationEquation &eqn) const override;
  void   initParameter(SlrNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const override;
  Double updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
