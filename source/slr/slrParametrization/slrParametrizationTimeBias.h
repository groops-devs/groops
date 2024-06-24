/***********************************************/
/**
* @file slrParametrizationTimeBias.h
*
* @brief Time biases.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONTIMEBIAS__
#define __GROOPS_SLRPARAMETRIZATIONTIMEBIAS__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationTimeBias = R"(
\subsection{TimeBias}\label{slrParametrizationType:timeBias}
Estimates a \configClass{temporal changing}{parametrizationTemporalType}
time bias in $[ms]$ for \configClass{selectStations}{platformSelectorType}.

The \file{parameter names}{parameterName} are \verb|<station>:timeBias:<temporal>:<interval>|.
)";
#endif

/***********************************************/

#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Range biases.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationTimeBias : public SlrParametrizationBase
{
  class Parameter
  {
  public:
    SlrStationPtr     station;
    SlrParameterIndex index;
    Vector            x;
  };

  Slr                       *slr;
  std::string                name;
  PlatformSelectorPtr        selectorStations;
  ParametrizationTemporalPtr parametrizationTemporal;
  std::vector<Parameter*>    paraStations;

public:
  SlrParametrizationTimeBias(Config &config);
 ~SlrParametrizationTimeBias();

  void   init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void   initParameter(SlrNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const override;
  Double updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
};

/***********************************************/

#endif
