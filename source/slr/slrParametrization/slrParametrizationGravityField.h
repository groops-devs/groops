/***********************************************/
/**
* @file slrParametrizationGravityField.h
*
* @brief GravityField.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONGRAVITYFIELD__
#define __GROOPS_SLRPARAMETRIZATIONGRAVITYFIELD__

// Latex documentation

#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationGravityField = R"(
\subsection{GravityField}\label{slrParametrizationType:gravityField}
Estimates a (time depending) gravity field together with at least one
\configClass{parametrization:dynamicOrbits}{slrParametrizationType:dynamicOrbits}.
The parametrization of the gravity field can be set with
\configClass{parametrization}{parametrizationGravityType}.

The \file{parameter names}{parameterName} are \verb|gravityfield:<parametrization>:*:*|.
)";
#endif

/***********************************************/

#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief GravityField.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationGravityField : public SlrParametrizationBase
{
  Slr                      *slr;
  std::string               name;

public:
  ParametrizationGravityPtr parametrization;
  SlrParameterIndex         indexParameter;
  Vector                    x;

  SlrParametrizationGravityField(Config &config);

  void   getParametrizationGravity(std::vector<const SlrParametrizationGravityField*> &paramGravityField) const override;
  void   init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void   initParameter(SlrNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  Double updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
