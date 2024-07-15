/***********************************************/
/**
* @file slrParametrizationEarthRotation.h
*
* @brief EarthRotation.
* @see SlrParametrization
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONEARTHROTATION__
#define __GROOPS_SLRPARAMETRIZATIONEARTHROTATION__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationEarthRotation = R"(
\subsection{EarthRotation}\label{slrParametrizationType:earthRotation}
Earth rotation parameters (ERPs) can be estimated by defining
\config{estimatePole} ($x_p$, $y_p$) and \config{estimateUT1} ($dUT1, LOD$).

Estimating length of day (LOD) with the sign according to IGS conventions requires a negative
value in \configClass{parametrizationTemporal:trend:timeStep}{parametrizationTemporalType:trend}.

Constraints on the defined parameters can be added via
\configClass{parametrization:constraints}{slrParametrizationType:constraints}.
An example would be to set up \configClass{estimateUT1:constant}{parametrizationTemporalType:constant}
so the $dUT1$ parameter is included in the normal equation system . Since $dUT1$ cannot be
determined by SLR, a hard constraint to its a priori value can then be added.

The \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|earth:polarMotion.xp:<temporal>:<interval>|,
\item \verb|earth:polarMotion.yp:<temporal>:<interval>|,
\item \verb|earth:UT1:<temporal>:<interval>|,
\item \verb|earth:nutation.X:<temporal>:<interval>|,
\item \verb|earth:nutation.>:<temporal>:<interval>|.
\end{itemize}
)";
#endif

/***********************************************/

#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief EarthRotation.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationEarthRotation : public SlrParametrizationBase
{
  Slr                       *slr;
  std::string                name;
  FileName                   fileNameEOP;
  ParametrizationTemporalPtr parametrizationPole, parametrizationUT1, parametrizationNutation;
  SlrParameterIndex          indexParameterPole,  indexParameterUT1, indexParameterNutation;
  Vector                     xPole, xUT1, xNutation;

public:
  SlrParametrizationEarthRotation(Config &config);

  void   init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void   initParameter(SlrNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const override;
  Double updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
