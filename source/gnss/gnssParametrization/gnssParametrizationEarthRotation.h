/***********************************************/
/**
* @file gnssParametrizationEarthRotation.h
*
* @brief EarthRotation.
* @see GnssParametrization
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONEARTHROTATION__
#define __GROOPS_GNSSPARAMETRIZATIONEARTHROTATION__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationEarthRotation = R"(
\subsection{EarthRotation}\label{gnssParametrizationType:earthRotation}
Earth rotation parameters (ERPs) can be estimated by defining
\config{estimatePole} ($x_p$, $y_p\, [mas]$) and \config{estimateUT1} ($dUT1\, [ms], LOD$).

Estimating length of day (LOD) with the sign according to IGS conventions requires a negative
value in \configClass{parametrizationTemporal:trend:timeStep}{parametrizationTemporalType:trend}.

Constraints on the defined parameters can be added via
\configClass{parametrization:constraints}{gnssParametrizationType:constraints}.
An example would be to set up \configClass{estimateUT1:constant}{parametrizationTemporalType:constant}
so the $dUT1$ parameter is included in the normal equation system . Since $dUT1$ cannot be
determined by GNSS, a hard constraint to its a priori value can then be added.

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
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief EarthRotation.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationEarthRotation : public GnssParametrizationBase
{
  Gnss                      *gnss;
  std::string                name;
  FileName                   fileNameEOP;
  ParametrizationTemporalPtr parametrizationPole, parametrizationUT1, parametrizationNutation;
  GnssParameterIndex         indexParameterPole,  indexParameterUT1, indexParameterNutation;
  Vector                     xPole, xUT1, xNutation;

public:
  GnssParametrizationEarthRotation(Config &config);

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
