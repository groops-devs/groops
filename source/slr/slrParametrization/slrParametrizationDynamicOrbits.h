/***********************************************/
/**
* @file slrParametrizationDynamicOrbits.h
*
* @brief Orbits by variational equations.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONDYNAMICORBITS__
#define __GROOPS_SLRPARAMETRIZATIONDYNAMICORBITS__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationDynamicOrbits = R"(
\subsection{DynamicOrbits}\label{slrParametrizationType:dynamicOrbits}
The estimation of (reduced) dynamic orbits is formulated as variational equations.
It is based on \configFile{inputfileVariational}{variationalEquation} calculated with \program{PreprocessingVariationalEquation}.
Necessary integrations are performed by integrating a moving interpolation polynomial of degree \config{integrationDegree}.
The \configClass{parametrizationAcceleration}{parametrizationAccelerationType} must include at least those
parameters that were estimated in \program{PreprocessingVariationalEquationOrbitFit}.
Additional \configClass{stochasticPulse}{timeSeriesType} parameters can be set up to reduce orbit mismodeling.

The parameters and \file{parameter names}{parameterName} are divided into global
\begin{itemize}
\item \verb|<satellite>:<parametrizationAcceleration>:*:*|,
\item \verb|<satellite>:stochasticPulse.x::<time>|,
\item \verb|<satellite>:stochasticPulse.y::<time>|,
\item \verb|<satellite>:stochasticPulse.z::<time>|,
\end{itemize}
and arc related parameters
\begin{itemize}
\item \verb|<satellite>:arc<no>.<parametrizationAcceleration>:*:*|,
\item \verb|<satellite>:arc<no>.position0.x::|,
\item \verb|<satellite>:arc<no>.position0.y::|,
\item \verb|<satellite>:arc<no>.position0.z::|.
\item \verb|<satellite>:arc<no>.velocity0.x::|,
\item \verb|<satellite>:arc<no>.velocity0.y::|,
\item \verb|<satellite>:arc<no>.velocity0.z::|.
\end{itemize}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slr.h"
#include "slr/slrParametrization/slrParametrization.h"
#include "slr/slrParametrization/slrParametrizationGravityField.h"

/***** CLASS ***********************************/

/** @brief Orbits by variational equations.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationDynamicOrbits : public SlrParametrizationBase
{
  class Parameter
  {
  public:
    SlrSatellitePtr            sat;
    SlrParameterIndex          index;
    UInt                       idxParameterSatellite;
    std::vector<ParameterName> parameterNames;
    std::vector<Time>          times;
    Matrix                     PosDesign, VelDesign;
    Vector                     pos, vel;
    Vector                     x;
    Polynomial                 polynomial;
  };

  Slr                           *slr;
  std::string                    name;
  PlatformSelectorPtr            selectorSatellites;
  std::vector<const SlrParametrizationGravityField*> paramGravityField;
  FileName                       fileNameOrbit, fileNameParameter, fileNameVariational;
  std::vector<Time>              pulses;
  ParametrizationAccelerationPtr parametrizationAcceleration;
  EphemeridesPtr                 ephemerides;
  UInt                           integrationDegree, interpolationDegree;
  std::vector<Parameter*>        parameters;

public:
  SlrParametrizationDynamicOrbits(Config &config);
 ~SlrParametrizationDynamicOrbits();

  void   init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void   initParameter(SlrNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix    (const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const override;
  Double updateParameter (const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults    (const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
