/***********************************************/
/**
* @file observationPodVariational.h
*
* @brief Precise Orbit Data (POD) observations (Variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2012-08-20
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONPODVARIATIONAL__
#define __GROOPS_OBSERVATIONPODVARIATIONAL__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationPodVariational = R"(
\subsection{PodVariational}\label{observationType:podVariational}
The observation equations for precise orbit data (POD) are formulated as variational equations.
It is based on \file{inputfileVariational}{variationalEquation} calculated with \program{PreprocessingVariationalEquation}.
Necessary integrations are performed by integrating a moving interpolation polynomial of degree \config{integrationDegree}.

The kinematic positions as pseudo observations are taken from
\config{rightHandSide} and should not given equally spaced in time. The observation
equations are interpolated to these times by a moving polynomial of degree \config{interpolationDegree}.

The accuracy or the full covariance matrix of the precise orbit data is provided in
\configClass{covariancePod}{covariancePodType} and can be estimated with \program{PreprocessingPod}.

\config{accelerateComputation}: In the event that the sampling of the kinematic orbit is much higher than the sampling
of the variational equations (e.g. 1 second vs. 5 seconds) the accumulation of the observation equations
can be accelerated by transforming the observation equations
\begin{equation}
  \M l = \M J \M A \M x + \M e,
\end{equation}
where $\M J$ describes the interpolation of the sampling of the variational design matrix~$\M A$
to the sampling of the observations $\M l$ with more rows than columns. The QR decomposition
\begin{equation}
  \M J = \begin{pmatrix} \M Q_1 & \M Q_2 \end{pmatrix}
         \begin{pmatrix} \M R \\ \M 0 \end{pmatrix}.
\end{equation}
can be used to transform the observation equations
\begin{equation}
  \begin{pmatrix} \M Q_1^T \M l \\ \M Q_2^T \M l \end{pmatrix} =
  \begin{pmatrix} \M Q_1^T \M R \\ \M 0 \end{pmatrix} \M A \M x +
  \begin{pmatrix} \M Q_1^T \M e \\ \M Q_2^T \M e \end{pmatrix}.
\end{equation}
As the zero lines should not be considered the computational time for the accumulation is reduced.
This option is not meaningful for evaluating the residuals such in \program{PreprocessingPod}.

The following parameters with \file{parameter names}{parameterName} are set up:
\begin{itemize}
\item \verb|*:<parametrizationGravity>:*:*|,
\item \verb|<satellite>:<parametrizationAcceleration>:*:*|,
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

#include "classes/observation/observation.h"
#include "misc/observation/observationMiscPodVariational.h"
#include "misc/observation/covariancePod.h"

/***** CLASS ***********************************/

/** @brief Precise Orbit Data (POD) observations (Variational equations).
* @ingroup observationGroup
* @see Observation */
class ObservationPodVariational : public Observation
{
  ObservationMiscPodVariationalPtr observationMisc;
  CovariancePodPtr                 covPod;

public:
  ObservationPodVariational(Config &config);
 ~ObservationPodVariational() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override {return observationMisc->setInterval(timeStart, timeEnd);}
  UInt parameterCount()          const override {return observationMisc->parameterCount();}
  UInt gravityParameterCount()   const override {return observationMisc->gravityParameterCount();}
  UInt rightSideCount()          const override {return observationMisc->rightSideCount();}
  UInt arcCount()                const override {return observationMisc->arcCount();}
  void parameterName(std::vector<ParameterName> &name) const override {observationMisc->parameterName(name);}

  void observation(UInt arc, Matrix &l, Matrix &A, Matrix &B) override;
};

/***********************************************/

#endif /* __GROOPS_OBSERVATION__ */
