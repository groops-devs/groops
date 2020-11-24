/***********************************************/
/**
* @file observationPodIntegral.h
*
* @brief Precise Orbit Data (POD) observations (short arc integral).
* Solution of the Fredholm integral.
*
* @author Torsten Mayer-Guerr
* @date 2003-1222
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONPODINTEGRAL__
#define __GROOPS_OBSERVATIONPODINTEGRAL__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationPodIntegral = R"(
\subsection{PodIntegral}\label{observationType:podIntegral}
The observation equations for precise orbit data (POD) of short arcs are given by
\begin{equation}
  {\M r}_\epsilon(\tau) = {\M r}_A(1-\tau) + {\M r}_B\tau - T^2\int_0^1 K(\tau,\tau')
  \left(\M f_0(\tau')+\nabla V(\tau')\right)\,d\tau'
\end{equation}
with the integral kernel
\begin{equation}
  K(\tau,\tau') = \begin{cases} \tau'(1-\tau) & \text{for }\tau'\le\tau \\
  \tau(1-\tau') & \text{for }\tau'>\tau \end{cases},
\end{equation}
and the normalized time
\begin{equation}
  \tau = \frac{t-t_A}{T}\qquad\text{with}\qquad T=t_B-t_A.
\end{equation}
The kinematic positions~${\M r}_\epsilon(\tau)$ as pseudo observations are taken from
\configClass{rightHandSide}{podRightSideType}. From these positions the influence of the reference forces $\M f_0(\tau)$
is subtracted which are computed with the background models in \configClass{rightHandSide}{podRightSideType}.
The integral is solved by the integration of a moving interpolation polynomial of degree \config{integrationDegree}.
The boundary values ${\M r}_A$ and ${\M r}_B$ (satellite's state vector) are estimated per arc
and are usally directly eliminated if \config{keepSatelliteStates} is not set.

The unknown gravity field $\nabla V(\M r, t)$ parametrized by \configClass{parametrizationGravity}{parametrizationGravityType}
is not evaluated at the observed positions but at the orbit given by \configFile{inputfileOrbit}{instrument}.
The same is true for the reference forces. The linearized effect of the gravity field change by the position
adjustment is taken into account by \config{gradientfield}. This may be a low order field up to a
spherical harmonics degree of $n=2$ or $n=3$.

The \configFile{inputfileOrbit}{instrument}, \configFile{inputfileStarCamera}{instrument}, and \configFile{inputfileAccelerometer}{instrument}
must be synchronous and must be given with a constant sampling and without any gaps in each short arc
(see \program{InstrumentSynchronize}).
The kinematic positions~${\M r}_\epsilon(\tau)$ should not given equally spaced in time
but must be divided into the same arcs as the other instrument data.
The observation equations are interpolated to this time by a polynomial interpolation
with degree \config{interpolationDegree}.

The accuracy or the full covariance matrix of the precise orbit data is provided in
\configClass{covariancePod}{covariancePodType} and can be estimated with \program{PreprocessingPod}.

For \config{accelerateComputation} see \configClass{observation:podVariational}{observationType:podVariational}.
)";
#endif

/***********************************************/

#include "classes/observation/observation.h"
#include "misc/observation/observationMiscPodIntegral.h"
#include "misc/observation/covariancePod.h"

/***** CLASS ***********************************/

/** @brief Precise Orbit Data (POD) observations (short arc integral).
* @ingroup observationGroup
* Solution of the Fredholm integral.
* @see Observation */
class ObservationPodIntegral : public Observation
{
  ObservationMiscPodIntegralPtr observationMisc;
  CovariancePodPtr              covPod;

public:
  ObservationPodIntegral(Config &config);
 ~ObservationPodIntegral() {}

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
