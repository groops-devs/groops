/***********************************************/
/**
* @file observationPodAcceleration.h
*
* @brief Precise Orbit Data (POD) observations (Acceleration approach).
*
* @author Torsten Mayer-Guerr
* @date 2006-02-02
* update 2011-06-09 Norbert Zehentner
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONPODACCELERATION__
#define __GROOPS_OBSERVATIONPODACCELERATION__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationPodAcceleration = R"(
\subsection{PodAcceleration}\label{observationType:podAcceleration}
The observation equations for precise orbit data (POD) are given by
\begin{equation}
\ddot{\M r}(t) - \M g_0(t) = \nabla V(\M r, t),
\end{equation}
where the accelerations of the satellite $\ddot{\M r}(t)$ are derived from the kinematic positions
in \configClass{rightHandSide}{podRightSideType}. The orbit differentation is performed by a moving
polynomial interpolation or approximation with degree \config{interpolationDegree}
and number of used epochs \config{numberOfEpochs}. The reference forces $\M g_0(t)$ are computed
with the background models in \configClass{rightHandSide}{podRightSideType}.

All instrument data \configFile{inputfileOrbit}{instrument}, \configFile{inputfileStarCamera}{instrument},
and \configFile{inputfileAccelerometer}{instrument} must be synchronous and be given
with a constant sampling without any gaps in each short arc (see \program{InstrumentSynchronize}).

The unknown gravity field $\nabla V(\M r, t)$ parametrized by \configClass{parametrizationGravity}{parametrizationGravityType}
is not evaluated at the observed positions but at the orbit given by \configFile{inputfileOrbit}{instrument}.
The same is true for the reference forces. This orbit may be a more accurate dynamical orbit but
in most cases the kinematic orbit provides good results.

The accuracy or the full covariance matrix of the precise orbit data is provided in
\configClass{covariancePod}{covariancePodType} and can be estimated with \program{PreprocessingPod}.

The following parameters with \file{parameter names}{parameterName} are set up:
\begin{itemize}
\item \verb|*:<parametrizationGravity>:*:*|,
\item \verb|<satellite>:<parametrizationAcceleration>:*:*|.
\end{itemize}
)";
#endif

/***********************************************/

#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/observation/observation.h"
#include "misc/observation/observationMisc.h"
#include "misc/observation/covariancePod.h"

/***** CLASS ***********************************/

/** @brief Precise Orbit Data (POD) observations (Acceleration approach).
* @ingroup observationGroup
* @see Observation */
class ObservationPodAcceleration : public Observation
{
  SatelliteModelPtr              satellite;
  std::vector<PodRightSidePtr>   rhs; // right hand sides
  CovariancePodPtr               covPod;
  InstrumentFile                 orbitFile;
  InstrumentFile                 starCameraFile;
  EarthRotationPtr               earthRotation;
  EphemeridesPtr                 ephemerides;
  ParametrizationGravityPtr      parametrization;
  ParametrizationAccelerationPtr parametrizationAcceleration;
  Vector                         coeff;

public:
  ObservationPodAcceleration(Config &config);
 ~ObservationPodAcceleration() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override;
  UInt parameterCount()        const override;
  UInt gravityParameterCount() const override {return parametrization->parameterCount();}
  UInt rightSideCount()        const override {return rhs.size();}
  UInt arcCount()              const override {return orbitFile.arcCount();}
  void parameterName(std::vector<ParameterName> &name) const override;

  void observation(UInt arc, Matrix &l, Matrix &A, Matrix &B) override;
};

/***********************************************/

#endif /* __GROOPS_OBSERVATION__ */
