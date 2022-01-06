/***********************************************/
/**
* @file observationPodEnergy.h
*
* @brief Precise Orbit Data (POD) observations (Energy integral).
*
* @author Torsten Mayer-Guerr
* @date 2006-02-02
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONPODENERGY__
#define __GROOPS_OBSERVATIONPODENERGY__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationPodEnergy = R"(
\subsection{PodEnergy}\label{observationType:podEnergy}
The observation equations for precise orbit data (POD) are given by
\begin{equation}
  \frac{1}{2}\dot{\M r}^2
  -\dot{\M r} \cdot (\M\Omega\times\M r)
  +\int_{t_0}^t(\dot{\M\Omega}\times\M r)\cdot \dot{\M r}\,dt
  - \int_{t_0}^t \M g_0 \cdot\dot{\M r}'\,dt
  = V + E.
\end{equation}
where the velocities of the satellite $\ddot{\M r}(t)$ are derived from
the kinematic positions in \configClass{rightHandSide}{podRightSideType} and the Earth's rotation vector~$\M\Omega(t)$ is modeled
within \configClass{earthRotation}{earthRotationType}. The orbit differentation is
performed by a polynomial interpolation with degree \config{interpolationDegree}.
The integrals are solved a polynomial interpolation with degree \config{integrationDegree}.
The reference forces $\M g_0(t)$ are computed with the background models in \configClass{rightHandSide}{podRightSideType}.

All instrument data \configFile{inputfileOrbit}{instrument}, \configFile{inputfileStarCamera}{instrument}, and \configFile{inputfileAccelerometer}{instrument}
must be synchronous and be given with a constant sampling without any gaps in each short arc
(see \program{InstrumentSynchronize}).

The unknown gravity potential $V(\M r)$ parametrized by \configClass{parametrizationGravity}{parametrizationGravityType}
is not evaluated at the observed positions but at the orbit given by \configFile{inputfileOrbit}{instrument}.
The same is true for the reference forces. This orbit may be a more accurate dynamical orbit but
in most cases the kinematic orbit provides good results.

An unknown energy bias~$E$ per arc is parametrized by \configClass{parametrizationBias}{parametrizationTemporalType}
and should be a constant in theory but temporal changes might help to absorb other unmodelled effects.

The accuracy or the full covariance matrix of the precise orbit data is provided in
\configClass{covariancePod}{covariancePodType} and can be estimated with \program{PreprocessingPod}.
)";
#endif

/***********************************************/

#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/observation/observation.h"
#include "misc/observation/observationMisc.h"
#include "misc/observation/covariancePod.h"

/***** CLASS ***********************************/

/** @brief Precise Orbit Data (POD) observations (Energy integral).
* @ingroup observationGroup
* @see Observation */
class ObservationPodEnergy : public Observation
{
  SatelliteModelPtr            satellite;
  std::vector<PodRightSidePtr> rhs; // right hand sides
  CovariancePodPtr             covPod;
  InstrumentFile               orbitFile;
  InstrumentFile               starCameraFile;
  EarthRotationPtr             earthRotation;
  EphemeridesPtr               ephemerides;
  ParametrizationGravityPtr    parametrization;
  ParametrizationTemporalPtr   bias;
  UInt                         interpolationDegree, integrationDegree;
  Vector                       coeff;
  Matrix                       integrationMatrix;

public:
  ObservationPodEnergy(Config &config);
 ~ObservationPodEnergy() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override {return parametrization->setInterval(timeStart, timeEnd);}
  UInt parameterCount()        const override {return parametrization->parameterCount();}
  UInt gravityParameterCount() const override {return parametrization->parameterCount();}
  UInt rightSideCount()        const override {return rhs.size();}
  UInt arcCount()              const override {return orbitFile.arcCount();}
  void parameterName(std::vector<ParameterName> &name) const override {parametrization->parameterName(name);}

  void observation(UInt arc, Matrix &l, Matrix &A, Matrix &B) override;
};

/***********************************************/

#endif /* __GROOPS_OBSERVATION__ */
