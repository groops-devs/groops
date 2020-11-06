/***********************************************/
/**
* @file observationGradiometer.h
*
* @brief GOCE gradiometer observations.
*
* @author Torsten Mayer-Guerr
* @date 2011-05-14
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONGRADIOMETER__
#define __GROOPS_OBSERVATIONGRADIOMETER__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationGradiometer = R"(
\subsection{Gradiometer}\label{observationType:gradiometer}
Observation equations for satellite gravity gradiometry (SGG)
\begin{equation}
  \nabla\nabla V(\M r) =
  \begin{pmatrix}
    \frac{\partial^2 V}{\partial x^2}         & \frac{\partial^2 V}{\partial x\partial y} & \frac{\partial^2 V}{\partial x\partial z} \\
    \frac{\partial^2 V}{\partial y\partial x} & \frac{\partial^2 V}{\partial y^2}         & \frac{\partial^2 V}{\partial y\partial z} \\
    \frac{\partial^2 V}{\partial z\partial x} & \frac{\partial^2 V}{\partial z\partial y} & \frac{\partial^2 V}{\partial z^2}
  \end{pmatrix}.
\end{equation}
From the \configFile{inputfileGradiometer}{instrument} observations precomputed \configFile{inputfileReferenceGradiometer}{instrument}
together with other background models are reduced, all given in \configClass{rightHandSide}{sggRightSideType}.

All instrument data \configFile{inputfileGradiometer}{instrument}, \configFile{inputfileOrbit}{instrument},
and \configFile{inputfileStarCamera}{instrument} must be synchronous and be diveded
into each short arcs (see \program{InstrumentSynchronize}).

Additional to the \configClass{parametrizationGravity}{parametrizationGravityType}
an (temporal changing) bias for each gradiometer component and arc can be estimated with
\configClass{parametrizationBias}{parametrizationTemporalType}.

The accuracy or the full covariance matrix of the gradiometer is provided in
\config{covarianceSgg} and can be estimated with \program{PreprocessingGradiometer}.

)";
#endif

/***********************************************/

#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/tides/tides.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/observation/observation.h"
#include "misc/observation/observationMisc.h"

/***** CLASS ***********************************/

/** @brief GOCE gradiometer observations.
* @ingroup observationGroup
* @see Observation */
class ObservationGradiometer : public Observation
{
  std::vector<SggRightSidePtr> rhs; // right hand sides
  InstrumentFile               orbitFile;
  InstrumentFile               starCameraFile;
  EarthRotationPtr             earthRotation;
  EphemeridesPtr               ephemerides;
  ParametrizationGravityPtr    parametrization;
  ParametrizationTemporalPtr   sggBias;
  UInt                         componentCount;
  Bool                         useXX, useYY, useZZ;
  Bool                         useXY, useXZ, useYZ;
  Double                       sigma;
  Vector                       sigmaArc;
  Matrix                       CovCholesky; // cholesky decomposition of the Covariance

public:
  ObservationGradiometer(Config &config);
 ~ObservationGradiometer() {}

  UInt parameterCount()          const {return parametrization->parameterCount();}
  UInt gravityParameterCount()   const {return parametrization->parameterCount();}
  UInt rightSideCount()          const {return rhs.size();}
  UInt arcCount()                const {return orbitFile.arcCount();}
  void parameterName(std::vector<ParameterName> &name) const {parametrization->parameterName(name);}

  void observation(UInt arc, Matrix &l, Matrix &A, Matrix &B);
};

/***********************************************/

#endif
