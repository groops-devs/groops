/***********************************************/
/**
* @file observationSstIntegral.h
*
* @brief Satellite to satellite tracking (Short Arc Integral).
*
* @author Torsten Mayer-Guerr
* @date 2009-11-01
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONSSTINTEGRAL__
#define __GROOPS_OBSERVATIONSSTINTEGRAL__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationSstIntegral = R"(
\subsection{SstIntegral}\label{observationType:sstIntegral}
Like \configClass{observation:podIntegral}{observationType:podIntegral} (see there for details)
but with two satellites and additional satellite-to-satellite (SST) observations.

If multiple \configFile{inputfileSatelliteTracking}{instrument} are given
all data are add together. So corrections in extra files like the light time correction
can easily be added. Empirical parameters for the SST observations can be setup with
\configClass{parametrizationSst}{parametrizationSatelliteTrackingType}.
The accuracy or the full covariance matrix of SST is provided in
\configClass{covarianceSst}{covarianceSstType}.

The following parameters with \file{parameter names}{parameterName} are set up:
\begin{itemize}
\item \verb|*:<parametrizationGravity>:*:*|,
\item \verb|<satellite1>:<parametrizationAcceleration>:*:*|,
\item \verb|<satellite2>:<parametrizationAcceleration>:*:*|,
\item \verb|<satellite1>.<satellite2>:<parametrizationSatelliteTracking>:*:*|,
\end{itemize}
and for each arc if \config{keepSatelliteStates} is set
\begin{itemize}
\item \verb|<satellite1>:arc<no>.position.start.x::|,
\item \verb|<satellite1>:arc<no>.position.start.y::|,
\item \verb|<satellite1>:arc<no>.position.start.z::|.
\item \verb|<satellite1>:arc<no>.position.end.x::|,
\item \verb|<satellite1>:arc<no>.position.end.y::|,
\item \verb|<satellite1>:arc<no>.position.end.z::|.
\item \verb|<satellite2>:arc<no>.position.start.x::|,
\item \verb|<satellite2>:arc<no>.position.start.y::|,
\item \verb|<satellite2>:arc<no>.position.start.z::|.
\item \verb|<satellite2>:arc<no>.position.end.x::|,
\item \verb|<satellite2>:arc<no>.position.end.y::|,
\item \verb|<satellite2>:arc<no>.position.end.z::|.
\end{itemize}
)";
#endif

/***********************************************/

#include "misc/observation/observationMiscSstIntegral.h"
#include "misc/observation/covariancePod.h"
#include "misc/observation/covarianceSst.h"

/***** CLASS ***********************************/

/** @brief Satellite to satellite tracking (Short Arc Integral).
* @ingroup observationGroup
* @see Observation */
class ObservationSstIntegral : public Observation
{
  ObservationMiscSstIntegralPtr observationMisc;
  CovarianceSstPtr covSst;
  CovariancePodPtr covPod1, covPod2;

public:
  ObservationSstIntegral(Config &config);
 ~ObservationSstIntegral() {}

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
