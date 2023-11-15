/***********************************************/
/**
* @file observationSstVariational.h
*
* @brief Satellite to satellite tracking (Variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2012-06-10
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONSSTVARIATIONAL__
#define __GROOPS_OBSERVATIONSSTVARIATIONAL__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationSstVariational = R"(
\subsection{SstVariational}\label{observationType:sstVariational}
Like \configClass{observation:podVariational}{observationType:podVariational} (see there for details)
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
\item \verb|<satellite1>:arc<no>.<parametrizationAcceleration>:*:*|,
\item \verb|<satellite1>:arc<no>.position0.x::|,
\item \verb|<satellite1>:arc<no>.position0.y::|,
\item \verb|<satellite1>:arc<no>.position0.z::|.
\item \verb|<satellite1>:arc<no>.velocity0.x::|,
\item \verb|<satellite1>:arc<no>.velocity0.y::|,
\item \verb|<satellite1>:arc<no>.velocity0.z::|.
\item \verb|<satellite2>:<parametrizationAcceleration>:*:*|,
\item \verb|<satellite2>:arc<no>.<parametrizationAcceleration>:*:*|,
\item \verb|<satellite2>:arc<no>.position0.x::|,
\item \verb|<satellite2>:arc<no>.position0.y::|,
\item \verb|<satellite2>:arc<no>.position0.z::|.
\item \verb|<satellite2>:arc<no>.velocity0.x::|,
\item \verb|<satellite2>:arc<no>.velocity0.y::|,
\item \verb|<satellite2>:arc<no>.velocity0.z::|.
\item \verb|<satellite1>.<satellite2>:<parametrizationSatelliteTracking>:*:*|.
\end{itemize}
)";
#endif


/***********************************************/

#include "misc/observation/observationMiscSstVariational.h"
#include "misc/observation/covariancePod.h"
#include "misc/observation/covarianceSst.h"

/***** CLASS ***********************************/

/** @brief Satellite to satellite tracking (Variational equations).
* @ingroup observationGroup
* @see Observation */
class ObservationSstVariational : public Observation
{
  ObservationMiscSstVariationalPtr observationMisc;
  CovarianceSstPtr covSst;
  CovariancePodPtr covPod1, covPod2;

public:
  ObservationSstVariational(Config &config);
 ~ObservationSstVariational() {}

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
