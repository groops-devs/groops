/***********************************************/
/**
* @file observationDualSstVariational.h
*
* @brief Satellite to satellite tracking with two simultaneous ranging observations (Variational equations).
*
* @author Andreas Kvas
* @date 2020-07-24
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONDUALSSTVARIATIONAL__
#define __GROOPS_OBSERVATIONDUALSSTVARIATIONAL__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationDualSstVariational = R"(
\subsection{DualSstVariational}\label{observationType:dualSstVariational}
Like \configClass{observation:sstVariational}{observationType:sstVariational} (see there for details)
but with two simultaneous satellite-to-satellite (SST) observations.

This class reads two SST observation files (\configFile{inputfileSatelliteTracking1}{instrument} and
\configFile{inputfileSatelliteTracking2}{instrument}).
Empirical parameters for the SST observations can be setup independently for both SST observation
types with \configClass{parametrizationSst1}{parametrizationSatelliteTrackingType} and
\configClass{parametrizationSst2}{parametrizationSatelliteTrackingType}.

Both SST observation types are reduced by the same background models and the same impact
of accelerometer measurements. The covariance matrix of the reduced observations should not consider
the the instrument noise only (\configClass{covarianceSst1/2}{covarianceSstType}) but must
take the cross correlations \configClass{covarianceAcc}{covarianceSstType} into account.
The covariance matrix of the reduced observations is given by
\begin{equation}
  \M\Sigma(\begin{bmatrix} \Delta l_{SST1} \\ \Delta l_{SST2} \end{bmatrix})
  = \begin{bmatrix} \M\Sigma_{SST1} + \M\Sigma_{ACC} & \M\Sigma_{ACC} \\
                   \M\Sigma_{ACC} & \M\Sigma_{SST2} + \M\Sigma_{ACC}
    \end{bmatrix}.
\end{equation}
)";
#endif

/***********************************************/

#include "misc/observation/observationMiscDualSstVariational.h"
#include "misc/observation/covariancePod.h"
#include "misc/observation/covarianceSst.h"

/***** CLASS ***********************************/

/** @brief Satellite to satellite tracking with two simultaneous ranging observations (Variational equations).
* @ingroup observationGroup
* @see Observation */
class ObservationDualSstVariational : public Observation
{
  ObservationMiscDualSstVariationalPtr observationMisc;
  CovarianceSstPtr covSst1, covSst2, covAcc;
  CovariancePodPtr covPod1, covPod2;

public:
  ObservationDualSstVariational(Config &config);
 ~ObservationDualSstVariational() {}

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
