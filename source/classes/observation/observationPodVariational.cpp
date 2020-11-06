/***********************************************/
/**
* @file observationPodVariational.cpp
*
* @brief Precise Orbit Data (POD) observations (Variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2012-08-20
*
*/
/***********************************************/

#include "base/import.h"
#include "misc/observation/observationMiscPodVariational.h"
#include "misc/observation/covariancePod.h"
#include "classes/observation/observation.h"
#include "classes/observation/observationPodVariational.h"

/***********************************************/

ObservationPodVariational::ObservationPodVariational(Config &config)
{
  try
  {
    observationMisc = ObservationMiscPodVariationalPtr(new ObservationMiscPodVariational(config));
    readConfig(config, "covariancePod", covPod, Config::OPTIONAL, "", "covariance matrix of kinematic orbits");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationPodVariational::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    ObservationMiscPod::Arc arc = observationMisc->computeArc(arcNo, covPod);
    l = arc.l;
    A = arc.A;
    B = arc.B;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
