/***********************************************/
/**
* @file observationPodIntegral.cpp
*
* @brief Precise Orbit Data (POD) observations (short arc integral).
* Solution of the Fredholm integral.
*
* @author Torsten Mayer-Guerr
* @date 2003-12-22
*
*/
/***********************************************/

#include "base/import.h"
#include "misc/observation/observationMiscPodIntegral.h"
#include "misc/observation/covariancePod.h"
#include "classes/observation/observation.h"
#include "classes/observation/observationPodIntegral.h"

/***********************************************/

ObservationPodIntegral::ObservationPodIntegral(Config &config)
{
  try
  {
    observationMisc = ObservationMiscPodIntegralPtr(new ObservationMiscPodIntegral(config));
    readConfig(config, "covariancePod", covPod, Config::OPTIONAL, "", "covariance matrix of kinematic orbits");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationPodIntegral::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
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
