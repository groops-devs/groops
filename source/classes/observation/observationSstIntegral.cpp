/***********************************************/
/**
* @file observationSstIntegral.cpp
*
* @brief Satellite to satellite tracking (Short Arc Integral).
*
* @author Torsten Mayer-Guerr
* @date 2009-11-01
*
*/
/***********************************************/

#include "base/import.h"
#include "misc/observation/observationMiscSstIntegral.h"
#include "misc/observation/covarianceSst.h"
#include "misc/observation/covariancePod.h"
#include "classes/observation/observation.h"
#include "observationSstIntegral.h"

/***********************************************/

ObservationSstIntegral::ObservationSstIntegral(Config &config)
{
  try
  {
    observationMisc = ObservationMiscSstIntegralPtr(new ObservationMiscSstIntegral(config));
    readConfig(config, "covarianceSst",  covSst,  Config::MUSTSET, "", "covariance matrix of satellite to satellite tracking observations");
    readConfig(config, "covariancePod1", covPod1, Config::MUSTSET, "", "covariance matrix of kinematic orbits (satellite 1)");
    readConfig(config, "covariancePod2", covPod2, Config::MUSTSET, "", "covariance matrix of kinematic orbits (satellite 2)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationSstIntegral::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    ObservationMiscSstIntegral::Arc arc = observationMisc->computeArc(arcNo, covSst, covPod1, covPod2);
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
