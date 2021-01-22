/***********************************************/
/**
* @file observationDualSstVariational.cpp
*
* @brief Satellite to satellite tracking with two simultaneous ranging observations (Variational equations).
*
* @author Andreas
* @date 2020-07-24
*
*/
/***********************************************/

#include "base/import.h"
#include "misc/observation/observationMiscDualSstVariational.h"
#include "misc/observation/covarianceSst.h"
#include "misc/observation/covariancePod.h"
#include "classes/observation/observation.h"
#include "observationDualSstVariational.h"

/***********************************************/

ObservationDualSstVariational::ObservationDualSstVariational(Config &config)
{
  try
  {
    observationMisc = ObservationMiscDualSstVariationalPtr(new ObservationMiscDualSstVariational(config));
    readConfig(config, "covarianceSst1", covSst1, Config::MUSTSET,  "", "covariance matrix of first satellite to satellite tracking observations");
    readConfig(config, "covarianceSst2", covSst2, Config::MUSTSET,  "", "covariance matrix of second satellite to satellite tracking observations");
    readConfig(config, "covarianceAcc",  covAcc,  Config::OPTIONAL, "", "common covariance matrix of reduced satellite to satellite tracking observations");
    readConfig(config, "covariancePod1", covPod1, Config::MUSTSET,  "", "covariance matrix of kinematic orbits (satellite 1)");
    readConfig(config, "covariancePod2", covPod2, Config::MUSTSET,  "", "covariance matrix of kinematic orbits (satellite 2)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationDualSstVariational::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    ObservationMiscDualSstVariational::Arc arc = observationMisc->computeArc(arcNo, covSst1, covSst2, covAcc, covPod1, covPod2);
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
