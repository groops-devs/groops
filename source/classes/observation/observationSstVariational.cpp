/***********************************************/
/**
* @file observationSstVariational.cpp
*
* @brief Satellite to satellite tracking (Variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2012-06-10
*
*/
/***********************************************/

#include "base/import.h"
#include "misc/observation/observationMiscSstVariational.h"
#include "misc/observation/covarianceSst.h"
#include "misc/observation/covariancePod.h"
#include "classes/observation/observation.h"
#include "observationSstVariational.h"

/***********************************************/

ObservationSstVariational::ObservationSstVariational(Config &config)
{
  try
  {
    observationMisc = ObservationMiscSstVariationalPtr(new ObservationMiscSstVariational(config));
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

void ObservationSstVariational::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    ObservationMiscSstVariational::Arc arc = observationMisc->computeArc(arcNo, covSst, covPod1, covPod2);
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
