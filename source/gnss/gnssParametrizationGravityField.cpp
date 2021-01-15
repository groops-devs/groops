/***********************************************/
/**
* @file gnssParametrizationGravityField.cpp
*
* @brief GNSS gravity field.
*
* @author Torsten Mayer-Guerr
* @date 2014-05-26
*
*/
/***********************************************/

#define DOCSTRING_GnssParametrizationGravityField

#include "base/import.h"
#include "config/config.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "gnss/gnss.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssParametrizationGravityField.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(GnssParametrizationGravityField, "gnssParametrizationGravityFieldType")
GROOPS_READCONFIG_CLASS(GnssParametrizationGravityField, "gnssParametrizationGravityFieldType")

/***********************************************/

GnssParametrizationGravityField::GnssParametrizationGravityField(Config &config, const std::string &name)
{
  try
  {
    readConfigSequence(config, name, Config::MUSTSET, "", "");
    renameDeprecatedConfig(config, "representation", "parametrization", date2time(2020, 6, 3));
    readConfig(config, "parametrization", parametrizationPtr, Config::MUSTSET, "", "gravity field parametrization");
    endSequence(config);
    if(isCreateSchema(config)) return;

    x = Vector(parametrizationPtr->parameterCount());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationGravityField::initParameter(Gnss::NormalEquationInfo &normalEquationInfo)
{
  try
  {
    index = Gnss::ParameterIndex();

    if(!(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_GRAVITY))
      return;

    std::vector<ParameterName> parameterNames;
    parametrizationPtr->parameterName(parameterNames);
    index = normalEquationInfo.parameterNamesOther(parameterNames);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationGravityField::aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(index)
      copy(x, x0.row(normalEquationInfo.index(index), parametrizationPtr->parameterCount()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationGravityField::updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/, Bool /*printStatistics*/)
{
  try
  {
    if(index)
      this->x += x.row(normalEquationInfo.index(index), parametrizationPtr->parameterCount());
    return 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationGravityField::writeResults(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, const std::string &/*suffix*/)
{
  try
  {
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
