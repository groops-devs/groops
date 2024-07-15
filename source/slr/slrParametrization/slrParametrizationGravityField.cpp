/***********************************************/
/**
* @file slrParametrizationGravityField.cpp
*
* @brief GravityField.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "base/planets.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "slr/slrParametrization/slrParametrizationGravityField.h"

/***** CLASS ***********************************/

/** @brief GravityField.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */

/***********************************************/

SlrParametrizationGravityField::SlrParametrizationGravityField(Config &config)
{
  try
  {
    readConfig(config, "name",                name,            Config::OPTIONAL, "parameter.gravityField", "used for parameter selection");
//     readConfig(config, "outputfileParameter", fileName,        Config::OPTIONAL, "", "");
    readConfig(config, "parametrization",     parametrization, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationGravityField::getParametrizationGravity(std::vector<const SlrParametrizationGravityField*> &paramGravityField) const
{
  paramGravityField.push_back(this);
}

/***********************************************/

void SlrParametrizationGravityField::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    this->slr = slr;
    this->x   = Vector(parametrization->parameterCount());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationGravityField::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    indexParameter = SlrParameterIndex();
    if(!isEnabled(normalEquationInfo, name))
      return;

    std::vector<ParameterName> parameterNames;
    parametrization->parameterName(parameterNames);
    for(auto &parameterName : parameterNames)
      parameterName.object = "gravityfield";
    indexParameter = normalEquationInfo.parameterNamesOther(parameterNames);

    logInfo<<parametrization->parameterCount()%"%9i gravity field parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationGravityField::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(indexParameter)
      copy(x, x0.row(normalEquationInfo.index(indexParameter), x.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SlrParametrizationGravityField::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    if(indexParameter)
      this->x += x.row(normalEquationInfo.index(indexParameter), parametrization->parameterCount());
    return 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationGravityField::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &/*suffix*/) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

//     if(!fileName.empty())
//     {
//       logStatus<<"write sphercial haramonics to files <"<<fileName.appendBaseName(suffix)<<">"<<Log::endl;
//       writeFileSphericalHarmonics(fileName.appendBaseName(suffix), parametrization->sphericalHarmonics(Time(), x, sigma2x));
//     }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
