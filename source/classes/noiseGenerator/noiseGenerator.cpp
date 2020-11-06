/***********************************************/
/**
* @file noiseGenerator.cpp
*
* @brief Generate different types of noise.
*
* @author Matthias Ellmer
* @date 2013-09-18
*
*/
/***********************************************/

#define DOCSTRING_NoiseGenerator

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "classes/noiseGenerator/noiseGeneratorWhite.h"
#include "classes/noiseGenerator/noiseGeneratorExpressionPSD.h"
#include "classes/noiseGenerator/noiseGeneratorDigitalFilter.h"
#include "classes/noiseGenerator/noiseGeneratorPowerLaw.h"

/***********************************************/

GROOPS_REGISTER_CLASS(NoiseGenerator, "noiseGeneratorType",
                      NoiseGeneratorWhite,
                      NoiseGeneratorExpressionPSD,
                      NoiseGeneratorDigitalFilter,
                      NoiseGeneratorPowerLaw)

GROOPS_READCONFIG_UNBOUNDED_CLASS(NoiseGenerator, "noiseGeneratorType")

/***********************************************/

NoiseGenerator::NoiseGenerator(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "generate different types of noise"))
    {
      if(readConfigChoiceElement(config, "white",         type, "white noise"))
        noiseGenerator.push_back(new NoiseGeneratorWhite(config));
      if(readConfigChoiceElement(config, "expressionPSD", type, "defined by one sided PSD"))
        noiseGenerator.push_back(new NoiseGeneratorExpressionPSD(config));
      if(readConfigChoiceElement(config, "filter",        type, "noise conforms to psd of filter"))
        noiseGenerator.push_back(new NoiseGeneratorDigitalFilter(config));
      if(readConfigChoiceElement(config, "powerLaw",      type, "noise follows power law relationship (f^alpha)"))
        noiseGenerator.push_back(new NoiseGeneratorPowerLaw(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

NoiseGenerator::~NoiseGenerator()
{
  for(UInt i=0; i<noiseGenerator.size(); i++)
    delete noiseGenerator.at(i);
}

/***********************************************/

Matrix NoiseGenerator::noise(UInt samples, UInt series) const
{
  try
  {
    Matrix sum(samples, series);
    for(UInt i=0; i<noiseGenerator.size(); i++)
      sum += noiseGenerator.at(i)->noise(samples, series);
    return sum;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix NoiseGenerator::covarianceFunction(UInt length, Double sampling) const
{
  try
  {
    Vector sum(length);

    for(UInt i=0; i<noiseGenerator.size(); i++)
      sum += noiseGenerator.at(i)->covarianceFunction(length, sampling);

    Matrix covariance(length,2);
    for(UInt i=0; i<length; i++)
      covariance(i,0) = i*sampling;

    copy(sum, covariance.column(1));
    return covariance;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
