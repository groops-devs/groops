/***********************************************/
/**
* @file sphericalHarmonicsFilter.cpp
*
* @brief Filtering spherical harmonics.
*
* @author Torsten Mayer-Guerr
* @date 2008-06-08
*
*/
/***********************************************/

#define DOCSTRING_SphericalHarmonicsFilter

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilterDdk.h"
#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilterGauss.h"
#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilterMatrix.h"
#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilter.h"

/***********************************************/

GROOPS_REGISTER_CLASS(SphericalHarmonicsFilter, "sphericalHarmonicsFilterType",
                      SphericalHarmonicsFilterDdk,
                      SphericalHarmonicsFilterGauss,
                      SphericalHarmonicsFilterMatrix)

GROOPS_READCONFIG_UNBOUNDED_CLASS(SphericalHarmonicsFilter, "sphericalHarmonicsFilterType")

/***********************************************/

SphericalHarmonicsFilter::SphericalHarmonicsFilter(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "filtering of spherical harmonics"))
    {
      if(readConfigChoiceElement(config, "filterDdk",       type, "smoothing by a DDK filter"))
        filters.push_back(new SphericalHarmonicsFilterDdk(config));
      if(readConfigChoiceElement(config, "filterGauss",       type, "smoothing by a isotropic gaussian filter"))
        filters.push_back(new SphericalHarmonicsFilterGauss(config));
      if(readConfigChoiceElement(config, "filterMatrix",      type, "smoothing by a filter matrix"))
        filters.push_back(new SphericalHarmonicsFilterMatrix(config));
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

SphericalHarmonicsFilter::~SphericalHarmonicsFilter()
{
  for(UInt i=0; i<filters.size(); i++)
    delete filters.at(i);
}

/***********************************************/

SphericalHarmonics SphericalHarmonicsFilter::filter(const SphericalHarmonics &harm) const
{
  try
  {
    SphericalHarmonics harm2 = harm;
    for(UInt i=0; i<filters.size(); i++)
      harm2 = filters.at(i)->filter(harm2);
    return harm2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
