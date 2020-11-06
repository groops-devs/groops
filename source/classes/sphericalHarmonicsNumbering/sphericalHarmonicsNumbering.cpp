/***********************************************/
/**
* @file sphericalHarmonicsNumbering.cpp
*
* @brief Numbering schema of spherical harmonics coefficients.
*
* @author Torsten Mayer-Guerr
* @date 2009-01-23
*
*/
/***********************************************/

#define DOCSTRING_SphericalHarmonicsNumbering

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumberingDegree.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumberingOrder.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumberingOrderNonAlternating.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumberingFile.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***********************************************/

GROOPS_REGISTER_CLASS(SphericalHarmonicsNumbering, "sphericalHarmonicsNumberingType",
                      SphericalHarmonicsNumberingDegree,
                      SphericalHarmonicsNumberingOrder,
                      SphericalHarmonicsNumberingOrderNonAlternating,
                      SphericalHarmonicsNumberingFile)

GROOPS_READCONFIG_CLASS(SphericalHarmonicsNumbering, "sphericalHarmonicsNumberingType")

/***********************************************/

SphericalHarmonicsNumberingPtr SphericalHarmonicsNumbering::create(Config &config, const std::string &name)
{
  try
  {
    SphericalHarmonicsNumberingPtr numbering;
    std::string type;
    readConfigChoice(config, name, type, Config::MUSTSET, "", "Numbering schema of spherical harmonics coefficients");

    if(readConfigChoiceElement(config, "degreewise",  type, "sort degree by degree"))
      numbering = SphericalHarmonicsNumberingPtr(new SphericalHarmonicsNumberingDegree(config));
    if(readConfigChoiceElement(config, "orderwise",   type, "sort order by order"))
      numbering = SphericalHarmonicsNumberingPtr(new SphericalHarmonicsNumberingOrder(config));
    if(readConfigChoiceElement(config, "orderwiseNonAlternating",  type, "sort order by order with cnm, snm non-alternating"))
      numbering = SphericalHarmonicsNumberingPtr(new SphericalHarmonicsNumberingOrderNonAlternating(config));
    if(readConfigChoiceElement(config, "file",  type, "sort as specified in the chosen file"))
      numbering = SphericalHarmonicsNumberingPtr(new SphericalHarmonicsNumberingFile(config));

    endChoice(config);

    return numbering;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SphericalHarmonicsNumbering::numbering(UInt maxDegree, UInt minDegree, std::vector<UInt> &degree, std::vector<UInt> &order, std::vector<UInt> &cs) const
{
  try
  {
    const UInt paraCount = parameterCount(maxDegree, minDegree);
    degree.resize(paraCount);
    order.resize(paraCount);
    cs.resize(paraCount);

    std::vector<std::vector<UInt>> idxC, idxS;
    numbering(maxDegree, minDegree, idxC, idxS);

    for(UInt n=minDegree; n<=maxDegree; n++)
      for(UInt m=0; m<=n; m++)
      {
        if(idxC[n][m] != NULLINDEX)
        {
          degree[idxC[n][m]] = n;
          order[idxC[n][m]]  = m;
          cs[idxC[n][m]]     = 0;
        }

        if(idxS[n][m] != NULLINDEX)
        {
          degree[idxS[n][m]] = n;
          order[idxS[n][m]]  = m;
          cs[idxS[n][m]]     = 1;
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
