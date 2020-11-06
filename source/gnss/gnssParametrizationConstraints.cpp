/***********************************************/
/**
* @file gnssParametrizationConstraints.cpp
*
* @brief GNSS parameter constraints.
*
* @author Torsten Mayer-Guerr
* @date 2019-05-28
*
*/
/***********************************************/

#define DOCSTRING_GnssParametrizationConstraints

#include "base/import.h"
#include "config/configRegister.h"
#include "parallel/matrixDistributed.h"
#include "gnss/gnss.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssParametrizationConstraints.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(GnssParametrizationConstraints, "gnssParametrizationConstraintsType")
GROOPS_READCONFIG_CLASS(GnssParametrizationConstraints, "gnssParametrizationConstraintsType")

/***********************************************/

GnssParametrizationConstraints::GnssParametrizationConstraints(Config &config, const std::string &name)
{
  try
  {
    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "comment",           comment,           Config::OPTIONAL, "",  "");
    readConfig(config, "parameters",        parameterSelector, Config::MUSTSET,  "",  "parameter to constrain");
    readConfig(config, "sigma",             sigma,             Config::MUSTSET,  "",  "sigma of the constraint (same unit as parameter)");
    readConfig(config, "bias",              bias,              Config::DEFAULT,  "0", "constrain all selected parameters towards this value");
    readConfig(config, "relativeToApriori", relativeToApriori, Config::DEFAULT,  "0", "constrain only dx and not full x=dx+x0");
    endSequence(config);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationConstraints::observationEquation(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_OTHER))
      return;

    Vector x0 = Vector(normalEquationInfo.parameterCount());
    if(!relativeToApriori)
      x0 = gnss().aprioriParameter(normalEquationInfo);

    const Double weight = 1./std::pow(sigma, 2);
    const std::vector<UInt> indices = parameterSelector->indexVector(normalEquationInfo.parameterNames());
    for(UInt index : indices)
      if(index != NULLINDEX)
      {
        const UInt idBlock    = normals.index2block(index);
        const UInt blockIndex = normals.blockIndex(idBlock);
        normals.setBlock(idBlock, idBlock);
        if(normals.isMyRank(idBlock, idBlock))
          normals.N(idBlock, idBlock)(index-blockIndex, index-blockIndex) += weight;
        if(Parallel::isMaster(normalEquationInfo.comm))
        {
          n.at(idBlock)(index-blockIndex, 0) += weight * (bias - x0.at(index));
          lPl += weight * std::pow(bias - x0.at(index), 2);
          obsCount++;
        }
      }

    if(indices.size())
      logStatus<<"constrain "<<comment<<" ("<<indices.size()<<" parameters)"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
