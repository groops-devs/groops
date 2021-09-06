/***********************************************/
/**
* @file gnssParametrizationConstraints.h
*
* @brief Parameter constraints.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTS__
#define __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationConstraints = R"(
\subsection{Constraints}\label{gnssParametrizationType:constraints}
Add a pseudo observation equation (constraint)
for each selected \configClass{parameters}{parameterSelectorType}
\begin{equation}
  b-x_0 = 1 \cdot dx + \epsilon,
\end{equation}
where $b$ is the \config{bias} and $x_0$ is the a priori value of the parameter
if \config{relativeToApriori} is not set.
The standard deviation \config{sigma} is used to weight the observation equations.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Parameter constraints.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationConstraints : public GnssParametrizationBase
{
  std::string          name;
  ParameterSelectorPtr parameterSelector;
  Double               sigma;
  Double               bias;
  Bool                 relativeToApriori;
  Gnss                *gnss;

public:
  GnssParametrizationConstraints(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/) {this->gnss = gnss;}
  void constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
};

/***********************************************/

inline GnssParametrizationConstraints::GnssParametrizationConstraints(Config &config)
{
  try
  {
    readConfig(config, "name",              name,              Config::OPTIONAL, "constraint.name",  "");
    readConfig(config, "parameters",        parameterSelector, Config::MUSTSET,  "",  "parameter to constrain");
    readConfig(config, "sigma",             sigma,             Config::MUSTSET,  "",  "sigma of the constraint (same unit as parameter)");
    readConfig(config, "bias",              bias,              Config::DEFAULT,  "0", "constrain all selected parameters towards this value");
    readConfig(config, "relativeToApriori", relativeToApriori, Config::DEFAULT,  "0", "constrain only dx and not full x=dx+x0");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssParametrizationConstraints::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    Vector x0 = Vector(normalEquationInfo.parameterCount());
    if(!relativeToApriori)
      x0 = gnss->aprioriParameter(normalEquationInfo);

    const Double weight = 1./std::pow(sigma, 2);
    const std::vector<UInt> indices = parameterSelector->indexVector(normalEquationInfo.parameterNames());
    UInt count = 0;
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
        count++;
      }

    if(count)
      logStatus<<"constrain "<<name<<" ("<<count<<" parameters)"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
