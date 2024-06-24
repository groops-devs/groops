/***********************************************/
/**
* @file slrParametrizationConstraints.h
*
* @brief Parameter constraints.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONCONSTRAINTS__
#define __GROOPS_SLRPARAMETRIZATIONCONSTRAINTS__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationConstraints = R"(
\subsection{Constraints}\label{slrParametrizationType:constraints}
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
#include "slr/slr.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Parameter constraints.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationConstraints : public SlrParametrizationBase
{
  std::string          name;
  ParameterSelectorPtr parameterSelector;
  Double               sigma;
  Double               bias;
  Bool                 relativeToApriori;
  Slr                 *slr;

public:
  SlrParametrizationConstraints(Config &config);

  void init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/) {this->slr = slr;}
  void constraints(const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
};

/***********************************************/

inline SlrParametrizationConstraints::SlrParametrizationConstraints(Config &config)
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

inline void SlrParametrizationConstraints::constraints(const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    Vector x0 = Vector(normalEquationInfo.parameterCount());
    if(!relativeToApriori)
      x0 = slr->aprioriParameter(normalEquationInfo);

    const Double weight = 1./std::pow(sigma, 2);
    const std::vector<UInt> indices = parameterSelector->indexVector(normalEquationInfo.parameterNames());
    UInt count = 0;
    for(UInt index : indices)
      if(index != NULLINDEX)
      {
        const UInt idBlock    = normals.index2block(index);
        const UInt blockIndex = normals.blockIndex(idBlock);
        normals.setBlock(idBlock, idBlock);
        normals.N(idBlock, idBlock)(index-blockIndex, index-blockIndex) += weight;
        n.at(idBlock)(index-blockIndex, 0) += weight * (bias - x0.at(index));
        lPl += weight * std::pow(bias - x0.at(index), 2);
        obsCount++;
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
