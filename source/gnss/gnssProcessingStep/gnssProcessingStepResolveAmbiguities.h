/***********************************************/
/**
* @file gnssProcessingStepResolveAmbiguities.h
*
* @brief GNSS processing step: ResolveAmbiguities.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPRESOLVEAMBIGUITIES__
#define __GROOPS_GNSSPROCESSINGSTEPRESOLVEAMBIGUITIES__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepResolveAmbiguities = R"(
\subsection{ResolveAmbiguities}\label{gnssProcessingStepType:resolveAmbiguities}
Performs a least squares adjustment like \configClass{processingStep:estimate}{gnssProcessingStepType:estimate}
but with additional integer phase ambiguity resolution.
After this step all resolved ambiguities are removed from the normal equation system.

Integer ambiguity resolution is performed based on the least squares ambiguity decorrelation adjustment
(LAMBDA) method (Teunissen 1995, DOI \href{https://doi.org/10.1007/BF00863419}{10.1007/BF00863419}), specifically
the modified algorithm (MLAMBDA) by Chang et al. (2005, DOI \href{https://doi.org/10.1007/s00190-005-0004-x}{10.1007/s00190-005-0004-x}).
First the covariance matrix of the integer ambiguity parameters is computed by eliminating all but those parameters
from the full normal equation matrix and inverting it. Then, a Z-transformation is performed as described by
Chang et al. (2005) to decorrelate the ambiguity parameters without losing their integer nature.

The search process follows MLAMBDA and uses integer minimization of the weighted sum of squared residuals.
It is computationally infeasible to search a hyper-ellipsoid with a dimension of ten thousand or more.
Instead, a blocked search algorithm is performed by moving a window with a length of, for example,
\config{searchBlockSize}=\verb|200| parameters over the decorrelated ambiguities, starting from the most accurate.
In each step, the window is moved by half of its length and the overlapping parts are compared to each other.
If all fixed ambiguities in the overlap agree, the algorithm continues.
Otherwise, both windows are combined and the search is repeated using the combined window, again comparing with the overlapping
part of the preceding window. If not all solutions could be checked for a block after \config{maxSearchSteps},
the selected \config{incompleteAction} is performed.
If the algorithm reaches ambiguities with a standard deviation higher than \config{sigmaMaxResolve},
ambiguity resolution stops and the remaining ambiguities are left as float values.
Otherwise, all ambiguity parameters are fixed to integer values.

In contrast to an integer least squares solution over the full ambiguity vector, it is not guaranteed that the resulting solution
is optimal in the sense of minimal variance with given covariance.
This trade-off is necessary to cope with large numbers of ambiguities.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileMatrix.h"
#include "gnss/gnssLambda.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: ResolveAmbiguities.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepResolveAmbiguities : public GnssProcessingStepBase
{
  FileName                      fileNameAmbiguities;
  Double                        sigmaMaxResolve;
  UInt                          searchBlockSize, maxSearchSteps;
  GnssLambda::IncompleteAction  incompleteAction;
  Bool                          computeResiduals, adjustSigma0, computeWeights;
  Double                        huber, huberPower;

public:
  GnssProcessingStepResolveAmbiguities(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepResolveAmbiguities::GnssProcessingStepResolveAmbiguities(Config &config)
{
  try
  {
    std::string choice;

    readConfig(config, "outputfileAmbiguities", fileNameAmbiguities, Config::OPTIONAL, "",    "resolved ambiguities");
    readConfig(config, "sigmaMaxResolve",       sigmaMaxResolve,     Config::OPTIONAL, "0.2", "max. allowed std. dev. of ambiguity to resolve [cycles]");
    readConfig(config, "searchBlockSize",       searchBlockSize,     Config::DEFAULT,  "200", "block size for blocked integer search");
    readConfig(config, "maxSearchSteps",        maxSearchSteps,      Config::OPTIONAL, "200000000", "max. steps of integer search for each block");
    if(readConfigChoice(config, "incompleteAction", choice, Config::MUSTSET, "shrinkBlockSize", "if not all solutions tested after maxSearchSteps"))
    {
      if(readConfigChoiceElement(config, "stop",            choice, "stop searching, ambiguities remain float in this block")) incompleteAction = GnssLambda::IncompleteAction::STOP;
      if(readConfigChoiceElement(config, "resolve",         choice, "use best integer solution found so far"))                 incompleteAction = GnssLambda::IncompleteAction::IGNORE;
      if(readConfigChoiceElement(config, "shrinkBlockSize", choice, "try again with half block size"))                         incompleteAction = GnssLambda::IncompleteAction::SHRINKBLOCKSIZE;
      if(readConfigChoiceElement(config, "throwException",  choice, "stop and throw an exception"))                            incompleteAction = GnssLambda::IncompleteAction::EXCEPTION;
      endChoice(config);
    }
    readConfig(config, "computeResiduals",      computeResiduals,    Config::DEFAULT,  "1",   "");
    readConfig(config, "adjustSigma0",          adjustSigma0,        Config::DEFAULT,  "1",   "adjust sigma0 by scale factor (per receiver and type)");
    readConfig(config, "computeWeights",        computeWeights,      Config::DEFAULT,  "1",   "downweight outliers");
    readConfig(config, "huber",                 huber,               Config::DEFAULT,  "2.5", "residuals > huber*sigma0 are downweighted");
    readConfig(config, "huberPower",            huberPower,          Config::DEFAULT,  "1.5", "residuals > huber: sigma=(e/huber)^huberPower*sigma0");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessingStepResolveAmbiguities::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== resolve ambiguities  ===================================="<<Log::endl;
    if(state.changedNormalEquationInfo)
      state.gnss->initParameter(state.normalEquationInfo);
    Matrix solutionSteps;
    state.estimateSolution(std::bind(&GnssLambda::searchIntegerBlocked, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                                     sigmaMaxResolve, searchBlockSize, maxSearchSteps, incompleteAction, TRUE/*timing*/,
                                     std::placeholders::_4, std::placeholders::_5, solutionSteps),
                           computeResiduals, computeWeights, adjustSigma0, huber, huberPower);
    state.changedNormalEquationInfo = TRUE;

    if(!fileNameAmbiguities.empty())
    {
      logStatus<<"write ambiguities to file(s) <"<<fileNameAmbiguities<<">"<<Log::endl;
      writeFileMatrix(fileNameAmbiguities, solutionSteps);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
