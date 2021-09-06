/***********************************************/
/**
* @file gnssProcessingStepEstimate.h
*
* @brief GNSS processing step: Estimate.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPESTIMATE__
#define __GROOPS_GNSSPROCESSINGSTEPESTIMATE__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepEstimate = R"(
\subsection{Estimate}\label{gnssProcessingStepType:estimate}
Iterative non-linear least squares adjustment.
In every iteration it accumulates the system of normal equations, solves the system and updates the estimated parameters.
The estimated parameters serve as a priori values in the next iteration and the following processing steps.
Iterates until either every single parameter update (converted to an influence in meter)
is below a \config{convergenceThreshold} or \config{maxIterationCount} is reached.

With \config{computeResiduals} the observation equations are computed
again after each update to compute the observation residuals.

The overall standard deviation of a single observation used for the weighting
is composed of several factors
\begin{equation}
  \hat{\sigma}_i = \hat{\sigma}_i^{huber} \hat{\sigma}_{[\tau\nu a]}^{recv} \sigma_{[\tau\nu a]}^{recv}(E,A),
\end{equation}
where $[\tau\nu a]$ is the signal type, the azmiuth and elevation dependent $\sigma_{[\tau\nu a]}^{recv}(E,A)$ is given by
\configFile{receiver:inputfileAccuracyDefinition}{gnssAntennaDefinition} and the other factors are
estimated iteratively from the residuals.

With \config{computeWeights} a standardized variance $\hat{s}_i^2$
for each residual $\hat{\epsilon}_i$ is computed
\begin{equation}
  \hat{s}_i^2 = \frac{1}{\hat{\sigma}_{[\tau\nu a]}^{recv} \sigma_{[\tau\nu a]}^{recv}(E,A)}\frac{\hat{\epsilon}_i^2}{r_i}
  \qquad\text{with}\qquad
  r_i = \left(\M A\left(\M A^T\M A\right)^{-1}\M A^T\right)_{ii}
\end{equation}
taking the redundancy $r_i$ into account. If $\hat{s}_i$ is above a threshold \config{huber}
the observation gets a higher standard deviation used for weighting according to
\begin{equation}
  \hat{\sigma}_i^{huber} =
  \left\{ \begin{array}{ll}
    1                              & s < huber,\\
    (\hat{s}_i/huber)^{huberPower} & s \ge huber
  \end{array} \right.,
\end{equation}
similar to \reference{robust least squares adjustment}{fundamentals.robustLeastSquares}.

With \config{adjustSigma0} individual variance factors can be computed
for each station and all phases of a system and each code observation \reference{type}{gnssType}
(e.g. for each \verb|L**G|, \verb|L**E|, \verb|C1CG|, \verb|C2WG|, \verb|C1CE|, \ldots)
separately
\begin{equation}
  \hat{\sigma}_{[\tau\nu a]}^{recv} = \sqrt{\frac{\hat{\M\epsilon}^T\M P\hat{\M\epsilon}}{r}}.
\end{equation}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: Estimate.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepEstimate : public GnssProcessingStepBase
{
  Bool   computeResiduals, adjustSigma0, computeWeights;
  Double huber, huberPower;
  Double convergenceThreshold;
  UInt   iterCount;

public:
  GnssProcessingStepEstimate(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepEstimate::GnssProcessingStepEstimate(Config &config)
{
  try
  {
    readConfig(config, "computeResiduals",     computeResiduals,     Config::DEFAULT, "1",    "");
    readConfig(config, "adjustSigma0",         adjustSigma0,         Config::DEFAULT, "1",    "adjust sigma0 by scale factor (per receiver and type)");
    readConfig(config, "computeWeights",       computeWeights,       Config::DEFAULT, "1",    "downweight outliers");
    readConfig(config, "huber",                huber,                Config::DEFAULT, "2.5",  "residuals > huber*sigma0 are downweighted");
    readConfig(config, "huberPower",           huberPower,           Config::DEFAULT, "1.5",  "residuals > huber: sigma=(e/huber)^huberPower*sigma0");
    readConfig(config, "convergenceThreshold", convergenceThreshold, Config::DEFAULT, "0.01", "[m] stop iteration once full convergence is reached");
    readConfig(config, "maxIterationCount",    iterCount,            Config::DEFAULT, "3",    "maximum number of iterations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessingStepEstimate::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== estimate ================================================"<<Log::endl;
    if(state.changedNormalEquationInfo)
      state.gnss->initParameter(state.normalEquationInfo);
    state.changedNormalEquationInfo = FALSE;
    for(UInt iter=0; iter<iterCount; iter++)
    {
      logStatus<<iter+1<<". iteration  --------------------------"<<Log::endl;
      if(convergenceThreshold > state.estimateSolution(nullptr/*resolveAmbiguities*/, computeResiduals, computeWeights, adjustSigma0, huber, huberPower))
        break;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
