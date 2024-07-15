/***********************************************/
/**
* @file slrProcessingStep.cpp
*
* @brief Processing steps for SLR normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#define DOCSTRING_SlrProcessingStep

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "config/configRegister.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "misc/varianceComponentEstimation.h"
#include "slr/slr.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"
#include "slr/slrProcessingStep/slrProcessingStepEstimate.h"
#include "slr/slrProcessingStep/slrProcessingStepWriteResults.h"
#include "slr/slrProcessingStep/slrProcessingStepWriteNormalEquations.h"
#include "slr/slrProcessingStep/slrProcessingStepWriteAprioriSolution.h"
#include "slr/slrProcessingStep/slrProcessingStepWriteResiduals.h"
#include "slr/slrProcessingStep/slrProcessingStepWriteUsedStationList.h"
#include "slr/slrProcessingStep/slrProcessingStepWriteUsedSatelliteList.h"
#include "slr/slrProcessingStep/slrProcessingStepPrintResidualStatistics.h"
#include "slr/slrProcessingStep/slrProcessingStepSelectParametrizations.h"
#include "slr/slrProcessingStep/slrProcessingStepSelectSatellites.h"
#include "slr/slrProcessingStep/slrProcessingStepSelectStations.h"
#include "slr/slrProcessingStep/slrProcessingStepGroup.h"

/***********************************************/

GROOPS_REGISTER_CLASS(SlrProcessingStep, "slrProcessingStepType",
                      SlrProcessingStepEstimate,
                      SlrProcessingStepWriteResults,
                      SlrProcessingStepWriteNormalEquations,
                      SlrProcessingStepWriteAprioriSolution,
                      SlrProcessingStepWriteResiduals,
                      SlrProcessingStepWriteUsedStationList,
                      SlrProcessingStepWriteUsedSatelliteList,
                      SlrProcessingStepPrintResidualStatistics,
                      SlrProcessingStepSelectParametrizations,
                      SlrProcessingStepSelectSatellites,
                      SlrProcessingStepSelectStations,
                      SlrProcessingStepGroup)

GROOPS_READCONFIG_UNBOUNDED_CLASS(SlrProcessingStep, "slrProcessingStepType")

/***********************************************/

SlrProcessingStep::SlrProcessingStep(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", ""))
    {
      if(readConfigChoiceElement(config, "estimate",                       type, "least squares adjustment"))
        bases.push_back(new SlrProcessingStepEstimate(config));
      if(readConfigChoiceElement(config, "writeResults",                   type, "write all estimated parameters"))
        bases.push_back(new SlrProcessingStepWriteResults(config));
      if(readConfigChoiceElement(config, "writeNormalEquations",           type, "write unconstrained and constraint normal equations"))
        bases.push_back(new SlrProcessingStepWriteNormalEquations(config));
      if(readConfigChoiceElement(config, "writeAprioriSolution",           type, "write apriori solution vector"))
        bases.push_back(new SlrProcessingStepWriteAprioriSolution(config));
      if(readConfigChoiceElement(config, "writeResiduals",                 type, "write observation residuals"))
        bases.push_back(new SlrProcessingStepWriteResiduals(config));
      if(readConfigChoiceElement(config, "writeUsedStationList",           type, "write used stations"))
        bases.push_back(new SlrProcessingStepWriteUsedStationList(config));
      if(readConfigChoiceElement(config, "writeUsedSatelliteList",         type, "write used satellites"))
        bases.push_back(new SlrProcessingStepWriteUsedSatelliteList(config));
      if(readConfigChoiceElement(config, "printResidualStatistics",        type, "print residual statistics"))
        bases.push_back(new SlrProcessingStepPrintResidualStatistics(config));
      if(readConfigChoiceElement(config, "selectParametrizations",         type, "select parametrizations for all subsequent processing steps"))
        bases.push_back(new SlrProcessingStepSelectParametrizations(config));
      if(readConfigChoiceElement(config, "selectSatellites",               type, "use this subset of satellites in all subsequent processing steps"))
        bases.push_back(new SlrProcessingStepSelectSatellites(config));
      if(readConfigChoiceElement(config, "selectStations",                 type, "use this subset of stations in all subsequent processing steps"))
        bases.push_back(new SlrProcessingStepSelectStations(config));
      if(readConfigChoiceElement(config, "group",                          type, "group processing steps"))
        bases.push_back(new SlrProcessingStepGroup(config));
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

SlrProcessingStep::~SlrProcessingStep()
{
  for(auto base : bases)
    delete base;
}

/***********************************************/

void SlrProcessingStep::process(SlrProcessingStep::State &state)
{
  try
  {
    for(auto base : bases)
    {
      if(base->expectInitializedParameters() && state.changedNormalEquationInfo)
      {
        state.slr->initParameter(state.normalEquationInfo);
        state.changedNormalEquationInfo = FALSE;
      }
      base->process(state);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

SlrProcessingStep::State::State(SlrPtr slr) :
  slr(slr), normalEquationInfo(slr->stations.size(), slr->satellites.size()),
  changedNormalEquationInfo(TRUE), sigmaFactor(slr->stations.size(), 1.0)
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

void SlrProcessingStep::State::regularizeNotUsedParameters(UInt blockStart, UInt blockCount)
{
  try
  {
    UInt countZeros = 0;
    for(UInt i=blockStart; i<blockStart+blockCount; i++)
      if(normals.isMyRank(i,i))
      {
        Matrix &N = normals.N(i,i);
        for(UInt k=0; k<N.rows(); k++)
          if(N(k,k)==0.)
          {
            N(k,k) += 1.0;
            countZeros++;
            logWarning<<"    "<<normalEquationInfo.parameterNames().at(normals.blockIndex(i)+k).str()<<" has zero diagonal element"<<Log::endl;
          }
      }
    if(countZeros)
    {
      logWarning<<"  "<<countZeros<<" parameters have zero diagonal elements -> set to one"<<Log::endl;
      obsCount -= countZeros;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrProcessingStep::State::buildNormals(Bool constraintsOnly)
{
  try
  {
    normals.initEmpty(normalEquationInfo.blockIndices(), Parallel::selfCommunicator());
    n.resize(normals.blockCount());
    for(UInt i=0; i<normals.blockCount(); i++)
      n.at(i) = Vector(normals.blockSize(i));
    lPl = Vector(1);
    obsCount = 0;

    if(!normalEquationInfo.parameterCount())
      return;

    // Loop over all stations
    // ----------------------
    if(!constraintsOnly)
    {
      logStatus<<"accumulate normals"<<Log::endl;
      SlrDesignMatrix A(normalEquationInfo);
      SlrObservationEquation eqn;
      UInt idLoop = 0;
      Log::Timer timer(slr->stations.size());
      for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      {
        timer.loopStep(idLoop++);
        auto station = slr->stations.at(idStat);
        if(normalEquationInfo.estimateStation.at(idStat))
          for(UInt idSat=0; idSat<station->observations.size(); idSat++)
            if(normalEquationInfo.estimateSatellite.at(idSat))
              for(UInt idPass=0; idPass<station->observations.at(idSat).size(); idPass++)
                if(slr->basicObservationEquations(normalEquationInfo, idStat, idSat, idPass, eqn))
                {
                  A.init(eqn.l);
                  slr->designMatrix(normalEquationInfo, eqn, A);
                  A.accumulateNormals(normals, n, lPl(0), obsCount);
                }
      } // for(idStat)
      timer.loopEnd();
    }

    // other observations and constraints
    // ----------------------------------
    slr->constraints(normalEquationInfo, normals, n, lPl(0), obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SlrProcessingStep::State::estimateSolution(Bool computeResiduals, Bool computeWeights, Bool adjustSigma0, Double huber, Double huberPower)
{
  try
  {
    // setup normal equations
    buildNormals(FALSE/*constraintsOnly*/);

    if(normalEquationInfo.parameterCount())
    {
      logStatus<<"solve"<<Log::endl;
      regularizeNotUsedParameters(0, normals.blockCount());
      normals.cholesky(TRUE/*timing*/);
      normals.triangularTransSolve(n); // forward step
      for(UInt z=0; z<normals.blockCount(); z++)
      {
        lPl(0)   -= quadsum(n.at(z)); // lPl = lPl - n1' N1^(-1) n1
        obsCount -= normals.blockSize(z);
      }
      normals.triangularSolve(n); // Backward step

      // Compute sigma
      // -------------
      const Double sigma = Vce::standardDeviation(lPl(0), obsCount, huber, huberPower);
      logInfo<<"  sigma = "<<sigma%"%.2f"s<<Log::endl;
      if((sigma!=sigma) || (sigma<=0))
        logWarning<<"  Cannot compute sigma = sqrt("<<lPl(0)<<"/"<<obsCount<<")"<<Log::endl;
    }

    // generate monte carlo vectors
    // ----------------------------
    Matrix Wz;
    if(computeResiduals)
    {
      constexpr UInt monteCarloColumns = 100;
      std::vector<Matrix> monteCarlo(normals.blockCount()); // Monte Carlo Vector for redundancy computation
      for(UInt i=0; i<normals.blockCount(); i++)
        monteCarlo.at(i) = Vce::monteCarlo(normals.blockSize(i), monteCarloColumns);

      normals.triangularSolve(monteCarlo);

      Wz = Matrix(normalEquationInfo.parameterCount(), monteCarloColumns);
      for(UInt i=0; i<normalEquationInfo.blockCount(); i++)
        if(normalEquationInfo.blockSize(i))
          copy(monteCarlo.at(i), Wz.row(normalEquationInfo.blockIndex(i), normalEquationInfo.blockSize(i)));
    }

    // free memory
    normals.initEmpty(normals.blockIndex(), normals.communicator());

    // copy solution to one vector
    // ---------------------------
    Vector x = Vector(normalEquationInfo.parameterCount());
    for(UInt i=0; i<normalEquationInfo.blockCount(); i++)
      if(normalEquationInfo.blockSize(i))
        copy(n.at(i), x.row(normalEquationInfo.blockIndex(i), normalEquationInfo.blockSize(i)));
    n.clear();
    n.shrink_to_fit();
    if(!computeResiduals)
    {
      if(normalEquationInfo.parameterCount())
        logInfo<<"Parameter changes"<<Log::endl;
      return slr->updateParameter(normalEquationInfo, x, Matrix());
    }

    // Residual tracking
    // -----------------
    logStatus<<"Compute residuals"<<Log::endl;
    Slr::InfoParameterChange infosResiduals("mm");
    SlrDesignMatrix A(normalEquationInfo);
    SlrObservationEquation eqn;
    UInt idLoop = 0;
    Log::Timer timer(slr->stations.size());
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
    {
      timer.loopStep(idLoop++);
      auto station = slr->stations.at(idStat);
      if(normalEquationInfo.estimateStation.at(idStat))
      {
        for(UInt idSat=0; idSat<station->observations.size(); idSat++)
          if(normalEquationInfo.estimateSatellite.at(idSat))
            for(UInt idPass=0; idPass<station->observations.at(idSat).size(); idPass++)
              if(slr->basicObservationEquations(normalEquationInfo, idStat, idSat, idPass, eqn))
              {
                A.init(eqn.l);
                slr->designMatrix(normalEquationInfo, eqn, A);
                Vector We  = eqn.l - A.mult(x); // decorrelated residuals
                Matrix AWz = A.mult(Wz);        // redundancies

                // redundancies
                // ------------
                Vector r(We.rows());
                for(UInt i=0; i<We.rows(); i++)
                  r(i) = 1. - quadsum(AWz.row(i));

                // find max. residual (for statistics)
                // -----------------------------------
                for(UInt k=0; k<We.rows(); k++)
                  if(std::fabs(eqn.sigmas(k)-eqn.sigmas0(k)) < 1e-8) // without outlier
                    if(infosResiduals.update(1e3*(We(k)*eqn.sigmas(k) - station->observations.at(idSat).at(idPass)->residuals(k))))
                      infosResiduals.info = eqn.station->name()+" -> "+eqn.satellite->name()+", "+eqn.timesTrans.at(k).dateTimeStr();

                station->observations.at(idSat).at(idPass)->setDecorrelatedResiduals(We, r);
              }
      }
    } // for(idStat)
    timer.loopEnd();

    // new weights
    // -----------
    if(computeWeights || adjustSigma0)
    {
      if(computeWeights) logStatus<<"Downweight outliers"<<Log::endl;
      if(adjustSigma0)   logStatus<<"Estimate variance factors"<<Log::endl;
      for(auto station : slr->stations)
        if(normalEquationInfo.estimateStation.at(station->idStat()))
          for(UInt iter=0; iter<10; iter++)
          {
            if(adjustSigma0)
            {
              Double ePe = 0, redundancy = 0;
              for(UInt idSat=0; idSat<station->observations.size(); idSat++)
                if(normalEquationInfo.estimateSatellite.at(idSat))
                  for(auto &obs : station->observations.at(idSat))
                    for(UInt i=0; i<obs->residuals.rows(); i++)
                      if((obs->sigmas0(i) > 0) && (obs->redundancies(i) > 0))
                      {
                        ePe        += std::pow(obs->residuals(i)/obs->sigmas(i), 2);
                        redundancy += obs->redundancies(i);
                      }
              const Double factor = (redundancy > 3) ? Vce::standardDeviation(ePe, redundancy, huber, huberPower) : 1.;
              sigmaFactor.at(station->idStat()) *= factor;

              // adjust sigma and sigma0
              for(UInt idSat=0; idSat<station->observations.size(); idSat++)
                if(normalEquationInfo.estimateSatellite.at(idSat))
                  for(auto &obs : station->observations.at(idSat))
                    for(UInt i=0; i<obs->residuals.rows(); i++)
                    {
                      obs->sigmas0(i) *= factor;
                      obs->sigmas(i)  *= factor;
                    }
            } // if(adjustSigma0)

            if(computeWeights)
            {
              for(UInt idSat=0; idSat<station->observations.size(); idSat++)
                if(normalEquationInfo.estimateSatellite.at(idSat))
                  for(auto &obs : station->observations.at(idSat))
                    for(UInt i=0; i<obs->residuals.rows(); i++)
                    {
                      obs->sigmas(i) = obs->sigmas0(i);
                      if(obs->redundancies(i) > 0.1)
                      {
                        const Double s = std::fabs(obs->residuals(i))/std::sqrt(obs->redundancies(i))/obs->sigmas0(i);
                        obs->sigmas(i) *= (s > huber) ? std::pow(s/huber, huberPower) : 1.;
                      }
                    }
            } // if(computeWeights)

            if(!(computeWeights && adjustSigma0))
              break; // iteration not needed
          } // for(station, iter)
    } // if(computeWeights || adjustSigma0)

    // residual analysis
    // -----------------
    logInfo<<"Residuals changes (without outliers)"<<Log::endl;
    Double maxChange = 0;
    infosResiduals.print(1e-3, maxChange);

    if(normalEquationInfo.parameterCount())
      logInfo<<"Parameter changes"<<Log::endl;
    return slr->updateParameter(normalEquationInfo, x, Wz);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrProcessingStep::State::residualsStatistics(UInt idStat, UInt idSat, Double &ePe, Double &redundancy, UInt &obsCount, UInt &outlierCount)
{
  try
  {
    const UInt idSatStart = ((idSat!=NULLINDEX) ? idSat : 0);
    const UInt idSatEnd   = ((idSat!=NULLINDEX) ? idSat : MAX_UINT);
    for(auto station : slr->stations)
      if(normalEquationInfo.estimateStation.at(station->idStat()) && ((idStat == NULLINDEX) || (idStat == station->idStat())))
        for(UInt idSat=idSatStart; (idSat<=idSatEnd) && (idSat<station->observations.size()); idSat++)
          if(normalEquationInfo.estimateSatellite.at(idSat))
            for(auto &obs : station->observations.at(idSat))
              for(UInt i=0; i<obs->residuals.rows(); i++)
              {
                ePe        += std::pow(obs->residuals(i)/obs->sigmas(i), 2);
                redundancy += obs->redundancies(i);
                obsCount++;
                if(obs->sigmas(i) > obs->sigmas0(i))
                  outlierCount++;
              }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
