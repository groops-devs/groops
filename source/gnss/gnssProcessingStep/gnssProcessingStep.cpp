/***********************************************/
/**
* @file gnssProcessingStep.cpp
*
* @brief Processing steps for GNSS normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#define DOCSTRING_GnssProcessingStep

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "config/configRegister.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "misc/varianceComponentEstimation.h"
#include "gnss/gnss.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepEstimate.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepResolveAmbiguities.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepComputeCovarianceMatrix.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepWriteResults.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepWriteNormalEquations.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepWriteAprioriSolution.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepWriteResiduals.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepWriteUsedStationList.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepWriteUsedTransmitterList.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepPrintResidualStatistics.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepSelectParametrizations.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepSelectEpochs.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepSelectNormalsBlockStructure.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepSelectReceivers.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepForEachReceiverSeparately.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepGroup.h"
#include "gnss/gnssProcessingStep/gnssProcessingStepDisableTransmitterShadowEpochs.h"

/***********************************************/

GROOPS_REGISTER_CLASS(GnssProcessingStep, "gnssProcessingStepType",
                      GnssProcessingStepEstimate,
                      GnssProcessingStepResolveAmbiguities,
                      GnssProcessingStepComputeCovarianceMatrix,
                      GnssProcessingStepWriteResults,
                      GnssProcessingStepWriteNormalEquations,
                      GnssProcessingStepWriteAprioriSolution,
                      GnssProcessingStepWriteResiduals,
                      GnssProcessingStepWriteUsedStationList,
                      GnssProcessingStepWriteUsedTransmitterList,
                      GnssProcessingStepPrintResidualStatistics,
                      GnssProcessingStepSelectParametrizations,
                      GnssProcessingStepSelectEpochs,
                      GnssProcessingStepSelectNormalsBlockStructure,
                      GnssProcessingStepSelectReceivers,
                      GnssProcessingStepForEachReceiverSeparately,
                      GnssProcessingStepGroup,
                      GnssProcessingStepDisableTransmitterShadowEpochs)

GROOPS_READCONFIG_UNBOUNDED_CLASS(GnssProcessingStep, "gnssProcessingStepType")

/***********************************************/

GnssProcessingStep::GnssProcessingStep(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", ""))
    {
      if(readConfigChoiceElement(config, "estimate",                       type, "least squares adjustment"))
        bases.push_back(new GnssProcessingStepEstimate(config));
      if(readConfigChoiceElement(config, "resolveAmbiguities",             type, "resolve integer ambiguities"))
        bases.push_back(new GnssProcessingStepResolveAmbiguities(config));
      if(readConfigChoiceElement(config, "computeCovarianceMatrix",        type, "compute covariance matrix"))
        bases.push_back(new GnssProcessingStepComputeCovarianceMatrix(config));
      if(readConfigChoiceElement(config, "writeResults",                   type, "write all estimated parameters"))
        bases.push_back(new GnssProcessingStepWriteResults(config));
      if(readConfigChoiceElement(config, "writeNormalEquations",           type, "write unconstrained and constraint normal equations"))
        bases.push_back(new GnssProcessingStepWriteNormalEquations(config));
      if(readConfigChoiceElement(config, "writeAprioriSolution",           type, "write apriori solution vector"))
        bases.push_back(new GnssProcessingStepWriteAprioriSolution(config));
      if(readConfigChoiceElement(config, "writeResiduals",                 type, "write observation residuals"))
        bases.push_back(new GnssProcessingStepWriteResiduals(config));
      if(readConfigChoiceElement(config, "writeUsedStationList",           type, "write used stations"))
        bases.push_back(new GnssProcessingStepWriteUsedStationList(config));
      if(readConfigChoiceElement(config, "writeUsedTransmitterList",       type, "write used transmitters"))
        bases.push_back(new GnssProcessingStepWriteUsedTransmitterList(config));
      if(readConfigChoiceElement(config, "printResidualStatistics",        type, "print residual statistics"))
        bases.push_back(new GnssProcessingStepPrintResidualStatistics(config));
      if(readConfigChoiceElement(config, "selectParametrizations",         type, "select parametrizations for all subsequent processing steps"))
        bases.push_back(new GnssProcessingStepSelectParametrizations(config));
      if(readConfigChoiceElement(config, "selectEpochs",                   type, "select epochs to be used in all subsequent processing steps"))
        bases.push_back(new GnssProcessingStepSelectEpochs(config));
      if(readConfigChoiceElement(config, "selectNormalsBlockStructure",    type, "select block structure of the distributed normal equation matrix for all subsequent processing steps"))
        bases.push_back(new GnssProcessingStepSelectNormalsBlockStructure(config));
      if(readConfigChoiceElement(config, "selectReceivers",                type, "use this subset of stations in all subsequent processing steps"))
        bases.push_back(new GnssProcessingStepSelectReceivers(config));
      if(readConfigChoiceElement(config, "forEachReceiverSeparately",      type, "process receiver by receiver, transmitter parameters disabled"))
        bases.push_back(new GnssProcessingStepForEachReceiverSeparately(config));
      if(readConfigChoiceElement(config, "group",                          type, "group processing steps"))
        bases.push_back(new GnssProcessingStepGroup(config));
      if(readConfigChoiceElement(config, "disableTransmitterShadowEpochs", type, "disable transmitter epochs in eclipse"))
        bases.push_back(new GnssProcessingStepDisableTransmitterShadowEpochs(config));
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

GnssProcessingStep::~GnssProcessingStep()
{
  for(auto base : bases)
    delete base;
}

/***********************************************/

void GnssProcessingStep::process(GnssProcessingStep::State &state)
{
  try
  {
    for(auto base : bases)
    {
      Parallel::peek(state.normalEquationInfo.comm);
      if(base->expectInitializedParameters() && state.changedNormalEquationInfo)
        state.gnss->initParameter(state.normalEquationInfo);
      state.changedNormalEquationInfo = FALSE;
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

GnssProcessingStep::State::State(GnssPtr gnss, Parallel::CommunicatorPtr comm) :
  gnss(gnss), normalEquationInfo(gnss->times.size(), gnss->receivers.size(), gnss->transmitters.size(), comm),
  changedNormalEquationInfo(TRUE), sigmaType(gnss->receivers.size()), sigmaFactor(gnss->receivers.size())
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

void GnssProcessingStep::State::regularizeNotUsedParameters(UInt blockStart, UInt blockCount)
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
    Parallel::reduceSum(countZeros, 0, normals.communicator());
    if(countZeros && Parallel::isMaster(normals.communicator()))
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

void GnssProcessingStep::State::collectNormalsBlocks(UInt blockStart, UInt blockCount)
{
  try
  {
    if(Parallel::size(normals.communicator()) <= 1)
      return;

    auto countBlocks = [&]()
    {
      UInt countBlocks = 0;
      for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
        normals.loopBlockRow(idBlock, {idBlock, normals.blockCount()}, [&](UInt /*k*/, UInt /*ik*/) {countBlocks++;});
      return countBlocks;
    };

    // synchronize used blocks
    for(UInt idProcess=0; idProcess<Parallel::size(normals.communicator()); idProcess++)
    {
      Matrix rowColRank(countBlocks(), 3);
      UInt i = 0;
      for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
        normals.loopBlockRow(idBlock, {idBlock, normals.blockCount()}, [&](UInt k, UInt ik)
        {
          rowColRank(i, 0) = idBlock;
          rowColRank(i, 1) = k;
          rowColRank(i, 2) = normals._rank[ik];
          i++;
        });

      Parallel::broadCast(rowColRank, idProcess, normals.communicator());
      for(UInt i=0; i<rowColRank.rows(); i++)
        normals.setBlock(static_cast<UInt>(rowColRank(i,0)), static_cast<UInt>(rowColRank(i,1)), static_cast<UInt>(rowColRank(i,2)));
    }

    // for each block x which process contains
    Matrix usedRanks(countBlocks(), Parallel::size(normals.communicator()));
    UInt i = 0;
    for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
      normals.loopBlockRow(idBlock, {idBlock, normals.blockCount()}, [&](UInt /*k*/, UInt ik)
      {
        if(normals._N[ik].size())
          usedRanks(i, Parallel::myRank(normals.communicator())) = 1.;
        i++;
      });
    Parallel::reduceSum(usedRanks, 0, normals.communicator());
    Parallel::broadCast(usedRanks, 0, normals.communicator());

    // reduce sum normals
    i = 0;
    for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
      normals.loopBlockRow(idBlock, {idBlock, normals.blockCount()}, [&](UInt /*k*/, UInt ik)
      {
        std::vector<Bool> usedRank(usedRanks.columns());
        for(UInt idProcess=0; idProcess<usedRank.size(); idProcess++)
          usedRank[idProcess] = usedRanks(i, idProcess);
        normals.reduceSum(normals._N[ik], ik, usedRank, TRUE/*free*/);
        i++;
      });

    // right hand side
    for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
    {
      Parallel::reduceSum(n.at(idBlock), 0, normals.comm);
      if(!Parallel::isMaster(normals.comm))
        n.at(idBlock).setNull();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessingStep::State::buildNormals(Bool constraintsOnly, Bool solveEpochParameters)
{
  try
  {
    normals.initEmpty(normalEquationInfo.blockIndices(), normalEquationInfo.comm);
    n.resize(normals.blockCount());
    for(UInt i=0; i<normals.blockCount(); i++)
      n.at(i) = Vector(normals.blockSize(i));
    lPl = Vector(1);
    obsCount = 0;

    // Loop over all epochs
    // --------------------
    Parallel::barrier(normalEquationInfo.comm);
    logStatus<<"accumulate normals"<<Log::endl;
    GnssDesignMatrix A(normalEquationInfo);
    std::vector<GnssObservationEquation> eqns(gnss->transmitters.size());
    UInt blockStart = 0; // first block, which is not regularized and reduced
    UInt blockCount = 0;
    UInt idLoop     = 0;
    logTimerStart;
    for(UInt idEpoch : normalEquationInfo.idEpochs)
    {
      logTimerLoop(idLoop++, normalEquationInfo.idEpochs.size());

      gnss->constraintsEpoch(normalEquationInfo, idEpoch, normals, n, lPl(0), obsCount);

      // loop over all receivers
      if(!constraintsOnly)
        for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
          if(normalEquationInfo.estimateReceiver.at(idRecv) && gnss->receivers.at(idRecv)->isMyRank())
          {
            // all observation equations for this epoch
            UInt countEqn = 0;
            for(UInt idTrans=0; idTrans<gnss->receivers.at(idRecv)->idTransmitterSize(idEpoch); idTrans++)
              if(gnss->basicObservationEquations(normalEquationInfo, idRecv, idTrans, idEpoch, eqns.at(countEqn)))
              {
                if(eqns.at(countEqn).B.size()) // eliminate STEC
                  eliminationParameter(Matrix(eqns.at(countEqn).B), {eqns.at(countEqn).A, eqns.at(countEqn).l});
                countEqn++;
              }

            if(!normalEquationInfo.accumulateEpochObservations)
            {
              for(UInt i=0; i<countEqn; i++)
              {
                A.init(eqns.at(i).l);
                gnss->designMatrix(normalEquationInfo, eqns.at(i), A);
                A.accumulateNormals(normals, n, lPl(0), obsCount);
              }
            }
            else
            {
              // copy all observations to a single vector
              Vector l(std::accumulate(eqns.begin(), eqns.end(), UInt(0), [](UInt count, auto &e) {return count+e.l.rows();}));
              UInt idx=0;
              for(UInt i=0; i<countEqn; i++)
              {
                copy(eqns.at(i).l, l.row(idx, eqns.at(i).l.rows()));
                idx += eqns.at(i).l.rows();
              }
              A.init(l);

              idx=0;
              for(UInt i=0; i<countEqn; i++)
              {
                gnss->designMatrix(normalEquationInfo, eqns.at(i), A.selectRows(idx, eqns.at(i).l.rows()));
                idx += eqns.at(i).l.rows();
              }

              A.selectRows(0, 0); // select all
              A.accumulateNormals(normals, n, lPl(0), obsCount);
            }
         } // for(idRecv)

      // perform following steps not every epoch
      blockCount += normalEquationInfo.blockCountEpoch(idEpoch);
      if((blockCount < normalEquationInfo.defaultBlockCountReduction) && (idEpoch != normalEquationInfo.idEpochs.back()))
        continue;

      collectNormalsBlocks(blockStart, blockCount);

      if(solveEpochParameters && !constraintsOnly)
      {
        regularizeNotUsedParameters(blockStart, blockCount);
        normals.cholesky(FALSE, blockStart, blockCount, FALSE);
        normals.triangularTransSolve(n, blockStart, blockCount, FALSE);
        if(Parallel::isMaster(normalEquationInfo.comm))
          for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
          {
            lPl(0)   -= quadsum(n.at(idBlock)); // lPl = lPl - n1' N1^(-1) n1
            obsCount -= normals.blockSize(idBlock);
          }
      }

      if(solveEpochParameters)
      {
        // remove N12 (epoch <-> other parameters)
        for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
          normals.loopBlockRow(idBlock, {normalEquationInfo.blockInterval(), normals.blockCount()}, [&](UInt /*k*/, UInt ik)
          {
            if(normals._N[ik].size())
              normals._N[ik] = Matrix();
          });
      }

      blockStart += blockCount;
      blockCount  = 0;
    } // for(idEpoch)
    Parallel::barrier(normalEquationInfo.comm);
    logTimerLoopEnd(normalEquationInfo.idEpochs.size());

    // other observations and constraints
    // ----------------------------------
    gnss->constraints(normalEquationInfo, normals, n, lPl(0), obsCount);

    // collect normal equations
    // ------------------------
    if(Parallel::size(normals.comm)>1)
    {
      logStatus<<"collect normal equations"<<Log::endl;
      collectNormalsBlocks(blockStart, normals.blockCount()-blockStart);
      Parallel::reduceSum(lPl,      0, normalEquationInfo.comm);
      Parallel::reduceSum(obsCount, 0, normalEquationInfo.comm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssProcessingStep::State::estimateSolution(const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger,
                                                   Bool computeResiduals, Bool computeWeights, Bool adjustSigma0, Double huber, Double huberPower)
{
  try
  {
    // setup normal equations
    const Bool solveEpochParameters = !normalEquationInfo.keepEpochNormalsInMemory && (normalEquationInfo.blockInterval() < normalEquationInfo.blockCount());
    buildNormals(FALSE/*constraintsOnly*/, solveEpochParameters);

    // eliminate all other parameters from the ambiguity normals
    // ---------------------------------------------------------
    Parallel::barrier(normalEquationInfo.comm);
    logStatus<<"solve"<<Log::endl;
    const UInt blockStart = (solveEpochParameters) ? normalEquationInfo.blockInterval() : 0; // epoch parameters already eliminated?
    const UInt blockCount = ((searchInteger) ? normalEquationInfo.blockAmbiguity() : normals.blockCount()) - blockStart;
    regularizeNotUsedParameters(blockStart, normals.blockCount()-blockStart);
    normals.cholesky(TRUE/*timing*/, blockStart, blockCount, TRUE/*collect*/);
    normals.triangularTransSolve(n, blockStart, blockCount, TRUE/*collect*/); // forward step
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(UInt z=blockStart; z<blockStart+blockCount; z++)
      {
        lPl(0)   -= quadsum(n.at(z)); // lPl = lPl - n1' N1^(-1) n1
        obsCount -= normals.blockSize(z);
      }

    // resolve integer ambiguities
    // ---------------------------
    if(searchInteger)
    {
      logStatus<<"Resolve integer ambiguities (may take a while)"<<Log::endl;
      Double sigmaFloat = gnss->ambiguityResolve(normalEquationInfo, normals, n, lPl(0), obsCount, searchInteger);
      logInfo<<"  sigma(float) = "<<sigmaFloat%"%.2f"s<<Log::endl;
      logInfo<<"  sigma(fixed) = "<<std::sqrt(lPl(0)/obsCount)%"%.2f"s<<Log::endl;
      Parallel::barrier(normalEquationInfo.comm);
    } // if(parameterCountAmbiguities)

    // Compute sigma
    // -------------
    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      const Double sigma = Vce::standardDeviation(lPl(0), obsCount, huber, huberPower);
      logInfo<<"  sigma = "<<sigma%"%.2f"s<<Log::endl;
      if((sigma!=sigma) || (sigma<=0))
        logWarning<<"  Cannot compute sigma = sqrt("<<lPl(0)<<"/"<<obsCount<<")"<<Log::endl;
    }

    constexpr UInt monteCarloColumns = 100;
    std::vector<Matrix> monteCarlo; // Monte Carlo Vector for redundancy computation

    // ========================================================================================

    // solve (Backward step)
    // ---------------------
    if(!solveEpochParameters)
    {
      normals.triangularSolve(n);
      if(computeResiduals)
      {
        monteCarlo.resize(normals.blockCount());
        for(UInt i=0; i<normals.blockCount(); i++)
          monteCarlo.at(i) = Matrix(normals.blockSize(i), monteCarloColumns);
        if(Parallel::isMaster(normalEquationInfo.comm))
          for(UInt i=0; i<blockCount; i++) // possible without resolved ambiguity parameters
            monteCarlo.at(i) = Vce::monteCarlo(monteCarlo.at(i).rows(), monteCarlo.at(i).columns());
        normals.triangularSolve(monteCarlo, 0, blockCount);
      }
    }

    // ========================================================================================

    // Reconstruct epoch parameters
    // ----------------------------
    if(solveEpochParameters)
    {
      logStatus<<"Reconstruct epoch parameters"<<Log::endl;

      // solve other parameters (Backward step)
      // --------------------------------------
      normals.triangularSolve(n, blockStart, normals.blockCount()-blockStart);

      if(computeResiduals)
      {
        monteCarlo.resize(normals.blockCount());
        for(UInt i=blockStart; i<normals.blockCount(); i++)
          monteCarlo.at(i) = Matrix(normals.blockSize(i), monteCarloColumns);
        if(Parallel::isMaster(normalEquationInfo.comm))
          for(UInt i=blockStart; i<blockStart+blockCount; i++) // possible without resolved ambiguity parameters
            monteCarlo.at(i) = Vce::monteCarlo(monteCarlo.at(i).rows(), monteCarlo.at(i).columns());
        normals.triangularSolve(monteCarlo, blockStart, blockCount);
      }

      // free N12, N22
      normals.eraseBlocks(blockStart, normals.blockCount()-blockStart);

      // broadcast n, Wz
      // ---------------
      for(UInt i=0; i<blockStart; i++)
        n.at(i).setNull();
      for(UInt i=blockStart; i<normalEquationInfo.blockCount(); i++)
        Parallel::broadCast(n.at(i), 0, normalEquationInfo.comm);
      if(computeResiduals)
      {
        for(UInt i=0; i<blockStart; i++)
          monteCarlo.at(i) = Matrix(normalEquationInfo.blockSize(i), monteCarloColumns);
        for(UInt i=blockStart; i<blockStart+blockCount; i++)
          Parallel::broadCast(monteCarlo.at(i), 0, normalEquationInfo.comm);
      }

      // reconstruct N12
      // ---------------
      GnssDesignMatrix A(normalEquationInfo);
      UInt idLoop = 0;
      logTimerStart;
      for(UInt idEpoch : normalEquationInfo.idEpochs)
      {
        logTimerLoop(idLoop++, normalEquationInfo.idEpochs.size());

        // loop over all receivers
        GnssObservationEquation eqn;
        for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
          if(normalEquationInfo.estimateReceiver.at(idRecv) && gnss->receivers.at(idRecv)->isMyRank())
            for(UInt idTrans=0; idTrans<gnss->receivers.at(idRecv)->idTransmitterSize(idEpoch); idTrans++)
              if(gnss->basicObservationEquations(normalEquationInfo, idRecv, idTrans, idEpoch, eqn))
              {
                if(eqn.B.size()) // eliminate STEC
                  eliminationParameter(Matrix(eqn.B), {eqn.A, eqn.l});
                A.init(eqn.l);
                gnss->designMatrix(normalEquationInfo, eqn, A);
                A.transMult(eqn.l-A.mult(n, blockStart, normalEquationInfo.blockCount()-blockStart), n, 0, blockStart);
                A.transMult(-A.mult(monteCarlo, blockStart, blockCount), monteCarlo, 0, blockStart);
              }
      } // for(idEpoch)
      Parallel::barrier(normalEquationInfo.comm);
      logTimerLoopEnd(normalEquationInfo.idEpochs.size());

      // collect
      // -------
      auto collect = [&](std::vector<Matrix> &n, UInt columns)
      {
        Matrix tmp(normalEquationInfo.blockIndex(blockStart), columns);
        for(UInt i=0; i<blockStart; i++) // copy to one block (better performance of reduceSum)
          if(normalEquationInfo.blockSize(i))
            copy(n.at(i), tmp.row(normalEquationInfo.blockIndex(i), normalEquationInfo.blockSize(i)));
        Parallel::reduceSum(tmp, 0, normalEquationInfo.comm);
        if(Parallel::isMaster(normalEquationInfo.comm))
          for(UInt i=0; i<blockStart; i++)
            if(normalEquationInfo.blockSize(i))
              copy(tmp.row(normalEquationInfo.blockIndex(i), normalEquationInfo.blockSize(i)), n.at(i));
        if(!Parallel::isMaster(normalEquationInfo.comm))
          for(UInt i=0; i<blockStart; i++)
            n.at(i).setNull();
      };

      collect(n, 1);
      if(computeResiduals)
        collect(monteCarlo, monteCarloColumns);

      // solve epoch (Forward step)
      // --------------------------
      normals.triangularTransSolve(n, 0, blockStart, FALSE);
      if(computeResiduals)
      {
        normals.triangularTransSolve(monteCarlo, 0, blockStart, FALSE);
        if(Parallel::isMaster(normalEquationInfo.comm))
          for(UInt i=0; i<blockStart; i++)
            axpy(1, Vce::monteCarlo(monteCarlo.at(i).rows(), monteCarlo.at(i).columns()), monteCarlo.at(i));
      }

      // solve epoch (Backward step)
      // ---------------------------
      normals.triangularSolve(n, 0, blockStart);
      if(computeResiduals)
        normals.triangularSolve(monteCarlo, 0, blockStart);
    }

    // ========================================================================================

    // free memory
    normals.initEmpty(normals.blockIndex(), normals.communicator());

    // copy solution to one vector
    // ---------------------------
    Vector x;
    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      x = Vector(normalEquationInfo.parameterCount());
      for(UInt i=0; i<normalEquationInfo.blockCount(); i++)
        if(normalEquationInfo.blockSize(i))
          copy(n.at(i), x.row(normalEquationInfo.blockIndex(i), normalEquationInfo.blockSize(i)));
    }
    n.clear();
    n.shrink_to_fit();
    Parallel::broadCast(x, 0, normalEquationInfo.comm);
    if(!computeResiduals)
    {
      logInfo<<"Parameter changes"<<Log::endl;
      return gnss->updateParameter(normalEquationInfo, x, Matrix());
    }

    // copy monte carlo to one vector
    // ------------------------------
    Matrix Wz;
    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      Wz = Matrix(normalEquationInfo.parameterCount(), monteCarloColumns);
      for(UInt i=0; i<normalEquationInfo.blockCount(); i++)
        if(normalEquationInfo.blockSize(i))
          copy(monteCarlo.at(i), Wz.row(normalEquationInfo.blockIndex(i), normalEquationInfo.blockSize(i)));
    }
    monteCarlo.clear();
    monteCarlo.shrink_to_fit();
    Parallel::broadCast(Wz, 0, normalEquationInfo.comm);

    // Residual tracking
    // -----------------
    Gnss::InfoParameterChange infoTec("tec");
    std::vector<GnssType>     typesResiduals = gnss->types(~(GnssType::PRN + GnssType::FREQ_NO));
    std::vector<Gnss::InfoParameterChange> infosResiduals(typesResiduals.size(), Gnss::InfoParameterChange("mm"));

    Double minSTEC   =  std::numeric_limits<Double>::infinity();
    Double maxSTEC   = -std::numeric_limits<Double>::infinity();
    Double meanSTEC  = 0;
    Double stdSTEC   = 0;
    UInt   countSTEC = 0;

    Parallel::barrier(normalEquationInfo.comm);
    logStatus<<"Compute residuals"<<Log::endl;
    GnssDesignMatrix A(normalEquationInfo);
    UInt idLoop = 0;
    logTimerStart;
    for(UInt idEpoch : normalEquationInfo.idEpochs)
    {
      logTimerLoop(idLoop++, normalEquationInfo.idEpochs.size());

      // loop over all receivers
      for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
        if(normalEquationInfo.estimateReceiver.at(idRecv) && gnss->receivers.at(idRecv)->isMyRank())
        {
          GnssObservationEquation eqn;
          for(UInt idTrans=0; idTrans<gnss->receivers.at(idRecv)->idTransmitterSize(idEpoch); idTrans++)
            if(gnss->basicObservationEquations(normalEquationInfo, idRecv, idTrans, idEpoch, eqn))
            {
              // setup observation equations
              A.init(eqn.l);
              gnss->designMatrix(normalEquationInfo, eqn, A);
              Vector We  = eqn.l - A.mult(x); // decorrelated residuals
              Matrix AWz = A.mult(Wz);        // redundancies

              // estimate & reduce B parameters (ionosphere)
              // -------------------------------------------
              if(eqn.B.size())
              {
                // estimate and eliminate STEC parameter
                Matrix B = eqn.B;
                Vector tau = QR_decomposition(B);
                QTransMult(B, tau, We);
                triangularSolve(1., B.row(0, B.columns()), We.row(0, B.columns()));
                eqn.receiver->observation(eqn.transmitter->idTrans(), eqn.idEpoch)->updateParameter(We.row(0,B.columns())); // ionosphere parameter

                if((norm(eqn.sigma-eqn.sigma0) < 1e-8) && infoTec.update(We(0)))
                  infoTec.info = "STEC ("+eqn.receiver->name()+", "+eqn.transmitter->name()+", "+eqn.timeRecv.dateTimeStr()+")";

                // STEC statistics
                const Double STEC = gnss->receivers.at(idRecv)->observation(idTrans, idEpoch)->dSTEC;
                minSTEC   = std::min(STEC, minSTEC);
                maxSTEC   = std::max(STEC, maxSTEC);
                meanSTEC += STEC;
                stdSTEC  += STEC*STEC;
                countSTEC++;

                // remove STEC from residuals
                We.row(0, B.columns()).setNull();
                QMult(B, tau, We);

                // influence of B parameters = B(B'B)^(-1)B'
                QTransMult(B, tau, AWz);
                AWz.row(0, B.columns()).setNull();
                for(UInt i=0; i<B.columns(); i++)
                  AWz(i,i) = 1.0;
                QMult(B, tau, AWz);
              }

              // redundancies
              // ------------
              Vector r(We.rows());
              for(UInt i=0; i<We.rows(); i++)
                r(i) = 1. - quadsum(AWz.row(i));

              // find max. residual (for statistics)
              // -----------------------------------
              if(norm(eqn.sigma-eqn.sigma0) < 1e-8) // without outlier
                for(UInt k=0; k<eqn.types.size(); k++)
                {
                  const UInt idType = GnssType::index(typesResiduals, eqn.types.at(k));
                  if(infosResiduals.at(idType).update(1e3*(We(k)*eqn.sigma(k) - gnss->receivers.at(idRecv)->observation(idTrans, idEpoch)->at(eqn.types.at(k)).residuals)))
                    infosResiduals.at(idType).info = (typesResiduals.at(idType)+eqn.transmitter->PRN()).str()+", ("+ eqn.receiver->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
                } // for(k)

              gnss->receivers.at(idRecv)->observation(idTrans, idEpoch)->setDecorrelatedResiduals(eqn.types, We, r);
            } // for(idTrans)
        } // for(idRecv)
    } // for(idEpoch)
    Parallel::barrier(normalEquationInfo.comm);
    logTimerLoopEnd(normalEquationInfo.idEpochs.size());

    // new weights
    // -----------
    if(computeWeights || adjustSigma0)
    {
      if(computeWeights) logStatus<<"Downweight outliers"<<Log::endl;
      if(adjustSigma0)   logStatus<<"Estimate variance factors"<<Log::endl;
      for(auto recv : gnss->receivers)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->isMyRank())
          for(UInt iter=0; iter<10; iter++)
          {
            // --- lambda -------------------------------------
            auto getResiduals = [&](UInt idEpoch, UInt idTrans, Bool useSigma0, std::vector<GnssType> &types, std::vector<Double> &ePe, std::vector<Double> &redundancy)
            {
              GnssObservation &obs = *(recv->observation(idTrans, idEpoch));
              for(UInt idType=0; idType<obs.size(); idType++)
                if((obs.at(idType).sigma0 > 0) && (obs.at(idType).redundancy > 0))
                {
                  UInt idx = GnssType::index(types, obs.at(idType).type);
                  if(idx == NULLINDEX) // new type?
                  {
                    idx = types.size();
                    if(obs.at(idType).type == GnssType::PHASE)
                      types.push_back(obs.at(idType).type & (GnssType::TYPE + GnssType::SYSTEM)); // only PHASE and System
                    else
                      types.push_back(obs.at(idType).type & (GnssType::TYPE + GnssType::FREQUENCY + GnssType::ATTRIBUTE + GnssType::SYSTEM));
                    ePe.push_back(0);
                    redundancy.push_back(0);
                  }
                  ePe.at(idx)        += std::pow(obs.at(idType).residuals/(useSigma0 ? obs.at(idType).sigma0 : obs.at(idType).sigma), 2);
                  redundancy.at(idx) += obs.at(idType).redundancy;
                }
            };
            // ------------------------------------------------

            if(adjustSigma0)
            {
              std::vector<GnssType> types;
              std::vector<Double>   ePe, redundancy;
              for(UInt idEpoch : normalEquationInfo.idEpochs)
                if(recv->useable(idEpoch))
                  for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
                    if(recv->observation(idTrans, idEpoch))
                      getResiduals(idEpoch, idTrans, FALSE/*useSigma0*/, types, ePe, redundancy);

              std::vector<Double> factors(types.size(), 1.);
              for(UInt idType=0; idType<types.size(); idType++)
                if(redundancy.at(idType) > 3)
                  factors.at(idType) = Vce::standardDeviation(ePe.at(idType), redundancy.at(idType), huber, huberPower);

              // adjust sigma and sigma0
              UInt idx;
              for(UInt idEpoch=0; idEpoch<recv->idEpochSize(); idEpoch++)
                for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
                  if(recv->observation(idTrans, idEpoch))
                  {
                    GnssObservation &obs = *(recv->observation(idTrans, idEpoch));
                    for(UInt idType=0; idType<obs.size(); idType++)
                      if((obs.at(idType).sigma0 > 0) && obs.at(idType).type.isInList(types, idx))
                      {
                        obs.at(idType).sigma0 *= factors.at(idx);
                        obs.at(idType).sigma  *= factors.at(idx);
                      } // for(idType)
                  } // for(idTrans, idEpoch)

              // apply old factors to store total factor
              for(UInt idType=0; idType<types.size(); idType++)
                if(types.at(idType).isInList(sigmaType.at(recv->idRecv()), idx))
                  factors.at(idType) *= sigmaFactor.at(recv->idRecv()).at(idx);
              sigmaType.at(recv->idRecv())   = types;
              sigmaFactor.at(recv->idRecv()) = factors;
            } // if(adjustSigma0)

            if(computeWeights)
            {
              for(UInt idEpoch : normalEquationInfo.idEpochs)
                if(recv->useable(idEpoch))
                  for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
                    if(recv->observation(idTrans, idEpoch))
                    {
                      std::vector<GnssType> types;
                      std::vector<Double> ePe, redundancy;
                      getResiduals(idEpoch, idTrans, TRUE/*useSigma0*/, types, ePe, redundancy);

                      Vector factors(types.size(), 1.);
                      for(UInt i=0; i<types.size(); i++)
                        if(redundancy.at(i) > 0.1)
                        {
                          const Double s = std::sqrt(ePe.at(i)/redundancy.at(i));
                          factors(i) = ((s > huber) ? std::pow(s/huber, huberPower) : 1.);
                        }

                      GnssObservation &obs = *(recv->observation(idTrans, idEpoch));
                      UInt idx;
                      for(UInt idType=0; idType<obs.size(); idType++)
                        if(obs.at(idType).type.isInList(types, idx))
                          obs.at(idType).sigma = obs.at(idType).sigma0 * factors(idx);
                    }
            } // if(computeWeights)

            if(!(computeWeights && adjustSigma0))
              break; // iteration not needed
          } // for(recv, iter)
    }

    // residual analysis
    // -----------------
    std::vector<GnssType> types = typesResiduals;
    std::vector<Double>   ePe(types.size(), 0), redundancy(types.size(), 0);
    std::vector<UInt>     obsCount(types.size(), 0), outlierCount(types.size(), 0);
    residualsStatistics(NULLINDEX/*idRecv*/, NULLINDEX/*idTrans*/, types, ePe, redundancy, obsCount, outlierCount);
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(UInt i=0; i<types.size(); i++)
        if(obsCount.at(i))
        {
          logInfo<<"  "<<types.at(i).str()
                <<": sigma0 = "    <<Vce::standardDeviation(ePe.at(i), redundancy.at(i), huber, huberPower)%"%4.2f"s
                <<", redundancy = "<<(redundancy.at(i)/obsCount.at(i))%"%4.2f"s
                <<", count = "     <<obsCount.at(i)%"%7i"s
                <<", outliers = "  <<outlierCount.at(i)%"%6i"s<<" ("<<(100.*outlierCount.at(i)/obsCount.at(i))%"%4.2f"s<<" %)"
                <<Log::endl;
        }
    logInfo<<"Residuals changes (without outliers)"<<Log::endl;
    Double maxChange = 0;
    for(auto &info : infosResiduals)
      info.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    // STEC analysis
    // -------------
    Parallel::reduceMin(minSTEC,   0, normalEquationInfo.comm);
    Parallel::reduceMax(maxSTEC,   0, normalEquationInfo.comm);
    Parallel::reduceSum(meanSTEC,  0, normalEquationInfo.comm);
    Parallel::reduceSum(stdSTEC,   0, normalEquationInfo.comm);
    Parallel::reduceSum(countSTEC, 0, normalEquationInfo.comm);
    stdSTEC   = std::sqrt((stdSTEC-meanSTEC*meanSTEC/countSTEC)/(countSTEC-1));
    meanSTEC /= countSTEC;
    std::string infoTecStr = " (total: "+meanSTEC%"%.2f +- "s+stdSTEC%"%.2f ["s+minSTEC%"%.2f -- "s+maxSTEC%"%.2f])"s;
    Parallel::broadCast(infoTecStr, 0, normalEquationInfo.comm);
    infoTec.info += infoTecStr;

    logInfo<<"Parameter changes"<<Log::endl;
    infoTec.synchronizeAndPrint(normalEquationInfo.comm, 0, maxChange);
    return gnss->updateParameter(normalEquationInfo, x, Wz);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessingStep::State::residualsStatistics(UInt idRecv, UInt idTrans,
                                                    std::vector<GnssType> &types, std::vector<Double> &ePe, std::vector<Double> &redundancy,
                                                    std::vector<UInt> &obsCount, std::vector<UInt> &outlierCount)
{
  try
  {
    const UInt idTransStart = ((idTrans!=NULLINDEX) ? idTrans : 0);
    const UInt idTransEnd   = ((idTrans!=NULLINDEX) ? idTrans : MAX_UINT);
    for(auto recv : gnss->receivers)
      if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->isMyRank() && ((idRecv == NULLINDEX) || (idRecv == recv->idRecv())))
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(recv->useable(idEpoch))
            for(UInt idTrans=idTransStart; (idTrans<=idTransEnd) && (idTrans<recv->idTransmitterSize(idEpoch)); idTrans++)
              if(recv->observation(idTrans, idEpoch))
              {
                const GnssObservation &obs = *recv->observation(idTrans, idEpoch);
                UInt idx;
                for(UInt idType=0; idType<obs.size(); idType++)
                  if((obs.at(idType).sigma0 > 0) && (obs.at(idType).type.isInList(types, idx)))
                  {
                    ePe.at(idx)        += std::pow(obs.at(idType).residuals/obs.at(idType).sigma, 2);
                    redundancy.at(idx) += obs.at(idType).redundancy;
                    obsCount.at(idx)++;
                    if(obs.at(idType).sigma > obs.at(idType).sigma0)
                      outlierCount.at(idx)++;
                  }
              } // for(idTrans, idEpoch)

    for(UInt i=0; i<types.size(); i++)
    {
      Parallel::reduceSum(ePe.at(i),          0, normalEquationInfo.comm);
      Parallel::reduceSum(redundancy.at(i),   0, normalEquationInfo.comm);
      Parallel::reduceSum(obsCount.at(i),     0, normalEquationInfo.comm);
      Parallel::reduceSum(outlierCount.at(i), 0, normalEquationInfo.comm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
