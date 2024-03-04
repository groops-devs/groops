/***********************************************/
/**
* @file gnssProcessingStepWriteNormalEquations.h
*
* @brief GNSS processing step: WriteNormalEquations.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPWRITENORMALEQUATIONS__
#define __GROOPS_GNSSPROCESSINGSTEPWRITENORMALEQUATIONS__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepWriteNormalEquations = R"(
\subsection{WriteNormalEquations}\label{gnssProcessingStepType:writeNormalEquations}
Accumulates the normal equations matrix and writes it.
If \configClass{remainingParameters}{parameterSelectorType}
is set only the selected parameters are written to the normal equations
and all other parameters are eliminated beforehand (implicitly solved).

The solution of the normals would results in $\Delta\M x$
(see \configClass{parametrizations}{gnssParametrizationType}). To write the
appropriate apriori vector $\M x_0$ use
\configClass{processingStep:writeAprioriSolution}{gnssProcessingStepType:writeAprioriSolution}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileNormalEquation.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "gnss/gnssNormalEquationInfo.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: WriteNormalEquations.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepWriteNormalEquations : public GnssProcessingStepBase
{
  FileName             fileNameNormals;
  ParameterSelectorPtr parameterSelector;
  Bool                 constraintsOnly;
  UInt                 defaultBlockSize;

public:
  GnssProcessingStepWriteNormalEquations(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepWriteNormalEquations::GnssProcessingStepWriteNormalEquations(Config &config)
{
  try
  {
    defaultBlockSize = NULLINDEX;
    readConfig(config, "outputfileNormalEquations", fileNameNormals,   Config::MUSTSET,  "output/normals_{loopTime:%D}.dat", "normals");
    readConfig(config, "remainingParameters",       parameterSelector, Config::OPTIONAL, "",  "parameter order/selection of output normal equations");
    readConfig(config, "constraintsOnly",           constraintsOnly,   Config::DEFAULT,  "0", "write only normals of constraints without observations");
    readConfig(config, "defaultNormalsBlockSize",   defaultBlockSize,  Config::OPTIONAL, "",  "block size for distributing the normal equations, 0: one block, empty: original block size");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssProcessingStepWriteNormalEquations::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== write normal equations =================================="<<Log::endl;
    if(!state.normalEquationInfo.parameterCount())
    {
      logWarningOnce<<"Empty normal equations matrix"<<Log::endl;
      return;
    }
    Bool eliminateEpochParameters = FALSE;
    std::vector<UInt> blockIndex, indexVector(state.normalEquationInfo.parameterCount());
    std::iota(indexVector.begin(), indexVector.end(), 0);

    // compute index vectors and block structure for remaining parameters
    UInt eliminationBlocks = 0, eliminationCount = 0;
    if(parameterSelector)
    {
      indexVector = parameterSelector->indexVector(state.normalEquationInfo.parameterNames());

      const UInt countEpochParameter = state.normalEquationInfo.blockIndex(state.normalEquationInfo.blockInterval());
      eliminateEpochParameters = (*std::min_element(indexVector.begin(), indexVector.end()) >= countEpochParameter);
      std::vector<UInt> eliminationIndexVector = ParameterSelector::indexVectorComplement(indexVector, state.normalEquationInfo.parameterCount());
      if(eliminateEpochParameters)
        eliminationIndexVector.erase(std::remove_if(eliminationIndexVector.begin(), eliminationIndexVector.end(),
                                                    [&](UInt x){return x < countEpochParameter;}), eliminationIndexVector.end());
      eliminationCount = eliminationIndexVector.size();
      indexVector.insert(indexVector.begin(), eliminationIndexVector.begin(), eliminationIndexVector.end());

      blockIndex = MatrixDistributed::computeBlockIndex(eliminationCount, 2048);
      eliminationBlocks = blockIndex.size()-1;
      if(defaultBlockSize == 0)
        blockIndex.push_back(indexVector.size());
      else
        while(blockIndex.back() < indexVector.size())
          blockIndex.push_back(std::min(blockIndex.back()+((defaultBlockSize == NULLINDEX) ? 2048 : defaultBlockSize), indexVector.size()));
    }
    else if(defaultBlockSize != NULLINDEX)
      blockIndex = MatrixDistributed::computeBlockIndex(state.normalEquationInfo.parameterCount(), defaultBlockSize);

    // accumulate normals
    state.buildNormals(constraintsOnly, eliminateEpochParameters);
    // merge right hand side to one matrix
    Matrix rhs(state.normals.parameterCount(), state.lPl.rows());
    if(Parallel::isMaster(state.normalEquationInfo.comm))
      for(UInt i=0; i<state.normals.blockCount(); i++)
        copy(state.n.at(i), rhs.row(state.normals.blockIndex(i), state.normals.blockSize(i)));

    // reorder normals
    if((defaultBlockSize != NULLINDEX) || parameterSelector)
    {
      state.normals.reorder(indexVector, blockIndex);
      rhs = reorder(rhs, indexVector);

      if((eliminationCount > 0) && !constraintsOnly)
      {
        std::vector<ParameterName> parameterNames;
        for(UInt i=0; i<eliminationCount; i++)
          parameterNames.push_back((indexVector.at(i) != NULLINDEX) ? state.normalEquationInfo.parameterNames().at(indexVector.at(i)) : ParameterName());
        state.regularizeNotUsedParameters(0, eliminationBlocks, parameterNames);
        state.normals.cholesky(TRUE, 0, eliminationBlocks, TRUE);
        state.normals.triangularTransSolve(rhs, 0, eliminationBlocks);
        state.obsCount -= eliminationCount;
        for(UInt i=0; i<state.lPl.rows(); i++)
          state.lPl(i) -= quadsum(rhs.slice(0, i, eliminationCount, 1)); // lPl = lPl - n1' N1^(-1) n1
      }
      state.normals.eraseBlocks(0, eliminationBlocks);
      rhs = rhs.row(eliminationCount, rhs.rows()-eliminationCount);
    }

    // parameter names
    std::vector<ParameterName> parameterNames;
    for(UInt i=eliminationCount; i<indexVector.size(); i++)
      parameterNames.push_back((indexVector.at(i) != NULLINDEX) ? state.normalEquationInfo.parameterNames().at(indexVector.at(i)) : ParameterName());

    logStatus<<"write normal equations to <"<<fileNameNormals<<">"<<Log::endl;
    writeFileNormalEquation(fileNameNormals, NormalEquationInfo(parameterNames, state.lPl, state.obsCount), state.normals, rhs);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
