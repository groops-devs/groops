/***********************************************/
/**
* @file slrProcessingStepWriteNormalEquations.h
*
* @brief SLR processing step: WriteNormalEquations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPWRITENORMALEQUATIONS__
#define __GROOPS_SLRPROCESSINGSTEPWRITENORMALEQUATIONS__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepWriteNormalEquations = R"(
\subsection{WriteNormalEquations}\label{slrProcessingStepType:writeNormalEquations}
Accumulates the normal equations matrix and writes it.
If \configClass{remainingParameters}{parameterSelectorType}
is set only the selected parameters are written to the normal equations
and all other parameters are eliminated beforehand (implicitly solved).

The solution of the normals would results in $\Delta\M x$
(see \configClass{parametrizations}{slrParametrizationType}). To write the
appropriate apriori vector $\M x_0$ use
\configClass{processingStep:writeAprioriSolution}{slrProcessingStepType:writeAprioriSolution}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileNormalEquation.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "slr/slrNormalEquationInfo.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: WriteNormalEquations.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepWriteNormalEquations : public SlrProcessingStepBase
{
  FileName             fileNameNormals;
  ParameterSelectorPtr parameterSelector;
  Bool                 constraintsOnly;
  UInt                 defaultBlockSize;

public:
  SlrProcessingStepWriteNormalEquations(Config &config);
  void process(SlrProcessingStep::State &state) override;
};

/***********************************************/

inline SlrProcessingStepWriteNormalEquations::SlrProcessingStepWriteNormalEquations(Config &config)
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

inline void SlrProcessingStepWriteNormalEquations::process(SlrProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== write normal equations =================================="<<Log::endl;
    std::vector<UInt> blockIndex, indexVector(state.normalEquationInfo.parameterCount());
    std::iota(indexVector.begin(), indexVector.end(), 0);

    // compute index vectors and block structure for remaining parameters
    UInt eliminationBlocks = 0, eliminationCount = 0;
    if(parameterSelector)
    {
      indexVector = parameterSelector->indexVector(state.normalEquationInfo.parameterNames());
      std::vector<UInt> eliminationIndexVector = ParameterSelector::indexVectorComplement(indexVector, state.normalEquationInfo.parameterCount());
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
    state.buildNormals(constraintsOnly);
    // merge right hand side to one matrix
    Matrix rhs(state.normals.parameterCount(), state.lPl.rows());
    for(UInt i=0; i<state.normals.blockCount(); i++)
      copy(state.n.at(i), rhs.row(state.normals.blockIndex(i), state.normals.blockSize(i)));
    Vector x0 = state.slr->aprioriParameter(state.normalEquationInfo);

    // reorder normals
    if((defaultBlockSize != NULLINDEX) || parameterSelector)
    {
      state.normals.reorder(indexVector, blockIndex);
      rhs = reorder(rhs, indexVector);
      x0  = reorder(x0,  indexVector);
      if((eliminationCount > 0) && !constraintsOnly)
      {
        state.regularizeNotUsedParameters(0, eliminationBlocks);
        state.normals.cholesky(TRUE, 0, eliminationBlocks, TRUE);
        state.normals.triangularTransSolve(rhs, 0, eliminationBlocks);
        state.obsCount -= eliminationCount;
        for(UInt i=0; i<state.lPl.rows(); i++)
          state.lPl(i) -= quadsum(rhs.slice(0, i, eliminationCount, 1)); // lPl = lPl - n1' N1^(-1) n1
      }
      state.normals.eraseBlocks(0, eliminationBlocks);
      rhs = rhs.row(eliminationCount, rhs.rows()-eliminationCount);
      x0  = x0.row (eliminationCount, x0.rows()-eliminationCount);
    }

    // add a priori solution: n += N*x0
    Matrix Nx(x0.rows(), x0.columns());
    for(UInt i=0; i<state.normals.blockCount(); i++)
    {
      if(state.normals.isMyRank(i,i))
        matMult(1., state.normals.N(i,i), x0.row(state.normals.blockIndex(i), state.normals.blockSize(i)), Nx.row(state.normals.blockIndex(i), state.normals.blockSize(i)));
      for(UInt k=i+1; k<state.normals.blockCount(); k++)
        if(state.normals.isMyRank(i,k))
        {
          matMult(1., state.normals.N(i,k),         x0.row(state.normals.blockIndex(k), state.normals.blockSize(k)), Nx.row(state.normals.blockIndex(i), state.normals.blockSize(i)));
          matMult(1., state.normals.N(i,k).trans(), x0.row(state.normals.blockIndex(i), state.normals.blockSize(i)), Nx.row(state.normals.blockIndex(k), state.normals.blockSize(k)));
        }
    }
    for(UInt i=0; i<state.lPl.rows(); i++)
      state.lPl(i) += 2.*inner(x0.column(i), rhs.column(i)) + inner(x0.column(i), Nx.column(i));
    rhs += Nx;

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
