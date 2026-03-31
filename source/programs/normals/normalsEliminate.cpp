/***********************************************/
/**
* @file normalsEliminate.cpp
*
* @brief Eliminate parameters from a system of normal equations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2003-03-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program eliminates parameters from a system of \configFile{inputfileNormalEquation}{normalEquation}s.
To just remove (cutting out) parameters use \program{NormalsReorder}.

The \configClass{remainingParameters}{parameterSelectorType} allows the selection
of parameters that will remain, all others will be eliminated. The order of remaining parameters
can be modified via the parameter selection. Block size of the output normal matrix can be adjusted with
\config{outBlockSize}. If it is set to zero, the \configFile{outputfileNormalEquation}{normalEquation}
is written to a single block file.

For example the normal equations are divided into two groups of
parameters $\hat{\M x}_1$ and $\hat{\M x}_2$ according to
\begin{equation}
\begin{pmatrix}
  \M N_{11} & \M N_{12} \\
  \M N_{21} & \M N_{22}
\end{pmatrix}
\begin{pmatrix} \hat{\M x}_1 \\ \hat{\M x}_2 \end{pmatrix}
=
\begin{pmatrix}
  \M n_1 \\
  \M n_2
\end{pmatrix}.
\end{equation}
and $\hat{\M x}_2$ shall be eliminated, the reduced system of normal equations is given by
\begin{equation}
\bar{\M N}\hat{\M x} = \bar{\M n}
\qquad\text{with}\qquad
\bar{\M N}=\M N_{11}-\M N_{12}\M N_{22}^{-1}\M N_{12}^T
\qquad\text{and}\qquad\bar{\M n} =  \M n_1 - \M N_{12}\M N_{22}^{-1}\M n_2.
\end{equation}

See also \program{NormalsReorder}.
)";

/***********************************************/

#include "programs/program.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"

/***** CLASS ***********************************/

/** @brief Eliminate parameters from a system of normal equations.
* @ingroup programsGroup */
class NormalsEliminate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsEliminate, PARALLEL, "Eliminate parameters from a system of normal equations.", NormalEquation)

/***********************************************/

void NormalsEliminate::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName outName, inName;
    ParameterSelectorPtr parameterSelector;
    UInt blockSize;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "inputfileNormalequation",  "inputfileNormalEquation",  date2time(2020, 6, 3));

    readConfig(config, "outputfileNormalEquation", outName,            Config::MUSTSET,  "", "");
    readConfig(config, "inputfileNormalEquation",  inName,             Config::MUSTSET,  "", "");
    readConfig(config, "remainingParameters",      parameterSelector,  Config::MUSTSET,  "",     "parameter order/selection of output normal equations");
    readConfig(config, "outBlockSize",             blockSize,          Config::DEFAULT,  "2048", "block size for distributing the normal equations, 0: one block");
    if(isCreateSchema(config)) return;

    // ==================================

    logStatus<<"init normal equations"<<Log::endl;
    MatrixDistributed normal;
    Matrix rhs;
    NormalEquationInfo info;
    readFileNormalEquation(inName, info, normal, rhs, comm);

    // compute index vectors and block structure for remaining parameters
    std::vector<UInt> indexVector = parameterSelector->indexVector(info.parameterName);
    const UInt parameterCountOld  = normal.parameterCount();
    const UInt parameterCountNew  = indexVector.size();
    std::vector<UInt> blockIndex  = MatrixDistributed::computeBlockIndex(parameterCountNew, blockSize);

    // compute index vectors and block structure for to-be-eliminated parameters
    std::vector<UInt> eliminationIndexVector = ParameterSelector::indexVectorComplement(indexVector, parameterCountOld);
    const UInt eliminationCount = eliminationIndexVector.size();
    std::vector<UInt> eliminationBlockIndex  = MatrixDistributed::computeBlockIndex(eliminationCount, blockSize);
    if(eliminationCount == parameterCountOld)
      logWarning<<"eliminating all original parameters"<<Log::endl;

    logInfo<<"  number of unknown parameters (old): "<<parameterCountOld<<Log::endl;
    logInfo<<"  number of unknown parameters (new): "<<parameterCountNew<<Log::endl;
    logInfo<<"  number of right hand sides:         "<<info.lPl.rows()<<Log::endl;

    // create list of remaining parameter names
    std::vector<ParameterName> parameterNamesOld = info.parameterName;
    info.parameterName.resize(parameterCountNew);
    for(UInt i=0; i<parameterCountNew; i++)
      info.parameterName.at(i) = (indexVector.at(i) != NULLINDEX) ? parameterNamesOld.at(indexVector.at(i)) : ParameterName();

    // prepend to-be-eliminated parameters to (remaining) index vector and block structure
    if(eliminationCount > 0)
    {
      for(auto &&index : blockIndex)
        index += eliminationCount;
      eliminationBlockIndex.pop_back();
      blockIndex.insert(blockIndex.begin(), eliminationBlockIndex.begin(), eliminationBlockIndex.end());
      indexVector.insert(indexVector.begin(), eliminationIndexVector.begin(), eliminationIndexVector.end());
    }

    logStatus<<"reorder normal matrix"<<Log::endl;
    normal.reorder(indexVector, blockIndex);
    rhs = reorder(rhs, indexVector);

    if(eliminationCount > 0)
    {
      logStatus<<"eliminate parameters from normal equations"<<Log::endl;
      const UInt eliminationBlocks = eliminationBlockIndex.size();
      for(UInt i=0; i<eliminationBlocks; i++)
      {
        normal.setBlock(i, i);
        if(normal.isMyRank(i,i))
        {
          Matrix &N = normal.N(i,i);
          for(UInt k=0; k<N.rows(); k++)
            if(N(k,k) == 0.)
              N(k,k) += 1.0;
        }
      }
      normal.cholesky(TRUE, 0, eliminationBlocks, TRUE);
      normal.triangularTransSolve(rhs, 0, eliminationBlocks);
      normal.eraseBlocks(0, eliminationBlocks);
      info.observationCount -= eliminationCount;
      for(UInt i=0; i<info.lPl.rows(); i++)
        info.lPl(i) -= quadsum(rhs.slice(0, i, eliminationCount, 1)); // lPl = lPl - n1' N1^(-1) n1
      rhs = rhs.row(eliminationCount, rhs.rows()-eliminationCount);
    }
    else
      logWarningOnce<<"no parameters eliminated"<<Log::endl;

    logStatus<<"write normal equations to <"<<outName<<">"<<Log::endl;
    writeFileNormalEquation(outName, info, normal, rhs);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
