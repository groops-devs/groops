/***********************************************/
/**
* @file normalsReorderAndAccumulate.cpp
*
* @brief Accumulate normal equations and write it to file.
*
* @author Torsten Mayer-Guerr
* @date 2025-10-30
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program accumulates \configFile{inputfileNormalEquation}{normalEquation}s with respect
to the parameter names and writes the total combined system to \configFile{outputfileNormalequation}{normalEquation}.

The combined normal equation is extended to include all parameter names uniquely from all input normals.
The input normals are sorted so that parameters with the same name are accumulated.
This requires that the names in each normal equation are unique.

The output can be written as multiple small block files with \config{outBlockSize},
or as single block with \config{outBlockSize}=0,
or blocked with respect to the first part of the parameter names (object), if \config{outBlockSize} left empty.

See also \program{NormalsBuild} and \program{NormalsAccumulate}.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"
#include "files/fileParameterName.h"

/***** CLASS ***********************************/

/** @brief Accumulate normal equations and write it to file.
* @ingroup programsGroup */
class NormalsReorderAndAccumulate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsReorderAndAccumulate, PARALLEL, "accumulate normal equations and write to file", NormalEquation)

/***********************************************/

void NormalsReorderAndAccumulate::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOut;
    std::vector<FileName> fileNamesInAll;
    UInt     blockSize = NULLINDEX;

    readConfig(config, "outputfileNormalEquation", fileNameOut,    Config::MUSTSET,  "", "");
    readConfig(config, "inputfileNormalEquation",  fileNamesInAll, Config::MUSTSET,  "", "");
    readConfig(config, "outBlockSize",             blockSize,      Config::OPTIONAL, "2048", "block size for distributing the normal equations, 0: one block, empty: blocking by objects");
    if(isCreateSchema(config)) return;

    // ==================================

    logStatus<<"read normal equations info"<<Log::endl;
    std::vector<FileName>      fileNamesIn;
    std::vector<ParameterName> parameterNames;
    std::set<ParameterName>    sortedNames;
    UInt                       rhsCount = 0;

    Single::forEach(fileNamesInAll.size(), [&](UInt i)
    {
      Matrix rhs;
      NormalEquationInfo info;
      try
      {
        readFileNormalEquation(fileNamesInAll.at(i), info, rhs);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<" continue..."<<Log::endl;
        return;
      }

      fileNamesIn.push_back(fileNamesInAll.at(i));
      rhsCount = std::max(rhsCount, info.lPl.size());
      for(const auto &parameterName : info.parameterName)
        if(sortedNames.find(parameterName) == sortedNames.end())
        {
          sortedNames.insert(parameterName);
          parameterNames.push_back(parameterName);
        }
    });

    logInfo<<"  number of unknown parameters (new): "<<parameterNames.size()<<Log::endl;
    logInfo<<"  number of right hand sides:         "<<rhsCount<<Log::endl;
    if(!parameterNames.size())
      return;

    // ==================================

    std::vector<UInt> blockIndex;
    if(blockSize != NULLINDEX)
    {
      blockIndex = MatrixDistributed::computeBlockIndex(parameterNames.size(), blockSize);
    }
    else
    {
      std::stable_sort(parameterNames.begin(), parameterNames.end(), [](auto &p1, auto &p2){return p1.object < p2.object;});
      blockIndex.push_back(0);
      auto blockObject = parameterNames.front().object;
      for(UInt i=0; i<parameterNames.size(); i++)
        if(parameterNames.at(i).object != blockObject)
        {
          blockIndex.push_back(i);
          blockObject = parameterNames.at(i).object;
        }
      blockIndex.push_back(parameterNames.size());
    }
    logInfo<<"  number of blocks:                   "<<blockIndex.size()-1<<Log::endl;

    // ==================================

    logStatus<<"accumulate normals"<<Log::endl;
    MatrixDistributed  normalsSum;
    Matrix             rhsSum(blockIndex.back(), rhsCount);
    Vector             lPl(rhsCount);
    UInt               observationCount = 0;
    normalsSum.initEmpty(blockIndex, comm);
    // write parameter name file
    writeFileParameterName(fileNameOut.replaceFullExtension(".parameterNames.txt"), parameterNames);

    Single::forEach(fileNamesIn.size(), [&](UInt idxFile)
    {
      MatrixDistributed  normals;
      Matrix             rhs;
      NormalEquationInfo info;
      try
      {
        readFileNormalEquation(fileNamesIn.at(idxFile), info, normals, rhs, comm);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<" continue..."<<Log::endl;
        return;
      }

      // create index vector and reorder
      std::vector<UInt>  indexVector;
      auto iter = info.parameterName.begin(); // assume ordered list to accelerate search
      for(const auto &name : parameterNames)
      {
        iter = std::find(iter, info.parameterName.end(), name);
        if(iter == info.parameterName.end()) // not found? -> restart search from begin
          iter = std::find(info.parameterName.begin(), info.parameterName.end(), name);
        indexVector.push_back((iter != info.parameterName.end()) ? std::distance(info.parameterName.begin(), iter) : NULLINDEX);
      }
      normals.reorder(indexVector, blockIndex);
      rhs = reorder(rhs, indexVector);

      // accumulate
      rhsSum.column(0, rhs.columns()) += rhs;
      lPl.row(0, rhs.columns())       += info.lPl;
      observationCount                += info.observationCount;
      for(UInt i=0; i<normals.blockCount(); i++)
        for(UInt k=i; k<normals.blockCount(); k++)
          if(normals.isMyRank(i, k))
          {
            normalsSum.setBlock(i, k);
            normalsSum.N(i, k) += normals.N(i, k);
          }
    });

    // ==================================

    logStatus<<"write normal equations to <"<<fileNameOut<<">"<<Log::endl;
    writeFileNormalEquation(fileNameOut, NormalEquationInfo(parameterNames, lPl, observationCount), normalsSum, rhsSum);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
