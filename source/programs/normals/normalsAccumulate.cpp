/***********************************************/
/**
* @file normalsAccumulate.cpp
*
* @brief Accumulate normal equations and write it to file.
*
* @author Andreas Kvas
* @date 2013-10-26
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program accumulates normal equations and writes the total combined system to
\configFile{outputfileNormalequation}{normalEquation}.
The \configFile{inputfileNormalEquation}{normalEquation}s must have all the same size and the same block structure.
This program is the simplified and fast version of the more general program \program{NormalsBuild}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"

/***** CLASS ***********************************/

/** @brief Accumulate normal equations and write it to file.
* @ingroup programsGroup */
class NormalsAccumulate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsAccumulate, SINGLEPROCESS, "accumulate normal equations and write to file", NormalEquation)

/***********************************************/

void NormalsAccumulate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut;
    std::vector<FileName> fileNameInAll;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "inputfileNormalequation",  "inputfileNormalEquation",  date2time(2020, 6, 3));

    readConfig(config, "outputfileNormalEquation", fileNameOut,   Config::MUSTSET,  "", "");
    readConfig(config, "inputfileNormalEquation",  fileNameInAll, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // ==================================

    logStatus<<"read normal equations info"<<Log::endl;
    Matrix nOut;
    NormalEquationInfo infoOut;
    std::vector<FileName> fileNameIn;
    std::vector<Matrix>   usedBlocksIn;
    for(UInt i=0; i<fileNameInAll.size(); i++)
    {
      Matrix nIn;
      NormalEquationInfo infoIn;
      try
      {
        readFileNormalEquation(fileNameInAll.at(i), infoIn, nIn);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<" continue..."<<Log::endl;
        continue;
      }

      fileNameIn.push_back(fileNameInAll.at(i));
      usedBlocksIn.push_back(infoIn.usedBlocks);
      if(!infoOut.blockIndex.size())
      {
        infoOut = infoIn;
        nOut    = nIn;
        continue;
      }

      if(infoOut.blockIndex.size() != infoIn.blockIndex.size())
        throw(Exception("normals must have the same size and block structure"));
      for(UInt z=0; z<infoOut.blockIndex.size(); z++)
        if(infoOut.blockIndex.at(z) != infoIn.blockIndex.at(z))
          throw(Exception("normals must have the same size and block structure"));
      for(UInt z=0; z<infoOut.blockIndex.size()-1; z++)
        for(UInt s=z; s<infoOut.blockIndex.size()-1; s++)
          if(infoIn.usedBlocks(z,s))
            infoOut.usedBlocks(z,s) = 1;
      for(UInt i=0; i<infoOut.parameterName.size(); i++)
        if(!infoOut.parameterName.at(i).combine(infoIn.parameterName.at(i)))
          logWarning << "Parameter names do not match at index " << i << ": '" << infoOut.parameterName.at(i).str() << "' != '" << infoIn.parameterName.at(i).str() << "'" << Log::endl;

      nOut += nIn;
      infoOut.lPl              += infoIn.lPl;
      infoOut.observationCount += infoIn.observationCount;
    }

    // ==================================

    logStatus<<"read and write block normals"<<Log::endl;
    Log::Timer timer(infoOut.blockIndex.size()*(infoOut.blockIndex.size()-1)/2);
    UInt idx = 0;
    for(UInt z=0; z<infoOut.blockIndex.size()-1; z++)
      for(UInt s=z; s<infoOut.blockIndex.size()-1; s++)
      {
        timer.loopStep(idx++);

        if(infoOut.usedBlocks(z,s) == 0)
          continue;

        std::string ext;
        if(infoOut.blockIndex.size()-1 > 1)
          ext = "."+z%"%02i-"s+s%"%02i"s;

        Matrix N;
        for(UInt i=0; i<fileNameIn.size(); i++)
          if(usedBlocksIn.at(i)(z,s))
          {
            logStatus<<"read matrix <"<<fileNameIn.at(i).appendBaseName(ext)<<">"<<Log::endl;
            Matrix N2;
            readFileMatrix(fileNameIn.at(i).appendBaseName(ext), N2);
            if(N.size())
              N += N2;
            else
              N = N2;
          }

        logStatus<<"write matrix <"<<fileNameOut.appendBaseName(ext)<<">"<<Log::endl;
        writeFileMatrix(fileNameOut.appendBaseName(ext), N);
      }
    timer.loopEnd();

    // ==================================

    logStatus<<"write normal equations to <"<<fileNameOut<<">"<<Log::endl;
    writeFileNormalEquation(fileNameOut, infoOut, nOut);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
