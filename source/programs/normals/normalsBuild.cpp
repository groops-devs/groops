/***********************************************/
/**
* @file normalsBuild.cpp
*
* @brief Accumulate normal equations and write it to file.
*
* @author Torsten Mayer-Guerr
* @date 2003-03-26
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program accumulates \configClass{normalEquation}{normalEquationType}s and
writes the total combined system to \configFile{outputfileNormalequation}{normalEquation}.
For a detailed description of the used algorithm see \configClass{normalEquation}{normalEquationType}.
Large normal equation systems can be divided into blocks with \config{normalsBlockSize}.

A simplifed and fast version of this program is \program{NormalsAccumulate}.
For input normals with different parameters see \program{NormalsReorderAndAccumulate}.
To solve the system of normal equations use \program{NormalsSolverVCE}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/normalEquation/normalEquation.h"

/***** CLASS ***********************************/

/** @brief Accumulate normal equations and write it to file.
* @ingroup programsGroup */
class NormalsBuild
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsBuild, PARALLEL, "accumulate normal equations and write it to file", NormalEquation)

/***********************************************/

void NormalsBuild::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName          fileNameNormals;
    NormalEquationPtr normals;
    UInt              blockSize;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "normalequation",           "normalEquation",           date2time(2020, 6, 3));

    readConfig(config, "outputfileNormalEquation", fileNameNormals, Config::MUSTSET,  "", "");
    readConfig(config, "normalEquation",           normals,         Config::MUSTSET,  "", "");
    readConfig(config, "normalsBlockSize",         blockSize,       Config::DEFAULT,  "2048", "block size for distributing the normal equations, 0: one block");
    if(isCreateSchema(config)) return;

    logStatus<<"init normal equations"<<Log::endl;
    normals->init(blockSize, comm);
    logInfo<<"  number of unknown parameters: "<<normals->parameterCount()<<Log::endl;
    logInfo<<"  number of right hand sides:   "<<normals->rightHandSideCount()<<Log::endl;

    logStatus<<"accumulate normal equations"<<Log::endl;
    normals->build();

    logStatus<<"write normal equations to <"<<fileNameNormals<<">"<<Log::endl;
    normals->write(fileNameNormals);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
