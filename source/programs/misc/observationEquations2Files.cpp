/***********************************************/
/**
* @file observationEquations2Files.cpp
*
* @brief Write the design matrix and observation vector.
*
* @author Torsten Mayer-Guerr
* @date 2022-12-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the linearized and decorrelated equation system for each arc $i$:
\begin{equation}
\M l_i  = \M A_i \M x + \M B_i \M y_i + \M e_i
\end{equation}
using class \configClass{observation}{observationType} and writes $\M A_i$, $\M B_i$ and $\M l_i$ as \file{matrix}{matrix} files.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"
#include "classes/observation/observation.h"

/***** CLASS ***********************************/

/** @brief Write the design matrix and observation vector.
* @ingroup programsGroup */
class ObservationEquations2Files
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(ObservationEquations2Files, PARALLEL, "write the design matrix and observation vector", Misc)

/***********************************************/

void ObservationEquations2Files::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName       fileNameL, fileNameA, fileNameB;
    FileName       fileNamePara;
    std::string    nameArc;
    ObservationPtr observation;

    readConfig(config, "outputfileObservationVector", fileNameL,    Config::OPTIONAL, "l_{arc:%03i}.dat",   "one file for each arc");
    readConfig(config, "outputfileDesignMatrix",      fileNameA,    Config::OPTIONAL, "A_{arc:%03i}.dat",   "one file for each arc, without arc related parameters");
    readConfig(config, "outputfileDesignMatrixArc",   fileNameB,    Config::OPTIONAL, "",                   "one file for each arc, arc related parameters");
    readConfig(config, "variableArc",                 nameArc,      Config::OPTIONAL, "arc",                "variable with arc number");
    readConfig(config, "outputfileParameterNames",    fileNamePara, Config::OPTIONAL, "parameterNames.txt", "without arc related parameters");
    readConfig(config, "observation",                 observation,  Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    if(Parallel::isMaster(comm) && !fileNamePara.empty())
    {
      logStatus<<"write parameter names to <"<<fileNamePara<<">"<<Log::endl;
      std::vector<ParameterName> names;
      observation->parameterName(names);
      writeFileParameterName(fileNamePara, names);
    }

    if(fileNameL.empty() && fileNameA.empty() && fileNameB.empty())
      return;

    logStatus<<"setup observation equations"<<Log::endl;
    Parallel::forEach(observation->arcCount(), [&] (UInt arcNo)
    {
      Matrix l, A, B;
      observation->observation(arcNo, l, A, B);

      VariableList fileNameVariableList;
      fileNameVariableList.setVariable(nameArc, arcNo);
      if(!fileNameL.empty()) writeFileMatrix(fileNameL(fileNameVariableList), l);
      if(!fileNameA.empty()) writeFileMatrix(fileNameA(fileNameVariableList), A);
      if(!fileNameB.empty()) writeFileMatrix(fileNameB(fileNameVariableList), B);
    }, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
