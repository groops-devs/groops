/***********************************************/
/**
* @file normalsReorder.cpp
*
* @brief Reorder normal equations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2011-02-15
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Reorder \configFile{inputfileNormalEquation}{normalEquation} by selecting parameters in a specific order.
The \configClass{parameterSelection}{parameterSelectorType} also allows to change dimension of the normal equations,
either by cutting parameters or by inserting zero rows/columns for additional parameters.
Without \configClass{parameterSelection}{parameterSelectorType} the order of parameters remains the same.
Additionally the block sizes of the files can be adjusted. If \config{outBlockSize} is set to zero,
the normal matrix is written to a single block file, which is needed by some programs.

To eliminate parameters without changing the result of the other parameters use \program{NormalsEliminate}.
)";

/***********************************************/

#include "programs/program.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"

/***** CLASS ***********************************/

/** @brief Reorder normal equations.
* @ingroup programsGroup */
class NormalsReorder
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(NormalsReorder, PARALLEL, "Reorder normal equations", NormalEquation)

/***********************************************/

void NormalsReorder::run(Config &config)
{
  try
  {
    FileName outName, inName;
    ParameterSelectorPtr parameterSelector;
    UInt     blockSize;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "inputfileNormalequation",  "inputfileNormalEquation",  date2time(2020, 6, 3));

    readConfig(config, "outputfileNormalEquation", outName,           Config::MUSTSET,  "", "");
    readConfig(config, "inputfileNormalEquation",  inName,            Config::MUSTSET,  "", "");
    readConfig(config, "parameterSelection",       parameterSelector, Config::OPTIONAL, "",     "parameter order/selection of output normal equations");
    readConfig(config, "outBlockSize",             blockSize,         Config::DEFAULT,   "2048", "block size for distributing the normal equations, 0: one block");
    if(isCreateSchema(config)) return;

    // ==================================

    logStatus<<"init normal equations"<<Log::endl;
    MatrixDistributed normal;
    Matrix rhs;
    NormalEquationInfo info;
    readFileNormalEquation(inName, info, normal, rhs);

    std::vector<UInt> indexVector(normal.parameterCount());
    if(parameterSelector)
      indexVector = parameterSelector->indexVector(info.parameterName);
    else
      std::iota(indexVector.begin(), indexVector.end(), 0);

    logInfo<<"  number of unknown parameters (old): "<<normal.parameterCount()<<Log::endl;
    logInfo<<"  number of unknown parameters (new): "<<indexVector.size()<<Log::endl;
    logInfo<<"  number of right hand sides:         "<<info.lPl.rows()<<Log::endl;

    logStatus<<"reorder normal matrix"<<Log::endl;
    normal.reorder(indexVector, MatrixDistributed::computeBlockIndex(indexVector.size(), blockSize));
    rhs = reorder(rhs, indexVector);
    std::vector<ParameterName> parameterNames(indexVector.size());
    for(UInt i = 0; i < indexVector.size(); i++)
      parameterNames.at(i) = (indexVector.at(i) != NULLINDEX ? info.parameterName.at(indexVector.at(i)) : ParameterName());

    logStatus<<"write normal equations to <"<<outName<<">"<<Log::endl;
    writeFileNormalEquation(outName, NormalEquationInfo(parameterNames, info.lPl, info.observationCount), normal, rhs);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
