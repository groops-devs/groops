/***********************************************/
/**
* @file parameterSelection2IndexVector.cpp
*
* @brief Generate index vector from parameter selection.
*
* @author Sebastian Strasser
* @date 2018-05-18
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Generate index vector from parameter selection in \file{matrix format}{matrix}.
This vector can be used in \program{MatrixCalculate}
with \configClass{matrix:reorder}{matrixGeneratorType:reorder}
to reorder arbitary vectors and matrices similar to \program{NormalsReorder}.

The \configClass{parameterSelection}{parameterSelectorType} allows reordering and dimension changes,
either by cutting parameters or by inserting additional parameters.
\configFile{outputfileIndexVector}{matrix} contains indices of parameters in
\configFile{inputfileParameterNames}{parameterName} or -1 for newly added parameters.
\configFile{outputfileParameterNames}{parameterName} contains the selected parameter names.
)";

/***********************************************/

#include "programs/program.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"

/***** CLASS ***********************************/

/** @brief Generate index vector from parameter selection.
* @ingroup programsGroup */
class ParameterSelection2IndexVector
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(ParameterSelection2IndexVector, SINGLEPROCESS, "Generate index vector from parameter selection.", NormalEquation)
GROOPS_RENAMED_PROGRAM(NormalsParameterSelection2IndexVector, ParameterSelection2IndexVector, date2time(2018, 7, 4))

/***********************************************/

void ParameterSelection2IndexVector::run(Config &config)
{
  try
  {
    FileName fileNameIndexVector, fileNameParameterNamesOut, fileNameParameterNamesIn;
    ParameterSelectorPtr selector;

    readConfig(config, "outputfileIndexVector",    fileNameIndexVector,       Config::OPTIONAL,  "", "indices of source parameters in target normal equations");
    readConfig(config, "outputfileParameterNames", fileNameParameterNamesOut, Config::OPTIONAL,  "", "output parameter names file");
    readConfig(config, "inputfileParameterNames",  fileNameParameterNamesIn,  Config::MUSTSET,   "", "parameter names file of source normal equations");
    readConfig(config, "parameterSelection",       selector,                  Config::MUSTSET,   "", "parameter order/selection of target normal equations");
    if(isCreateSchema(config)) return;

    std::vector<ParameterName> parameterNames;
    readFileParameterName(fileNameParameterNamesIn, parameterNames);
    const std::vector<UInt> indexVector = selector->indexVector(parameterNames);

    if(!fileNameIndexVector.empty())
    {
      Vector vector(indexVector.size());
      for(UInt i = 0; i < indexVector.size(); i++)
        vector(i) = (indexVector.at(i) != NULLINDEX ? static_cast<Double>(indexVector.at(i)) : -1);

      logStatus<<"writing index vector <"<<fileNameIndexVector<<">"<<Log::endl;
      writeFileMatrix(fileNameIndexVector, vector);
    }

    if(!fileNameParameterNamesOut.empty())
    {
      std::vector<ParameterName> parameterNamesOut(indexVector.size());
      for(UInt i = 0; i < indexVector.size(); i++)
        if(indexVector.at(i) != NULLINDEX)
          parameterNamesOut.at(i) = parameterNames.at(indexVector.at(i));

      logStatus<<"writing parameter names file <"<<fileNameParameterNamesOut<<">"<<Log::endl;
      writeFileParameterName(fileNameParameterNamesOut, parameterNamesOut);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
