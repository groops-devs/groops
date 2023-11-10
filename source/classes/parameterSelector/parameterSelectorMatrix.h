/***********************************************/
/**
* @file parameterSelectorMatrix.h
*
* @brief Parameter index vector from matrix.
*
* @author Sebastian Strasser
* @date 2018-05-08
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERSELECTORMATRIX__
#define __GROOPS_PARAMETERSELECTORMATRIX__

// Latex documentation
#ifdef DOCSTRING_ParameterSelector
static const char *docstringParameterSelectorMatrix = R"(
\subsection{Matrix}
Parameter index vector from matrix.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Parameter index vector from matrix.
* @ingroup parameterSelectorGroup
* @see ParameterSelector */
class ParameterSelectorMatrix : public ParameterSelectorBase
{
  FileName fileName;
  ExpressionVariablePtr exprColumn, exprStartRow, exprCountRows;

public:
  ParameterSelectorMatrix(Config &config);
  std::vector<UInt> indexVector(const std::vector<ParameterName> &parameterNames);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParameterSelectorMatrix::ParameterSelectorMatrix(Config &config)
{
  try
  {
    readConfig(config, "inputfileMatrix", fileName,      Config::MUSTSET,  "",  "index in old parameter list or -1 for new parameter");
    readConfig(config, "column",          exprColumn,    Config::DEFAULT,  "0", "use this column (counting from 0, variables: columns)");
    readConfig(config, "startRow",        exprStartRow,  Config::DEFAULT,  "0", "start at this row (counting from 0, variables: rows)");
    readConfig(config, "countRows",       exprCountRows, Config::OPTIONAL, "",  "use these many rows (default: use all, variables: rows)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<UInt> ParameterSelectorMatrix::indexVector(const std::vector<ParameterName> &parameterNames)
{
  try
  {
    Matrix matrix;
    readFileMatrix(fileName, matrix);

    VariableList varList;
    varList.setVariable("rows",    static_cast<Double>(matrix.rows()));
    varList.setVariable("columns", static_cast<Double>(matrix.columns()));
    const UInt   column    = static_cast<UInt>(exprColumn->evaluate(varList));
    const UInt   startRow  = static_cast<UInt>(exprStartRow->evaluate(varList));
    const UInt   countRows = (exprCountRows ? static_cast<UInt>(exprCountRows->evaluate(varList)) : matrix.rows()-startRow);
    const Vector index     = matrix.slice(startRow, column, countRows, 1);

    const UInt parameterCount = parameterNames.size();
    if(max(index) > static_cast<Double>(MAX_UINT))
      throw(Exception("unsigned integer overflow: "+max(index)%"%e > "s+MAX_UINT%"%.0f"s));
    if(max(index) >= parameterCount)
      throw(Exception("matrix contains index that exceeds parameter count: "+max(index)%"%i >= "s+parameterCount%"%i"s));

    std::vector<UInt> vector(index.size());
    for(UInt i=0; i<index.size(); i++)
      vector.at(i) = ((index(i) < 0) ? NULLINDEX : static_cast<UInt>(index(i)));

    return vector;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
