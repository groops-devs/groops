/***********************************************/
/**
* @file conditionMatrixEmpty.h
*
* @brief Evaluate if matrix (or instrument) file is empty/has zero size.
*
* @author Sebastian Strasser
* @date 2021-09-15
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONMATRIXEMPTY__
#define __GROOPS_CONDITIONMATRIXEMPTY__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionMatrixEmpty = R"(
\subsection{MatrixEmpty}
Evaluate if \file{matrix}{matrix} (or \file{instrument}{instrument}) file is empty/has zero size.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileMatrix.h"
#include "classes/condition/condition.h"

/***** CLASS ***********************************/

/** @brief Evaluate if matrix (or instrument) file is empty/has zero size.
* @ingroup ConditionGroup
* @see Condition */
class ConditionMatrixEmpty : public Condition
{
  FileName fileName;

public:
  ConditionMatrixEmpty(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionMatrixEmpty::ConditionMatrixEmpty(Config &config)
{
  try
  {
    readConfig(config, "inputfileMatrix", fileName, Config::MUSTSET, "",  "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool ConditionMatrixEmpty::condition(const VariableList &varList) const
{
  try
  {
    Matrix A;
    readFileMatrix(fileName(varList), A);
    return !A.size();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
