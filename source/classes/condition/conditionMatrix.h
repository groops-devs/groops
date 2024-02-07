/***********************************************/
/**
* @file conditionMatrix.h
*
* @brief Evaluate elements of a matrix based on an expression.
*
* @author Sebastian Strasser
* @date 2021-02-03
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONMATRIX__
#define __GROOPS_CONDITIONMATRIX__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionMatrix = R"(
\subsection{Matrix}
Evaluate elements of a \configClass{matrix}{matrixGeneratorType} based on an expression.
If \config{all}=\verb|yes|, all elements of the matrix must evaluate to true
for the condition to be fulfilled, otherwise any element evaluating to true is sufficient.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/condition/condition.h"
#include "classes/matrixGenerator/matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Evaluate elements of a matrix based on an expression.
* @ingroup conditionGroup
* @see Condition */
class ConditionMatrix : public Condition
{
  ExpressionVariablePtr expr;
  MatrixGeneratorPtr matrixGenerator;
  Bool all;

public:
  ConditionMatrix(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionMatrix::ConditionMatrix(Config &config)
{
  try
  {
    readConfig(config, "matrix",     matrixGenerator, Config::MUSTSET, "",  "expression is evaluated for each element of resulting matrix");
    readConfig(config, "expression", expr,            Config::MUSTSET, "",  "(variable: data) evaluated for each element");
    readConfig(config, "all",        all,             Config::DEFAULT, "0", "all (=yes)/any (=no) elements must evaluate to true");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool ConditionMatrix::condition(const VariableList &varList) const
{
  try
  {
    Matrix A = matrixGenerator->compute();

    VariableList varListMatrix = varList;
    for(UInt i=0; i<A.rows(); i++)
      for(UInt j=0; j<A.columns(); j++)
      {
        varListMatrix.setVariable("data", A(i,j));
        Bool result = (expr->evaluate(varListMatrix) != 0.);
        if(!all && result)
          return TRUE;
        else if(all && !result)
          return FALSE;
      }

    return all;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
