/***********************************************/
/**
* @file matrixGeneratorAppend.h
*
* @brief Append matrix to right or bottom.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXCGENERATORAPPEND__
#define __GROOPS_MATRIXCGENERATORAPPEND__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorAppend = R"(
\subsection{Append}
Append matrix to the right (first row) or bottom (first column).
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Append of a matrix.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorAppend : public MatrixGeneratorBase
{
  Bool               right, bottom;
  MatrixGeneratorPtr matrix;

public:
  MatrixGeneratorAppend(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorAppend::MatrixGeneratorAppend(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    std::string choice;
    readConfig(config, "matrix", matrix, Config::MUSTSET, "", "");
    if(readConfigChoice(config, "side", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "right",    choice, "")) {right = TRUE;  bottom = FALSE;}
      if(readConfigChoiceElement(config, "bottom",   choice, "")) {right = FALSE; bottom = TRUE;}
      if(readConfigChoiceElement(config, "diagonal", choice, "")) {right = TRUE;  bottom = TRUE;}
      endChoice(config);
    }
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorAppend::compute(Matrix &A, UInt &startRow, UInt &startCol)
{
  try
  {
    A = matrix->compute();
    if(right)
      startCol = static_cast<UInt>(varList["columnsBefore"]->evaluate(varList));
    if(bottom)
      startRow = static_cast<UInt>(varList["rowsBefore"]->evaluate(varList));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
