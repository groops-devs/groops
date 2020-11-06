/***********************************************/
/**
* @file matrixGeneratorFile.h
*
* @brief Matrix from file.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORFILE__
#define __GROOPS_MATRIXGENERATORFILE__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorFile = R"(
\subsection{File}
Matrix from \file{file}{matrix}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileMatrix.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Matrix from file.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorFile : public MatrixGeneratorBase
{
  FileName fileName;
  Double   factor;

public:
  MatrixGeneratorFile(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorFile::MatrixGeneratorFile(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    readConfig(config, "inputfileMatrix", fileName, Config::MUSTSET,  "",    "");
    readConfig(config, "factor",          factor,   Config::DEFAULT,  "1.0", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorFile::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    readFileMatrix(fileName, A);
    if(A.getType() == Matrix::SYMMETRIC)
      fillSymmetric(A);
    A *= factor;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
