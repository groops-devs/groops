/***********************************************/
/**
* @file matrixGeneratorNormalsFile.h
*
* @brief Symmetric matrix from a normal equation file.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORNORMALSFILE__
#define __GROOPS_MATRIXGENERATORNORMALSFILE__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorNormalsFile = R"(
\subsection{Normals file}
Matrix from a \file{normal equation file}{normalEquation}. The symmetric normal matrix,
the right hand side vector or the lPl vector can be selected.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileNormalEquation.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Symmetric matrix from a normal equation file.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorNormalsFile : public MatrixGeneratorBase
{
  enum Type {MATRIX, RHS, LPL};
  Type type;
  FileName fileName;
  Double   factor;

public:
  MatrixGeneratorNormalsFile(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorNormalsFile::MatrixGeneratorNormalsFile(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    type = MATRIX;
    std::string choice;

    renameDeprecatedConfig(config, "inputfileNormalequation", "inputfileNormalEquation", date2time(2020, 6, 3));

    readConfig(config, "inputfileNormalEquation", fileName, Config::MUSTSET,  "",    "");
    if(readConfigChoice(config, "type", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "normalMatrix",  choice, "")) type = MATRIX;
      if(readConfigChoiceElement(config, "rightHandSide", choice, "")) type = RHS;
      if(readConfigChoiceElement(config, "lPl",           choice, "")) type = LPL;
      endChoice(config);
    }
    readConfig(config, "factor", factor, Config::DEFAULT,  "1.0", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorNormalsFile::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    Matrix n;
    NormalEquationInfo info;

    if(type == MATRIX)
    {
      readFileNormalEquation(fileName, info, A, n);
      fillSymmetric(A);
    }
    else if(type == RHS)
    {
      readFileNormalEquation(fileName, info, n);
      A = n;
    }
    else if(type == LPL)
    {
      readFileNormalEquation(fileName, info, n);
      A = info.lPl;
    }

    A *= factor;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
