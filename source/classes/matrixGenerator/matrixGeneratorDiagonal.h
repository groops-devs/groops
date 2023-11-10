/***********************************************/
/**
* @file matrixGeneratorDiagonal.h
*
* @brief Exctract the diagonal of a matrix.
*
* @author Andreas Kvas
* @date 2018-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORDIAGONAL__
#define __GROOPS_MATRIXGENERATORDIAGONAL__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorDiagonal = R"(
\subsection{Diagonal}
Extract the diagonal or subdiagnoal ($n\times 1$ vector) of a matrix.
The zero \config{diagonal} means the main diagonal, a positive value the superdiagonal,
and a negative the subdiagonal.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Exctract the diagonal of a matrix.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorDiagonal : public MatrixGeneratorBase
{
  MatrixGeneratorPtr matrix;
  Int index;

public:
  MatrixGeneratorDiagonal(Config &config);
  void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorDiagonal::MatrixGeneratorDiagonal(Config &config)
{
  try
  {
    readConfig(config, "matrix",   matrix, Config::MUSTSET,  "",  "");
    readConfig(config, "diagonal", index,  Config::DEFAULT,  "0", "zero: main diagonal, positive: superdiagonal, negative: subdiagonal");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorDiagonal::compute(Matrix &A, UInt /*rowsBefore*/, UInt /*columnsBefore*/, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();

    Vector d(std::min(A.columns(), A.rows()) - std::abs(index));
    for(UInt k = 0; k<d.rows(); k++)
      d(k) = A(k-std::min(0, index), k+std::max(0, index));

    A = d;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
