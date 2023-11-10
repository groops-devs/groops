/***********************************************/
/**
* @file matrixGeneratorRankKUpdate.h
*
* @brief Rank k update (A^T*A).
*
* @author Torsten Mayer-Guerr
* @date 2019-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORRANKKUPDATE__
#define __GROOPS_MATRIXGENERATORRANKKUPDATE__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorRankKUpdate = R"(
\subsection{RankKUpdate}
Symmetric matrix from rank k update: $\M A^T\M A$.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Rank k update (A^T*A).
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorRankKUpdate : public MatrixGeneratorBase
{
  Double             factor;
  MatrixGeneratorPtr matrix;

public:
  MatrixGeneratorRankKUpdate(Config &config);
  void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorRankKUpdate::MatrixGeneratorRankKUpdate(Config &config)
{
  try
  {
    readConfig(config, "matrix", matrix, Config::MUSTSET, "",    "");
    readConfig(config, "factor", factor, Config::DEFAULT, "1.0", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorRankKUpdate::compute(Matrix &N, UInt /*rowsBefore*/, UInt /*columnsBefore*/, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    Matrix A = matrix->compute();
    N = Matrix(A.columns(), Matrix::SYMMETRIC);
    rankKUpdate(factor, A, N);
    fillSymmetric(N);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
