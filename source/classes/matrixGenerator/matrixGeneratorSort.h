/***********************************************/
/**
* @file matrixGeneratorSort.h
*
* @brief Sort matrix by column.
*
* @author Sebastian Strasser
* @date 2018-02-20
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORSORT__
#define __GROOPS_MATRIXGENERATORSORT__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorSort = R"(
\subsection{Sort}
Sort matrix by \config{column} in ascending order by default or in \config{descending} order.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Sort matrix by column.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorSort : public MatrixGeneratorBase
{
  MatrixGeneratorPtr matrix;
  std::vector<UInt>  column;
  Bool descending;

public:
  MatrixGeneratorSort(Config &config);
  void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorSort::MatrixGeneratorSort(Config &config)
{
  try
  {
    readConfig(config, "matrix",     matrix,     Config::MUSTSET, "",  "");
    readConfig(config, "column",     column,     Config::MUSTSET, "0", "sort by column, top = highest priority");
    readConfig(config, "descending", descending, Config::DEFAULT, "0", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorSort::compute(Matrix &A, UInt /*rowsBefore*/, UInt /*columnsBefore*/, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();

    for(auto c : column)
      if(c >= A.columns())
        throw(Exception("column index error: "+c%"&i"s+" >= "+A.columns()%"&i"s));

    std::vector<UInt> rowIndex(A.rows());
    std::iota(rowIndex.begin(), rowIndex.end(), 0);
    for(UInt idCol = column.size(); idCol-->0;)
      std::stable_sort(rowIndex.begin(), rowIndex.end(), [&](UInt i, UInt j){return A(i, column.at(idCol)) < A(j, column.at(idCol));});
    if(descending)
      std::reverse(rowIndex.begin(), rowIndex.end());

    A = reorder(matrix->compute(), rowIndex);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
