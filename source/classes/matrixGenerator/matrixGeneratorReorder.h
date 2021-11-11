/***********************************************/
/**
* @file matrixGeneratorReorder.h
*
* @brief Reorder matrix by index vectors.
*
* @author Andreas Kvas
* @date 2017-11-30
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORREORDER__
#define __GROOPS_MATRIXGENERATORREORDER__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorReorder = R"(
\subsection{Reorder}\label{matrixGeneratorType:reorder}
Reorder rows or columns of a matrix by an index vectors.
The index vector can be created with \program{ParameterSelection2IndexVector}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileMatrix.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Reorder matrix by index vectors.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorReorder : public MatrixGeneratorBase
{
  FileName fileNameRow, fileNameColumn;
  std::vector<UInt>  rowIndex;
  std::vector<UInt>  columnIndex;
  MatrixGeneratorPtr matrix;

public:
  MatrixGeneratorReorder(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorReorder::MatrixGeneratorReorder(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    renameDeprecatedConfig(config, "rowIndex",    "inputfileIndexVectorRow",    date2time(2018, 6, 6));
    renameDeprecatedConfig(config, "columnIndex", "inputfileIndexVectorColumn", date2time(2018, 6, 6));

    readConfig(config, "matrix",                     matrix,         Config::MUSTSET,  "", "");
    readConfig(config, "inputfileIndexVectorRow",    fileNameRow,    Config::OPTIONAL, "", "index in input matrix or -1 for new parameter.");
    readConfig(config, "inputfileIndexVectorColumn", fileNameColumn, Config::OPTIONAL, "", "index in input matrix or -1 for new parameter.");
    if(isCreateSchema(config)) return;

    // generate row index
    // ------------------
    if(!fileNameRow.empty())
    {
      Vector tmp;
      readFileMatrix(fileNameRow, tmp);

      rowIndex.resize(tmp.rows());
      for(UInt k = 0; k<rowIndex.size(); k++)
        rowIndex[k] = (tmp(k)<0.0) ? NULLINDEX : static_cast<UInt>(tmp(k));
    }

    // generate row index
    // ------------------
    if(!fileNameColumn.empty())
    {
      Vector tmp;
      readFileMatrix(fileNameColumn, tmp);

      columnIndex.resize(tmp.rows());
      for(UInt k = 0; k<columnIndex.size(); k++)
        columnIndex[k] = (tmp(k)<0.0) ? NULLINDEX : static_cast<UInt>(tmp(k));
    }

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorReorder::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();
    if(!fileNameRow.empty())
      A = reorder(A, rowIndex);
    if(!fileNameColumn.empty())
      A = reorder(A.trans(), columnIndex).trans();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
