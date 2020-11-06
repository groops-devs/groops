/***********************************************/
/**
* @file fileMatrix.h
*
* @brief Read/write Matrix.
*
* @author Torsten Mayer-Guerr
* @date 2013-02-07
*
*/
/***********************************************/

#ifndef __GROOPS_FILEMATRIX__
#define __GROOPS_FILEMATRIX__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_Matrix
static const char *docstringMatrix = R"(
Stores matrices and vectors. Only one triangle is written for symmetric or triangular matrices.

The header (the matrix definition) is optional.
Therefore a pure text with only numbers in columns are also allowed.
This simplifies the handling of external data.

Instead of a matrix file also an \file{instrument}{instrument} file is allowed.
The first column is the time [MJD], the other columns depends on the instrument type.

\begin{verbatim}
groops matrix version=20200123
LowerSymmetricMatrix( 4 x 4 )
  1.000000000000000000e+00
  0.000000000000000000e+00  1.000000000000000000e+00
  0.000000000000000000e+00  0.000000000000000000e+00  1.000000000000000000e+00
  0.000000000000000000e+00  0.000000000000000000e+00  0.000000000000000000e+00  1.000000000000000000e+00
\end{verbatim}

)";
#endif

/***********************************************/

#include "base/matrix.h"
#include "inputOutput/fileName.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_MATRIX_TYPE = "matrix";

/***** FUNCTIONS *******************************/

/** @brief Write into a Matrix file. */
void writeFileMatrix(const FileName &fileName, const const_MatrixSlice &x);

/** @brief Read from a Matrix file. */
void readFileMatrix(const FileName &fileName, Matrix &x);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
