/***********************************************/
/**
* @file fileNormalEquation.h
*
* @brief Read/write a system of normal equations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-05-27
*
*/
/***********************************************/

#ifndef __GROOPS_FILENORMALEQUATION__
#define __GROOPS_FILENORMALEQUATION__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_NormalEquation
static const char *docstringNormalEquation = R"(
Stores a \reference{system of normal equations}{normalEquationType}
\begin{equation}
  \M N \hat{\M x} = \M n.
\end{equation}.
This file format consists of multiple files.
The file name \verb|normals.dat.gz| corresponds to the following files:
\begin{itemize}
\item \verb|normals.dat.gz| or \verb|normals.00.00.dat.gz| ... \verb|normals.0n.0n.dat.gz|:
      the normal matrix $\M N$ as \file{matrix}{matrix},
\item \verb|normals.rightHandSide.dat.gz|:
      the right hand side(s) $\M n$ as \file{matrix}{matrix},
\item \verb|normals.parameterNames.txt|: \file{parameter names}{parameterName},
\item \verb|normals.info.xml|:
\     u.a. containing the number of observations and the quadratic sum of (reduced) observations $\M l^T\M P\M l$.
\end{itemize}
A large normal matrix may be splitted into blocks and stored in multiple files.
The block row/column number is indicated in the file name.
Only the upper blocks of the sysmmetric matrix are considered.
Matrix in blocks can be distributed on muliple nodes in parallel mode to efficiently use distributed memory.
)";
#endif

/***********************************************/

#include "base/parameterName.h"
#include "inputOutput/fileName.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_NORMALEQUATION_TYPE = "normalEquation";

/***** CLASSES **********************************/

class NormalEquationInfo
{
public:
  std::vector<ParameterName> parameterName;
  Vector lPl;
  UInt observationCount;
  std::vector<UInt> blockIndex;
  Matrix usedBlocks;

  NormalEquationInfo(UInt rhsCount = 1) : lPl(rhsCount), observationCount(0) {}
  NormalEquationInfo(const std::vector<ParameterName> &parameterName_, const Vector &lPl_ = Vector(1), UInt observationCount_ = 0)
    : parameterName(parameterName_), lPl(lPl_), observationCount(observationCount_) {}
};

/***** FUNCTIONS ********************************/

/** @brief Write a system of normal equations. */
void writeFileNormalEquation(const FileName &name, NormalEquationInfo info, const Matrix &N, const Matrix &n);

/** @brief Write a system of normal equations. */
void writeFileNormalEquation(const FileName &name, NormalEquationInfo info, const std::vector<std::vector<Matrix>> &N, const Matrix &n);

class MatrixDistributed;
/** @brief Write a system of normal equations.
* Must be called on every process.
* Only at master needed: @a parameterName, @a n, @a lPl, @a observationCount. */
void writeFileNormalEquation(const FileName &name, NormalEquationInfo info, const MatrixDistributed &N, const Matrix &n);

/** @brief Write a system of normal equations.
* Only the information file, parameter name file and the right hand sides are written. */
void writeFileNormalEquation(const FileName &name, NormalEquationInfo info, const Matrix &n);

/***********************************************/

/** @brief Read a system of normal equations. */
void readFileNormalEquation(const FileName &name, NormalEquationInfo &info, Matrix &N, Matrix &n);

/** @brief Write a system of normal equations.
* Must be called on every process. */
void readFileNormalEquation(const FileName &name, NormalEquationInfo &info, MatrixDistributed &normal, Matrix &n);

/** @brief Read a system of normal equations.
* Only the information file, parameter name file and the right hand sides are read. */
void readFileNormalEquation(const FileName &name, NormalEquationInfo &info, Matrix &n);

/***********************************************/

/// @}

#endif /* __GROOPS__ */
