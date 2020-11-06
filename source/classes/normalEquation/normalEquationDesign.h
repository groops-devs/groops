/***********************************************/
/**
* @file normalEquationDesign.h
*
* @brief Accumulate normals from observation equations.
* @f[ N = A^TPA,\quad n=A^TPl @f]
* @see NormalEquation
* @see Observation
*
* @author Torsten Mayer-Guerr
* @date 2004-12-10
*
*/
/***********************************************/

#ifndef __GROOPS_NORMALEQUATIONDESIGN__
#define __GROOPS_NORMALEQUATIONDESIGN__

// Latex documentation
#ifdef DOCSTRING_NormalEquation
static const char *docstringNormalEquationDesign = R"(
\subsection{Design}\label{normalEquationType:design}
This class acculumates normal equations from observation equations.
The class \configClass{observation}{observationType} computes
the linearized and decorrelated equation system for each arc $i$:
\begin{equation}
\M l_i  = \M A_i \M x + \M B_i \M y_i + \M e_i.
\end{equation}
The arc depending parameters~$\M y_i$ are eliminated and the system of normal
equations is acculumated according to
\begin{equation}
 \M N = \sum_{i=1}^m \M A_i^T  \M A_i
 \qquad\text{and}\qquad
\M n = \sum_{i=1}^m \M A_i^T \M l_i.
\end{equation}
)";
#endif

/***********************************************/

#include "classes/observation/observation.h"
#include "classes/normalEquation/normalEquation.h"

/***** CLASS ***********************************/

/** @brief Accumulate normals from observation equations.
* @ingroup normalEquationGroup
* @f[ N = A^TPA,\quad n=A^TPl @f]
* @see NormalEquation
* @see Observation */
class NormalEquationDesign : public NormalEquationBase
{
  UInt              startIndex;
  std::vector<UInt> intervals;
  ObservationPtr    observation;
  Double            sigma2, sigma2New;

public:
  NormalEquationDesign(Config &config);

  UInt   rightHandSideCount() const override {return observation->rightSideCount();}
  UInt   parameterCount()     const override {return observation->parameterCount() + startIndex;}
  void   parameterNames(std::vector<ParameterName> &names) const override;
  void   init(MatrixDistributed &normals, UInt rhsCount) override;
  Bool   addNormalEquation(UInt rhsNo, const const_MatrixSlice &x, const const_MatrixSlice &Wz,
                           MatrixDistributed &normals, Matrix &n, Vector &lPl, UInt &obsCount) override;
  Vector contribution(MatrixDistributed &Cov) override;
  std::vector<Double> varianceComponentFactors() const override {return std::vector<Double>({sigma2});}
};

/***********************************************/

#endif
