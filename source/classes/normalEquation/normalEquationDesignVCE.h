/***********************************************/
/**
* @file normalEquationDesignVCE.h
*
* @brief Accumulate normals from observation equations.
* With individual weights of each arc by variance component estimation.
* @f[ N = A^TPA,\quad n=A^TPl @f]
* @see NormalEquation
* @see Observation
*
* @author Torsten Mayer-Guerr
* @date 2004-12-10
*
*/
/***********************************************/

#ifndef __GROOPS_NORMALEQUATIONDESIGNVCE__
#define __GROOPS_NORMALEQUATIONDESIGNVCE__

// Latex documentation
#ifdef DOCSTRING_NormalEquation
static const char *docstringNormalEquationDesignVCE = R"(
\subsection{DesignVCE}\label{normalEquationType:designVCE}
This class acculumates normal equations from observation equations.
The class \configClass{observation}{observationType} computes
the linearized and decorrelated equation system for each arc $i$:
\begin{equation}
\M l_i  = \M A_i \M x + \M B_i \M y_i + \M e_i.
\end{equation}
The arc depending parameters~$\M y_i$ are eliminated and the system of normal
equations is acculumated according to
\begin{equation}
 \M N =  \sum_{i=1} \frac{1}{\sigma_i^2}\M A_i^T  \M A_i
 \qquad\text{and}\qquad
\M n = \sum_{i=1} \frac{1}{\sigma_i^2} \M A_i^T \M l_i.
\end{equation}
The variance $\sigma_i^2$ of each individual arc is determined by
\begin{equation}
\sigma_i^2 =
\frac{(\M l_i-\M A_i\M x)^T(\M l_i-\M A_i\M x)}
{n_i-\frac{1}{\sigma_i^2}\text{trace}\left(\M A_i^T  \M A_i\M N_{total}^{-1}\right)},
\end{equation}
where $n_i$ is the number of observations. If an apriori solution is not given at the first
iteration step a zero vector is assumed.
)";
#endif

/***********************************************/

#include "classes/observation/observation.h"
#include "classes/normalEquation/normalEquation.h"

/***** CLASS ***********************************/

/** @brief Accumulate normals from observation equations.
* @ingroup normalEquationGroup
* With individual weights of each arc by variance component estimation.
* @f[ N = A^TPA,\quad n=A^TPl @f]
* @see NormalEquation
* @see Observation */
class NormalEquationDesignVCE : public NormalEquationBase
{
  UInt                 iter;
  UInt                 startIndex;
  std::vector<Double>  sigma2;
  ObservationPtr       observation;
  std::vector<UInt>    intervals;

public:
  NormalEquationDesignVCE(Config &config);

  UInt   rightHandSideCount() const override {return observation->rightSideCount();}
  UInt   parameterCount()     const override {return observation->parameterCount() + startIndex;}
  void   parameterNames(std::vector<ParameterName> &names) const override;
  void   init(MatrixDistributed &normals, UInt rhsCount) override;
  Bool   addNormalEquation(UInt rhsNo, const const_MatrixSlice &x, const const_MatrixSlice &Wz,
                           MatrixDistributed &normals, Matrix &n, Vector &lPl, UInt &obsCount) override;
  Vector contribution(MatrixDistributed &Cov) override;
  std::vector<Double> varianceComponentFactors() const override {return sigma2;}
};

/***********************************************/

#endif
