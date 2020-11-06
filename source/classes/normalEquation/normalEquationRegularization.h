/***********************************************/
/**
* @file normalEquationRegularization.h
*
* @brief Regularization with a diagonal matrix.
* @see NormalEquation
*
* @author Torsten Mayer-Guerr
* @date 2004-12-10
*
*/
/***********************************************/

#ifndef __GROOPS_NORMALEQUATIONREGULARIZATION__
#define __GROOPS_NORMALEQUATIONREGULARIZATION__

// Latex documentation
#ifdef DOCSTRING_NormalEquation
static const char *docstringNormalEquationRegularization = R"(
\subsection{Regularization}\label{normalEquationType:regularization}
Set up a system of normal equations
\begin{equation}
\M N = \M R
\qquad\text{and}\qquad
\M n = \M R \M b,
\end{equation}
where $\M R$ is a diagonal matrix whose elements are given as a vector by
\configFile{inputfileDiagonalMatrix}{matrix} and $\M b$ is the right hand side towards which will
be regularized. It can be given by \configFile{inputfileBiasVector}{matrix}.
The diagonal matrix can be generated with \program{NormalsRegularizationBorders},
\program{NormalsRegularizationSphericalHarmonics}, or \program{MatrixCalculate}.
If $\M R$ is not given a unit matrix is assumed.
The right hand side $\M b$ may be generated with \program{Gravityfield2SphericalHarmonicsVector}.
If $\M b$ is not given a zero vector is assumed.
)";
#endif

/***********************************************/

#include "classes/normalEquation/normalEquation.h"

/***** CLASS ***********************************/

/** @brief Regularization with a diagonal matrix.
* @ingroup normalEquationGroup
* @see NormalEquation */
class NormalEquationRegularization : public NormalEquationBase
{
  Vector  K;
  Matrix  bias;
  Vector  lPl;
  UInt    startIndex;
  UInt    paraCount, rhsCount, obsCount;
  Double  sigma2;

public:
  NormalEquationRegularization(Config &config);

  UInt   rightHandSideCount() const override {return rhsCount;}
  UInt   parameterCount()     const override {return paraCount + startIndex;}
  void   parameterNames(std::vector<ParameterName> &/*names*/) const override {}
  void   init(MatrixDistributed &normals, UInt rhsCount) override;
  Bool   addNormalEquation(UInt rhsNo, const const_MatrixSlice &x, const const_MatrixSlice &Wz,
                           MatrixDistributed &normals, Matrix &n, Vector &lPl, UInt &obsCount) override;
  Vector contribution(MatrixDistributed &Cov) override;
  std::vector<Double> varianceComponentFactors() const override {return std::vector<Double>({sigma2});}
};

/***********************************************/

#endif
