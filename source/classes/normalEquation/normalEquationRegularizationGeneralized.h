/***********************************************/
/**
* @file normalEquationRegularizationGeneralized.h
*
* @brief Regularization with sum of symmetric matrices.
* @see NormalEquation
*
* @author Andreas Kvas
* @date 2010-02-10
*
*/
/***********************************************/

#ifndef __GROOPS_NORMALEQUATIONREGULARIZATIONGENERALIZED__
#define __GROOPS_NORMALEQUATIONREGULARIZATIONGENERALIZED__

// Latex documentation
#ifdef DOCSTRING_NormalEquation
static const char *docstringNormalEquationRegularizationGeneralized = R"(
\subsection{RegularizationGeneralized}

Generalized regularization which is represented by the observation equation
\begin{equation}
\mathbf{x}_0 = \mathbf{I} \mathbf{x} + \mathbf{v}, \mathbf{v} \sim \mathcal{N}(0, \sum_k \sigma^2_k \mathbf{V}_k).
\end{equation}

There are no requirements for partial covariance matrices $\mathbf{V}_k$ except for them being symmetric.
The accumulated covariance matrix $\sum_k \sigma^2_k \mathbf{V}_k$ must be positive definite however.
The variance components $\sigma^2_k$ are estimated during the adjustment process and are assumed to be positive.
All \configFile{inputfilePartialCovarianceMatrix}{matrix} must be of same size
and must match the dimension of \configFile{inputfileBiasMatrix}{matrix}
(if provided, otherwise a zero vector of appropriate dimensions is created).

The parameter \config{aprioriSigma} determines the initial variance factor for the partial covariance matrices. Either one $\sigma_0$ can be
supplied or one for each $\mathbf{V}_k$.

The regularization matrix can be applied to a subset of parameters by adjusting \config{startIndex}.
)";
#endif

/***********************************************/

#include "classes/normalEquation/normalEquation.h"

/***** CLASS ***********************************/

/** @brief Regularization with sum of symmetric matrices.
* @ingroup normalEquationGroup
* @see NormalEquation */
class NormalEquationRegularizationGeneralized : public NormalEquationBase
{
  std::vector<FileName> fileNamesCovariance;
  Matrix                bias;
  UInt                  startIndex, startBlock;
  UInt                  paraCount, rhsCount;
  std::vector<Double>   sigma2;
  MatrixDistributed     Sigma; // = sum_i sigma2_i * V_i
  std::vector<MatrixDistributed> V;

public:
  NormalEquationRegularizationGeneralized(Config &config);

  UInt   rightHandSideCount() const override {return rhsCount;}
  UInt   parameterCount()     const override {return paraCount + startIndex;}
  void   parameterNames(std::vector<ParameterName> &/*names*/) const override {}
  void   init(MatrixDistributed &normals, UInt rhsCount) override;
  Bool   addNormalEquation(UInt rhsNo, const const_MatrixSlice &x, const const_MatrixSlice &Wz,
                           MatrixDistributed &normals, Matrix &n, Vector &lPl, UInt &obsCount) override;
  Vector contribution(MatrixDistributed &Cov) override;
  std::vector<Double> varianceComponentFactors() const override {return sigma2;}
};

/***********************************************/

#endif
