/***********************************************/
/**
* @file kalmanProcessing.h
*
* @brief Miscellaneous functions for Kalman filter applications
*
* @author Andreas Kvas
* @date 2016-12-26
*
*/
/***********************************************/

#ifndef __GROOPS_KALMANPROCESSING__
#define __GROOPS_KALMANPROCESSING__

// Latex documentation
#ifdef DOCSTRING_AutoregressiveModelSequence
static const char *docstringAutoregressiveModelSequence = R"(
\section{AutoregressiveModelSequence}\label{autoregressiveModelSequenceType}
Represents a sequence of multivariate autoregressive (AR) models with increasing order $p$.
The AR models should be stored as \file{matrix file}{matrix} in the \reference{GROOPS definition of
AR models}{fundamentals.autoregressiveModel}.
The required AR models can be computed with \program{CovarianceMatrix2AutoregressiveModel},
and passed to this class through
\config{inputfileAutoregressiveModel} in increasing order.

The main purpose of AutoregressiveModelSequence is to use AR models of the form
\begin{equation}
  \label{eq:ar-model}
  \mathbf{y}_e(t_i) = \sum_{k=1}^p \mathbf{\Phi}^{(p)}_k\mathbf{y}_e(t_{i-k}) + \mathbf{w}(t_i),
  \hspace{5pt} \mathbf{w}(t_i) \sim \mathcal{N}(0, \mathbf{\Sigma}^{(p)}_\mathbf{w}),
\end{equation}
to create pseudo-observation equations
\begin{equation}
  \label{eq:pseudo-observations-transformed}
  0 = \bar{\mathbf{\Phi}} \Delta\mathbf{y} + \bar{\mathbf{w}}, \hspace{5pt} \bar{\mathbf{w}} \sim
  \mathcal{N}(0, \bar{\mathbf{\Sigma}}_{\bar{\mathbf{w}}}),
\end{equation}
with
\begin{equation}
  \label{eq:pseudo-observations-ar}
  \bar{\mathbf{\Phi}} =
  \begin{bmatrix}
    \mathbf{I} & & & & & \\
    -\mathbf{\Phi}^{(1)}_1 & \mathbf{I} & & & &  \\
    -\mathbf{\Phi}^{(2)}_2 & -\mathbf{\Phi}^{(2)}_1 & \mathbf{I} & & & \\
    -\mathbf{\Phi}^{(3)}_3 & -\mathbf{\Phi}^ {(3)}_2 & -\mathbf{\Phi}^ {(3)}_1 & \mathbf{I} & &  \\
    & -\mathbf{\Phi}^{(3)}_3 & -\mathbf{\Phi}^ {(3)}_2 & -\mathbf{\Phi}^ {(3)}_1 & \mathbf{I} &  \\
    & & \ddots & \ddots & \ddots & \ddots  \\
  \end{bmatrix},
  \hspace{15pt}
  \bar{\mathbf{\Sigma}}_{\bar{\mathbf{w}}} =
  \bar{\mathbf{\Sigma}}_{\bar{\mathbf{w}}} =
  \begin{bmatrix}
    \mathbf{\Sigma}^{(0)}_{\mathbf{w}} & & & & & \\
    & \mathbf{\Sigma}^{(1)}_{\mathbf{w}} & & & & \\
    & & \mathbf{\Sigma}^{(2)}_{\mathbf{w}} & & & \\
    & & & \mathbf{\Sigma}^{(3)}_{\mathbf{w}} & & \\
    & & & & \mathbf{\Sigma}^{(3)}_{\mathbf{w}} &  \\
    & & & & & \ddots \\
  \end{bmatrix}.
\end{equation}
used to constrain high-frequency temporal gravity field variations (see
\program{KalmanSmootherLeastSquares}, \program{NormalsBuildShortTimeStaticLongTime},
\program{PreprocessingSst}).

The corresponding normal equation coefficient matrix is given by
\begin{equation}
  \label{eq:ar-normals}
  \bar{\mathbf{\Phi}}^T\bar{\mathbf{\Sigma}}^{-1}_{\bar{\mathbf{w}}}\bar{\mathbf{\Phi}}
\end{equation}
and if all AR models are estimated from the same sample its inverse is a block-Toeplitz covariance matrix
\begin{equation}
  (\mathbf{\Sigma}_{\mathbf{y}_m})_{ij} =
  \begin{cases}
 \mathbf{\Sigma}(|j-i|) & \text{for } i \leq j \\
 \mathbf{\Sigma}(|j-i|))^T & \text{otherwise}
 \end{cases},
\end{equation}
which can be computed using \program{AutoregressiveModel2CovarianceMatrix}.

A detailed description with applications can be found in:
Kvas, A., Mayer-Gürr, T. GRACE gravity field recovery with background model uncertainties.
J Geod 93, 2543–2552 (2019). \url{https://doi.org/10.1007/s00190-019-01314-1}
)";
#endif

/***********************************************/

#include "base/matrix.h"
#include "config/config.h"

/***** CLASS ***********************************/

class AutoregressiveModel
{
private:
  std::vector< std::vector<Matrix> > N; // normalEquations
  Bool hasNormals;
  Matrix _model;

  /** @brief Normal equations from a multivariate AR Model
  * The AR model of order p \f$ x_t = B.at(0) x_{t-1} + B.at(1) x_{t-2} + ...
  * + B.at(p-1) + w \f$ with \f$ w \sim Q \f$ represents a constraint in
  * time domain. The corresponding observation equation reads:
  * \f$ 0 = [ B.at(p-1) B.at(p-2) ... B.at(0) -I] + w \f$ with \f$ w \sim Q\f$.
  * This function computes the upper triangle of the resulting normal equation. */
  void computeNormalEquation();

public:
  /** @brief Constructor. */
  AutoregressiveModel(const_MatrixSliceRef B) : hasNormals(FALSE), _model(B) {}

  /** @brief Return the normal equations for a single epoch as block matrix. */
  std::vector< std::vector<Matrix> > normalEquation() { if(!hasNormals) this->computeNormalEquation(); return N; }

  /** @brief Return the normal equation block (row, column) in a matrix with blockCount blocks. */
  Matrix distributedNormalsBlock(UInt blockCount, UInt row, UInt column);

  /** @brief Return the normal equation block (row, column) in a matrix with blockCount blocks. */
  void distributedNormalsBlock(UInt blockCount, UInt row, UInt column, Matrix &N);

  /** @brief Dimension (number of parameters) of the model. */
  UInt dimension() const { return _model.rows(); }

  /** @brief Model order (number of coefficients). */
  UInt order() const { return _model.columns()/_model.rows()-1; }

  /** @brief Return pseudo observation equation of the AR(p) model */
  Matrix pseudoObservationEquation() const { return _model; }

  /** @brief Assign a new prediction error covariance to the AR model. */
  void updatePredictionErrorCovariance(const Matrix &Q);

  /** @brief Compute the process residuals (e.g. Bx_t - x_t+1) from a solution vector. */
  Vector decorrelate(const Vector& x) const;

  /** @brief Return AR coefficients as std::vector */
  std::vector<Matrix> coefficients() const;

  /** @brief Return autocovariance of process noise as Matrix */
  Matrix whiteNoiseCovariance() const;

  /** @brief Return the AR(1) representatation of the AR model */
  void orderOneRepresentation(Matrix &B, Matrix &Q) const;

  /** @brief list of row, column pairs for each distributed normals block that will be filled. */
  std::vector<std::pair<UInt, UInt>> distributedNormalsBlockIndex(UInt epochCount) const;

  /** @brief Change the state representation of the AR model. */
  void rotate(const_MatrixSliceRef F);
};

/***** CLASS ***********************************/

class AutoregressiveModelSequence;
typedef std::shared_ptr<AutoregressiveModelSequence> AutoregressiveModelSequencePtr;

class AutoregressiveModelSequence
{
private:
  std::vector<AutoregressiveModel> _arModels;

public:
  AutoregressiveModelSequence(Config &config);

  /** @brief Maxmimum order of the AR models in the sequence. */
  UInt maximumOrder() const { return _arModels.back().order(); }

  /** @brief Dimension (number of parameters) of the model. */
  UInt dimension() const { return _arModels.front().dimension(); }

  /** @brief Return the normal equation block (row, column) in a matrix with blockCount blocks. */
  Matrix distributedNormalsBlock(UInt blockCount, UInt row, UInt column);

  /** @brief Return the normal equation block (row, column) in a matrix with blockCount blocks. */
  void distributedNormalsBlock(UInt blockCount, UInt row, UInt column, Matrix &N);

  /** @brief list of row, column pairs for each distributed normals block that will be filled. */
  std::vector<std::pair<UInt, UInt>> distributedNormalsBlockIndex(UInt epochCount) const { return _arModels.back().distributedNormalsBlockIndex(epochCount); }

  /** @brief Return blocked normal equations for each AR model as vector */
  std::vector< std::vector< std::vector<Matrix> > > normalEquationSequence();

  /** @brief creates an derived instance of this class. */
  static AutoregressiveModelSequencePtr create(Config &config, const std::string &name);
};

/** @brief Creates an instance of the class AutoregressiveModelSequence. */
template<> Bool readConfig(Config &config, const std::string &name, AutoregressiveModelSequence &arModelSequence, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***********************************************/

#endif /* __GROOPS_KALMANPROCESSING__ */
