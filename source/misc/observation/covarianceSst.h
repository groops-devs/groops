/***********************************************/
/**
* @file covarianceSst.h
*
* @brief Covariance matrix of satellite to satellite tracking observations.
*
* @author Torsten Mayer-Guerr
* @date 2010-07-18
*
*/
/***********************************************/

#ifndef __GROOPS_COVARIANCESST__
#define __GROOPS_COVARIANCESST__

// Latex documentation
#ifdef DOCSTRING_CovarianceSst
static const char *docstringCovarianceSst = R"(
\section{CovarianceSst}\label{covarianceSstType}
Provides arc wise covariance matrices for satellite-to-satellite observations SST).
The \configFile{inputfileCovarianceFunction}{matrix} provides a temporal covariance function.
From it the Toeplitz covariance matrix is constructed
\begin{equation}
  \M C = \begin{pmatrix}
    cov(t_0) & cov(t_1) & \cdots   &          &        &        \\
    cov(t_1) & cov(t_0) & cov(t_1) & \cdots   &        &        \\
    \cdots   & cov(t_1) & cov(t_0) & cov(t_1) & \cdots &        \\
             & \cdots   & \ddots   & \ddots   & \ddots & \cdots \\
  \end{pmatrix} \\
\end{equation}

The complete covariance matrix of an arc is given by
\begin{equation}
  \M C_{arc} = \sigma_0^2 \sigma_{arc}^2 \M C + \sigma_{S,arc}^2 \M S_{arc}+ \text{diag}(\sigma_1^2, \sigma_2^2, \ldots, \sigma_n^2)
\end{equation}
where \config{sigma}~$\sigma_0$ is an overall factor and the arc specific factors $\sigma_{arc}$
can be provided with \configFile{inputfileSigmasPerArc}{matrix}.
The second term describes general covariance matrices for each arc
\configFile{inputfileCovarianceMatrixArc}{matrix} together with the factors $\sigma_{S,arc}$ from \config{sigmasCovarianceMatrixArc}.
The last matrix can be used to downweight outliers in single epochs and will be added if
\configFile{inputfileSigmasPerEpoch}{instrument} is provided.

)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"

/**
* @defgroup CovarianceSstGroup CovarianceSst
* @brief Covariance function to decorrelate observation equations.
* @ingroup miscGroup
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class CovarianceSst;
typedef std::shared_ptr<CovarianceSst> CovarianceSstPtr;

/***** CLASS ***********************************/

/** @brief Covariance matrix of satellite to satellite observations. */
class CovarianceSst
{
  Double                sigma;
  Vector                sigmaArc;
  InstrumentFile        fileSigmaEpoch;
  Matrix                covFunction;
  std::vector<FileName> fileNamesCovarianceMatrix;
  Vector                covMatrixSigmas;

  static void testInput(const std::vector<Time> &times, const ObservationSigmaArc &sigmaEpoch, const_MatrixSliceRef covFunction, const_MatrixSliceRef W);

public:
  /// Constructor
  CovarianceSst() : sigma(1.0) {}

  /// Constructor
  CovarianceSst(Config &config, const std::string &name);

  /// Destructor.
  virtual ~CovarianceSst() {}

  /** @brief Full variance covariance matrix. */
  Matrix covariance(UInt arcNo, const std::vector<Time> &times);

  /** @brief Decorrelates observation equations.
  * The list of matrices, usually the observation vector and design matrices, are decorrelated. */
  void decorrelate(UInt arcNo, const std::vector<Time> &times, const std::list<MatrixSlice> &A);

  // --------------------

  /** @brief Decorrelates observation equations.
  * A list of matrices, usually the observation vector and design matrices, are decorrelated.
  * additionally, the Cholesky decomposition of the complete covariance matrix is returned.
  * @param times of the SST arcNo
  * @param sigmaArc variance factor for the toeplitz covariance function of this arc
  * @param sigmaEpoch Epoch-wise sigmas
  * @param covFunction Toeplitz covariance function
  * @param[in,out] Cov0 input: Empty or Non-Toeplitz covariance matrix for this arc. Output: Full covariance matrix.
  * @return reference to @a Cov0 */
  static MatrixSliceRef covariance(const std::vector<Time> &times, Double sigmaArc, const ObservationSigmaArc &sigmaEpoch,
                                   const_MatrixSliceRef covFunction, Matrix &Cov0);

  /** @brief Decorrelates observation equations.
  * A list of matrices, usually the observation vector and design matrices, are decorrelated.
  * additionally, the Cholesky decomposition of the complete covariance matrix is returned.
  * @param times of the SST arcNo
  * @param sigmaArc variance factor for the toeplitz covariance function of this arc
  * @param sigmaEpoch Epoch-wise sigmas
  * @param covFunction Toeplitz covariance function
  * @param[in,out] W input: Empty or Non-Toeplitz covariance matrix for this arc. Output: Cholesky decomposition of the full covariance matrix
  * @param[in,out] A Matrices to be decorrelated
  * @return reference to @a W */
  static MatrixSliceRef decorrelate(const std::vector<Time> &times, Double sigmaArc, const ObservationSigmaArc &sigmaEpoch,
                                    const_MatrixSliceRef covFunction, Matrix &W, const std::list<MatrixSlice> &A);

  /** @brief creates an derived instance of this class. */
  static CovarianceSstPtr create(Config &config, const std::string &name) {return std::make_shared<CovarianceSst>(config, name);}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class CovarianceSst.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates CovarianceSst */
template<> Bool readConfig(Config &config, const std::string &name, CovarianceSstPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif
