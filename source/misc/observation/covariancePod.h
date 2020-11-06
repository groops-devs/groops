/***********************************************/
/**
* @file covariancePod.h
*
* @brief Covariance matrix of kinematic orbits.
*
* @author Torsten Mayer-Guerr
* @date 2010-07-18
*
*/
/***********************************************/

#ifndef __GROOPS_COVARIANCEPOD__
#define __GROOPS_COVARIANCEPOD__

// Latex documentation
#ifdef DOCSTRING_CovariancePod
static const char *docstringCovariancePod = R"(
\section{CovariancePod}\label{covariancePodType}
Provides arc wise covariance matrices for precise orbit data.
Temporal correlations are modeled in the orbit system (along, cross, radial).
The \configFile{inputfileCovarianceFunction}{matrix} provides temporal covariance functions for each axis.
From the diagonal matrix for each time step
\begin{equation}
  Cov_{3\times3}(t) = \text{diag}(cov_x(t), cov_y(t), cov_z(t))
\end{equation}
the Toeplitz covariance matrix for an arc is constructed
\begin{equation}
  \M C = \begin{pmatrix}
    Cov(t_0) & Cov(t_1) & \cdots   &          &        &        \\
    Cov(t_1) & Cov(t_0) & Cov(t_1) & \cdots   &        &        \\
    \cdots   & Cov(t_1) & Cov(t_0) & Cov(t_1) & \cdots &        \\
             & \cdots   & \ddots   & \ddots   & \ddots & \cdots \\
  \end{pmatrix}
\end{equation}

The epoch wise $3\times3$ covariance matrices given by \configFile{inputfileCovariancePodEpoch}{instrument}
are eigen value decomposed
\begin{equation}
  \M C_{3\times3}(t_i) = \M Q \M\Lambda \M Q^T,
\end{equation}
where $\M Q$ is an orthgonal matrix and $\M\Lambda$ diagonal.
This used to split the covariances matrices
\begin{equation}
  \M C_{3\times3}(t_i) = \M D(t_i) \M D(t_i)^T = (\M Q \M\Lambda^{1/2} \M Q^T)(\M Q \M\Lambda^{1/2} \M Q^T)^T,
\end{equation}
and to compose a block diagonal matrix for an arc
\begin{equation}
  \M D = \text{diag}(\M D(t_1), \M D(t_2), \ldots, \M D(t_2)).
\end{equation}

The complete covariance matrix of an arc is given by
\begin{equation}
  \M C_{arc} = \sigma_0^2 \sigma_{arc}^2 \M D \M C \M D^T +
  \text{diag}(\sigma_1^2\M I_{3\times3}, \sigma_2^2\M I_{3\times3}, \ldots, \sigma_n^2\M I_{3\times3})
\end{equation}
where \config{sigma}~$\sigma_0$ is an overall factor
and the arc specific factors $\sigma_{arc}$ can be provided with \configFile{inputfileSigmasPerArc}{matrix}.
The last matrix can be used to downweight outliers in single epochs and will be added if
\configFile{inputfileSigmasPerEpoch}{instrument} is provided.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"

/**
* @defgroup CovariancePodGroup CovariancePod
* @brief Covariance function to decorrelate observation equations.
* @ingroup miscGroup
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class CovariancePod;
typedef std::shared_ptr<CovariancePod> CovariancePodPtr;

/***** CLASS ***********************************/

/** @brief Covariance matrix of kinematic orbit observations. */
class CovariancePod
{
  Double         sigma;
  Vector         sigmaArc;
  InstrumentFile fileSigmaEpoch;
  InstrumentFile fileCovPodEpoch;
  Matrix         covFunction;

  static void testInput(const OrbitArc &pod, const ObservationSigmaArc &sigmaEpoch, const Covariance3dArc &covPod, const_MatrixSliceRef covFunction);

public:
  /// Constructor
  CovariancePod() : sigma(1.0) {}

  /// Constructor
  CovariancePod(Config &config, const std::string &name);

  /// Destructor.
  virtual ~CovariancePod() {}

  /** @brief Full variance covariance matrix.
  * The matrix is related to observations given as x,y,z per epoch in CRF [m]. */
  Matrix covariance(UInt arcNo, const OrbitArc &pod);

  /** @brief Decorrelates observation equations.
  * The observations must be given as x,y,z per epoch in CRF [m].
  * The list of observation vector and design matrices are decorrelated. */
  void decorrelate(UInt arcNo, const OrbitArc &pod, const std::list<MatrixSlice> &A);

  /** @brief Decorrelates observation equations.
  * The observations must be given as x,y,z per epoch in CRF [m].
  * The list of observation vector and design matrices are decorrelated.
  * Computes additionally the Cholesky decomposition of the covariance matrix
  * (WARNING: Orbit rotation and epoch wise covariance matrix are not included in cholesky) */
  static void decorrelate(const OrbitArc &pod, Double sigmaArc,const ObservationSigmaArc &sigmaEpoch,
                          const Covariance3dArc &covPod, const_MatrixSliceRef covFunction, Matrix &W,
                          const std::list<MatrixSlice> &A);

  /** @brief creates an derived instance of this class. */
  static CovariancePodPtr create(Config &config, const std::string &name) {return std::make_shared<CovariancePod>(config, name);}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class CovariancePod.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates CovariancePod */
template<> Bool readConfig(Config &config, const std::string &name, CovariancePodPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif
