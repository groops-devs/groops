/***********************************************/
/**
* @file gnssParametrization.h
*
* @brief Parametrization of GNSS observations.
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSSPARAMETRIZATION__
#define __GROOPS_GNSSSPARAMETRIZATION__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrization = R"(
\section{GnssParametrization}\label{gnssParametrizationType}
This class defines the models and parameters of the linearized observation equations
for all phase and code measurements (see \program{GnssProcessing})
\begin{equation}\label{gnssParametrizationType:model}
  \M l - \M f(\M x_0) = \left.\frac{\partial \M f(\M x)}{\partial \M x}\right|_{\M x_0} \Delta\M x + \M\epsilon,
\end{equation}
where the left side is the observation vector minus the effects computed from the a priori models.
After each least squares adjustment
(see \configClass{GnssProcessing:processingStep:estimate}{gnssProcessingStepType:estimate})
the a priori parameters are updated
\begin{equation}\label{gnssParametrizationType:update}
  \M x_0 := \M x_0 + \Delta\hat{\M x}.
\end{equation}
The vector $\M x_0$ can be written with
\configClass{GnssProcessing:processingStep:writeAprioriSolution}{gnssProcessingStepType:writeAprioriSolution}.
Any \config{outputfiles} defined in the parametrizations are written with
\configClass{GnssProcessing:processingStep:writeResults}{gnssProcessingStepType:writeResults}.

Each parametrization (and possible constraint equations) has a \config{name} which enables
activating/deactivating the estimation of subsets of $\Delta\M x$ with
\configClass{GnssProcessing:processingStep:selectParametrizations}{gnssProcessingStepType:selectParametrizations}.
The a priori model $\M f(\M x_0)$ is unaffected and is always reduced.

The model for the different observation types can be described as
\begin{equation}\label{gnssParametrizationType:gnssFullModel}
\begin{split}
  f[\tau\nu a]_r^s(\M x) &= \text{geometry}(\M r_r^s) + \text{clock}^s(t) + \text{clock}_r(t) \\
               &+ \text{ionosphere}([\tau\nu],t,\M r_r^s) + \text{troposphere}(t,\M r_r^s) \\
               &+ \text{antenna}[\tau\nu a]^s  + \text{antenna}[\tau\nu a]_r \\
               &+ \text{bias}[\tau\nu a]^s + \text{bias}[\tau\nu a]_r
               + \lambda[L\nu] N[L\nu a]_r^s + \text{other}(\ldots) + \epsilon[\tau\nu a]_r^s
\end{split}
\end{equation}
The notation $[\tau\nu a]_r^s$ describes the
attribution to a signal type $\tau$ (i.e., C or L), frequency $\nu$,
signal attribute $a$ (e.g., C, W, Q, X), transmitting satellite $s$, and observing receiver $r$.
It follows the \href{https://files.igs.org/pub/data/format/rinex305.pdf}{RINEX 3 definition},
see \reference{GnssType}{gnssType}.

See also \program{GnssProcessing}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "parallel/matrixDistributed.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssNormalEquationInfo.h"


/**
* @defgroup gnssParametrizationGroup GnssParametrization
* @brief Parametrization of GNSS observations.
* @ingroup classesGroup
* The interface is given by @ref GnssParametrization. */
/// @{

/***** TYPES ***********************************/

class GnssParametrization;
class GnssParametrizationBase;
typedef std::shared_ptr<GnssParametrization> GnssParametrizationPtr;

/***** CLASS ***********************************/

/** @brief Parametrization of GNSS observations.
* An instance of this class can be created with @ref readConfig. */
class GnssParametrization
{
  std::vector<GnssParametrizationBase*> base;

public:
  /// Constructor.
  GnssParametrization(Config &config, const std::string &name);

 ~GnssParametrization();

  /** @brief init base on transmitters and receivers in @p gnss. */
  void init(Gnss *gnss, Parallel::CommunicatorPtr comm);

  /** @brief How many observations are needed to estimate parameters? */
  void requirements(GnssNormalEquationInfo &normalEquationInfo, std::vector<UInt> &transCount, std::vector<UInt> &transCountEpoch,
                    std::vector<UInt> &recvCount, std::vector<UInt> &recvCountEpoch);

  /** @brief Register parameters in @p normalEquationInfo. */
  void initParameter(GnssNormalEquationInfo &normalEquationInfo);

  /** @brief Correct observation equations/apply models. */
  void observationCorrections(GnssObservationEquation &eqn) const;

  /** @brief Total parameter vector used as priori Taylor point. */
  Vector aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo) const;

  /** @brief Design matrix for the basic observation equations @p eqn. */
  void designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const;

  /** @brief Add additional (pseudo-) observations equations to the normals for @p idEpoch. */
  void constraintsEpoch(const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;

  /** @brief Add additional (pseudo-) observations equations to the normals. */
  void constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;

  /** @brief Resolve ambiguities to integer. */
  Double ambiguityResolve(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                          const std::vector<Byte> &selectedTransmitters, const std::vector<Byte> &selectedReceivers,
                          const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger);

  /** @brief Update the values based on the passed estimated dx. */
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz);

  /** @brief Update the covariance information passed by @p covariance matrix. */
  void updateCovariance(const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance);

  /** @brief Write the output files defined in the parametrizations. */
  void writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const;

  /** @brief creates an derived instance of this class. */
  static GnssParametrizationPtr create(Config &config, const std::string &name) {return GnssParametrizationPtr(new GnssParametrization(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssParametrization.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and an class without points is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssParametrization */
template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class GnssParametrizationBase
{
public:
  virtual ~GnssParametrizationBase() {}

  static Bool isEnabled(const GnssNormalEquationInfo &normalEquationInfo, const std::string &name);

  virtual void   init(Gnss */*gnss*/, Parallel::CommunicatorPtr /*comm*/) {}
  virtual void   requirements(GnssNormalEquationInfo &/*normalEquationInfo*/, std::vector<UInt> &/*transCount*/, std::vector<UInt> &/*transCountEpoch*/,
                              std::vector<UInt> &/*recvCount*/, std::vector<UInt> &/*recvCountEpoch*/) {}
  virtual void   initParameter(GnssNormalEquationInfo &/*normalEquationInfo*/) {}
  virtual void   observationCorrections(GnssObservationEquation &/*eqn*/) const {}
  virtual void   aprioriParameter(const GnssNormalEquationInfo &/*normalEquationInfo*/, MatrixSliceRef /*x0*/) const {}
  virtual void   designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &/*eqn*/, GnssDesignMatrix &/*A*/) const {}
  virtual void   constraintsEpoch(const GnssNormalEquationInfo &/*normalEquationInfo*/, UInt /*idEpoch*/, MatrixDistributed &/*normals*/, std::vector<Matrix> &/*n*/, Double &/*lPl*/, UInt &/*obsCount*/) const {}
  virtual void   constraints(const GnssNormalEquationInfo &/*normalEquationInfo*/, MatrixDistributed &/*normals*/, std::vector<Matrix> &/*n*/, Double &/*lPl*/, UInt &/*obsCount*/) const {}
  virtual Double ambiguityResolve(const GnssNormalEquationInfo &/*normalEquationInfo*/, MatrixDistributed &/*normals*/, std::vector<Matrix> &/*n*/, Double &/*lPl*/, UInt &/*obsCount*/,
                                  const std::vector<Byte> &/*selectedTransmitters*/, const std::vector<Byte> &/*selectedReceivers*/,
                                  const std::function<Vector(const_MatrixSliceRef, MatrixSliceRef, const_MatrixSliceRef, Vector &, Double &)> &) {return 0;}
  virtual Double updateParameter(const GnssNormalEquationInfo &/*normalEquationInfo*/, const_MatrixSliceRef /*x*/, const_MatrixSliceRef /*Wz*/) {return 0;}
  virtual void   updateCovariance(const GnssNormalEquationInfo &/*normalEquationInfo*/, const MatrixDistributed &/*covariance*/) {}
  virtual void   writeResults(const GnssNormalEquationInfo &/*normalEquationInfo*/, const std::string &/*suffix*/) const {}
};

/***********************************************/

#endif
