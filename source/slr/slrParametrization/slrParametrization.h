/***********************************************/
/**
* @file slrParametrization.h
*
* @brief Parametrization of SLR observations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRSPARAMETRIZATION__
#define __GROOPS_SLRSPARAMETRIZATION__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrization = R"(
\section{SlrParametrization}\label{slrParametrizationType}
This class defines the models and parameters of the linearized observation equations
for normal points (see \program{SlrProcessing})
\begin{equation}\label{slrParametrizationType:model}
  \M l - \M f(\M x_0) = \left.\frac{\partial \M f(\M x)}{\partial \M x}\right|_{\M x_0} \Delta\M x + \M\epsilon,
\end{equation}
where the left side is the observation vector minus the effects computed from the a priori models.
After each least squares adjustment
(see \configClass{SlrProcessing:processingStep:estimate}{slrProcessingStepType:estimate})
the a priori parameters are updated
\begin{equation}\label{slrParametrizationType:update}
  \M x_0 := \M x_0 + \Delta\hat{\M x}.
\end{equation}
The vector $\M x_0$ can be written with
\configClass{SlrProcessing:processingStep:writeAprioriSolution}{slrProcessingStepType:writeAprioriSolution}.
Any \config{outputfiles} defined in the parametrizations are written with
\configClass{SlrProcessing:processingStep:writeResults}{slrProcessingStepType:writeResults}.

Each parametrization (and possible constraint equations) has a \config{name} which enables
activating/deactivating the estimation of subsets of $\Delta\M x$ with
\configClass{SlrProcessing:processingStep:selectParametrizations}{slrProcessingStepType:selectParametrizations}.
The a priori model $\M f(\M x_0)$ is unaffected and is always reduced.

The model for the one way range observations between station $s$ and reflector $r$
can be described as
\begin{equation}\label{slrParametrizationType:slrFullModel}
\begin{split}
  f_s^r(\M x) &= \frac{1}{2}\left(\left\lVert \M r^r(t_{bounce})-\M r_s(t_{trans}) \right\rVert
                          + \left\lVert \M r_s(t_{recv})-\M r^r(t_{bounce}) \right\rVert\right)  \\
              &+ \text{troposphere}(t,\M r_{ss}^r)
               + \text{bias}^r + \text{bias}_s + \text{bias}_s^r + \text{other}(\ldots) + \epsilon_r^s
\end{split}
\end{equation}

See also \program{SlrProcessing}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "parallel/matrixDistributed.h"
#include "slr/slr.h"
#include "slr/slrObservation.h"
#include "slr/slrDesignMatrix.h"
#include "slr/slrNormalEquationInfo.h"


/**
* @defgroup slrParametrizationGroup SlrParametrization
* @brief Parametrization of SLR observations.
* @ingroup classesGroup
* The interface is given by @ref SlrParametrization. */
/// @{

/***** TYPES ***********************************/

class SlrParametrization;
class SlrParametrizationBase;
class SlrParametrizationGravityField;
typedef std::shared_ptr<SlrParametrization> SlrParametrizationPtr;

/***** CLASS ***********************************/

/** @brief Parametrization of SLR observations.
* An instance of this class can be created with @ref readConfig. */
class SlrParametrization
{
  std::vector<SlrParametrizationBase*> base;

public:
  /// Constructor.
  SlrParametrization(Config &config, const std::string &name);

 ~SlrParametrization();

  std::vector<const SlrParametrizationGravityField*> getParametrizationGravity() const;

  /** @brief init base on satellites and stations in @p slr. */
  void init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField);

  /** @brief Correct observation equations/apply models. */
  void observationCorrections(SlrObservationEquation &eqn) const;

  /** @brief Register parameters in @p normalEquationInfo. */
  void initParameter(SlrNormalEquationInfo &normalEquationInfo);

  /** @brief Total parameter vector used as priori Taylor point. */
  Vector aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo) const;

  /** @brief Design matrix for the basic observation equations @p eqn. */
  void designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const;

  /** @brief Add additional (pseudo-) observations equations to the normals. */
  void constraints(const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;

  /** @brief Update the values based on the passed estimated dx. */
  Double updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz);

  /** @brief Update the covariance information passed by @p covariance matrix. */
  void updateCovariance(const SlrNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance);

  /** @brief Write the output files defined in the parametrizations. */
  void writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const;

  /** @brief creates an derived instance of this class. */
  static SlrParametrizationPtr create(Config &config, const std::string &name) {return SlrParametrizationPtr(new SlrParametrization(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class SlrParametrization.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and an class without points is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates SlrParametrization */
template<> Bool readConfig(Config &config, const std::string &name, SlrParametrizationPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class SlrParametrizationBase
{
public:
  virtual ~SlrParametrizationBase() {}

  static Bool isEnabled(const SlrNormalEquationInfo &normalEquationInfo, const std::string &name);

  virtual void   getParametrizationGravity(std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/) const {}
  virtual void   init(Slr */*slr*/, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/) {}
  virtual void   observationCorrections(SlrObservationEquation &/*eqn*/) const {}
  virtual void   initParameter(SlrNormalEquationInfo &/*normalEquationInfo*/) {}
  virtual void   aprioriParameter(const SlrNormalEquationInfo &/*normalEquationInfo*/, MatrixSliceRef /*x0*/) const {}
  virtual void   designMatrix(const SlrNormalEquationInfo &/*normalEquationInfo*/, const SlrObservationEquation &/*eqn*/, SlrDesignMatrix &/*A*/) const {}
  virtual void   constraints(const SlrNormalEquationInfo &/*normalEquationInfo*/, MatrixDistributed &/*normals*/, std::vector<Matrix> &/*n*/, Double &/*lPl*/, UInt &/*obsCount*/) const {}
  virtual Double updateParameter(const SlrNormalEquationInfo &/*normalEquationInfo*/, const_MatrixSliceRef /*x*/, const_MatrixSliceRef /*Wz*/) {return 0;}
  virtual void   updateCovariance(const SlrNormalEquationInfo &/*normalEquationInfo*/, const MatrixDistributed &/*covariance*/) {}
  virtual void   writeResults(const SlrNormalEquationInfo &/*normalEquationInfo*/, const std::string &/*suffix*/) const {}
};

/***********************************************/

#endif
