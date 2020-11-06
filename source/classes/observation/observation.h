/***********************************************/
/**
* @file observation.h
*
* @brief Observation equations.
* Set up linearized observation equations (design matrix)
* to connect unknown parameters with observations.
* It is used in an least squares adjustement.
*
* @author Torsten Mayer-Guerr
* @date 2001-11-12
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATION__
#define __GROOPS_OBSERVATION__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservation = R"(
\section{Observation}\label{observationType}
This class set up the oberservation equations in linearized Gauss-Makoff model
\begin{equation}\label{gmm}
\M l  = \M A \M x + \M e\qquad\text{and}\qquad\mathcal{C}(\M e) = \sigma^2\M P^{-1}.
\end{equation}
The observations are divided into short data blocks which can computed independently
and so easily can be parallized. Usually this data blocks are short arcs of a
satellites orbit. In most cases the unknown parameter vector contains coefficients
of a gravity field parametrization given by \configClass{parametrizationGravity}{parametrizationGravityType}.
Additional parameters like instrument calibrations parameters are appended at the
end of the vector~$\M x$.
It is possible to give several observation vectors in one model.

The observations within each arc are decorrelated in the following way:
In a first step a cholesky decomposition of the covariance matrix is performed
\begin{equation}
\M P^{-1} = \M W^T\M W,
\end{equation}
where $\M W$ is an upper regular triangular matrix.
In a second step the transformation
\begin{equation}\label{dekorrelierung}
\bar{\M A} = \M W^{-T}\M A\qquad\text{and}\qquad \bar{\M l} = \M W^{-T}\M l
\end{equation}
gives an estimation from decorrelated observations with equal variance
\begin{equation}\label{normal.GMM}
\bar{\M l} = \bar{\M A} \M x + \bar{\M e}
\qquad\text{and}\qquad
\mathcal{C}(\bar{\M e})= \sigma^2 \M I.
\end{equation}
Usually the arc depending parameters are eliminated in the next step.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/parameterName.h"
#include "config/config.h"

/**
* @defgroup observationGroup Observation
* @brief Observation equations (design matrix).
* @ingroup classesGroup
* The interface is given by @ref Observation.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Observation;
typedef std::shared_ptr<Observation> ObservationPtr;

/***** CLASS ***********************************/

/** @brief Observation equations.
* Set up linearized observation equations (design matrix)
* to connect unknown parameters with observations.
* It is used in an least squares adjustement.
* An Instance of this class can be created by @ref readConfig. */
class Observation
{
public:
  virtual ~Observation() {}

  /** @brief Number of unknown parameters.
  * This is the column count of the design matrix @a A. */
  virtual UInt parameterCount() const=0;

  /** @brief Number of Gravityfield parameters.
  * The parameterCount of Representation.
  * (Without parameters added by this class, e.g. instrumental calibration parameters). */
  virtual UInt gravityParameterCount() const = 0;

  /** @brief Number of observations.
  * This is the column count of the observation vector @a l. */
  virtual UInt rightSideCount() const=0;

  /** @brief Number of observation blocks (arcs). */
  virtual UInt arcCount() const = 0;

  /** @brief Name of parameters.
  * The names are appended to @a name. */
  virtual void parameterName(std::vector<ParameterName> &name) const = 0;

  /** @brief Observation equations for an Arc.
  * The desgin matrix @a A contains the common parameter (mostly gravity field parameters).
  * The design matrix for arc specific parameters (e.g. bias, staellite state vector) is separated in Matrix @a B.
  * @param arcNo Index of the arc to be computedNummer [0, arcCount)
  * @param[out] l Observation vector (multiple vectors/columns possible)
  * @param[out] A Design matrix for common parameters.
  * @param[out] B Design matrix for arc related parameters.
  */
  virtual void observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B) = 0;

  /** @brief creates an derived instance of this class. */
  static ObservationPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Observation.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a observation is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] observation Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Observation */
template<> Bool readConfig(Config &config, const std::string &name, ObservationPtr &observation, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS_OBSERVATION__ */
