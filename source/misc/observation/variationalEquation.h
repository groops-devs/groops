/***********************************************/
/**
* @file variationalEquation.h
*
* @brief Variational equations.
*
* @author Torsten Mayer-Guerr
* @date 2012-05-29
*
*/
/***********************************************/

#ifndef __GROOPS_VARIATIONALEQUATION__
#define __GROOPS_VARIATIONALEQUATION__

#include "files/fileSatelliteModel.h"
#include "files/fileVariationalEquation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"

/***** CLASS ***********************************/

/** @brief Variational equations.
* @ingroup miscGroup */
class VariationalEquation
{
public:
  VariationalEquation();

  void init(SatelliteModelPtr satellite, const std::vector<ParametrizationGravityPtr> &parameterGravity,
            ParametrizationAccelerationPtr parameterAcceleration, const std::vector<Time> &stochasticPulse,
            EphemeridesPtr ephemerides, UInt integrationDegree);

  void setArc(VariationalEquationArc &arc);

  const VariationalEquationArc &arc() const {return arc_;}

  /** @brief recompute parameter count.
  * Must only be called if parameters in @a parameterGravity
  * or @a parameterAcceleration are changed (via @a setInterval).
  * Current arc becomes invalid and initArc must be called */
  void computeIndices();

  /** @brief Total number of unknown parameters.
  * Design matrix (PosDesign, VelDesign) includes
  * 1. gravity field parameters
  * 2. satellite parameters
  * 3. arc related parameters (inclusive 6 satellite state param). */
  UInt parameterCount() const {return parameterCount_;}

  /** @brief Number of gravity field parameters.
  * The parameterCount of Representation.*/
  UInt parameterCountGravity() const {return gravityCount;}

  /** @brief Number of unknown satellite related parameters.
  Includes stochastic pulses (x,y,z). */
  UInt parameterCountSatellite() const {return satCount;}

  /** @brief Number of unknown arc related parameters.
  * Includes satellite state (6 parameter). */
  UInt parameterCountSatelliteArc() const {return satArcCount;}

  /** @brief Names of all parameters, includes
  * 1. gravity field parameters)
  * 2. satellite parameters
  * 3. arc related parameters (inclusive 6 satellite state param).
  * The names are appended to @a name. */
  void parameterName(std::vector<ParameterName> &name) const;

  /** @brief Names of gravity field parameters.
  * The names are appended to @a name. */
  void parameterNameGravity(std::vector<ParameterName> &name) const;

  /** @brief Names of satellite related parameters.
  * The names are appended to @a name. */
  void parameterNameSatellite(std::vector<ParameterName> &name) const;

  /** @brief Names of arc related parameters.
  * The names are appended to @a name. */
  void parameterNameSatelliteArc(std::vector<ParameterName> &name) const;

  /** @brief Setup observation equations.
  * @param idEpoch from VariationalEquationArc.
  * @param[out] pos0 (3 x 1) reference position.
  * @param[out] PosDesign (3 x parameterCount()) Design matrix. */
  void position(UInt idEpoch, MatrixSliceRef pos0, MatrixSliceRef PosDesign);

  /** @brief Setup observation equations.
  * @param idEpoch from VariationalEquationArc.
  * @param[out] vel0 (3 x 1) reference velocity.
  * @param[out] VelDesign (3 x parameterCount()) design matrix. */
  void velocity(UInt idEpoch, MatrixSliceRef vel0, MatrixSliceRef VelDesign);

private:
  SatelliteModelPtr                      satellite;
  VariationalEquationArc                 arc_;
  EphemeridesPtr                         ephemerides;
  std::vector<ParametrizationGravityPtr> parametersGravity;
  ParametrizationAccelerationPtr         parameterAcceleration;
  std::vector<Time>                      timePulse;
  std::vector<Matrix>                    stochasticPulse;


  UInt  parameterCount_;
  UInt  idxGravity, gravityCount;
  UInt  idxSat,     satCount;
  UInt  idxPulse;
  UInt  idxSatArc,  satArcCount;

  // integration
  UInt                integrationDegree;
  std::vector<Vector> coeffIntegral;
  UInt                idCoeff;
  Double              deltaT;
  UInt                idEpochAlpha;
  Matrix              Alpha;
  UInt                idIntegrand;
  std::vector<Matrix> Integrand;

  void   initIntegration();
  void   computeAlpha(UInt idEpoch);
  Matrix computeIntegrand(UInt idEpoch);
};

/***********************************************/

#endif
