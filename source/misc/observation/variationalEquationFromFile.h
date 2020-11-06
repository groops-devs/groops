/***********************************************/
/**
* @file variationalEquationFromFile.h
*
* @brief Variational equations.
*
* @author Torsten Mayer-Guerr
* @date 2014-03-24
*
*/
/***********************************************/

#ifndef __GROOPS_VARIATIONALEQUATIONFROMFILE__
#define __GROOPS_VARIATIONALEQUATIONFROMFILE__

#include "files/fileVariationalEquation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "misc/observation/variationalEquation.h"

/***** CLASS ***********************************/

/** @brief Variational equations.
* @ingroup miscGroup */
class VariationalEquationFromFile
{
public:
  /** @brief Linearized observation equations. */
  class ObservationEquation
  {
    public:
    std::vector<Time> times;

    Matrix pos0;      //!< Reference positions (Taylor point of linearization).
    Matrix PosDesign; //!< Design matrix, partial derivatives of positions with respect to parameters.

    Matrix vel0;      //!< Reference velocities (Taylor point of linearization).
    Matrix VelDesign; //!< Design matrix, partial derivatives of velocities with respect to parameters.

    std::vector<Rotary3d> rotSat;
    std::vector<Rotary3d> rotEarth;
  };

  /// Default Constructor.
  VariationalEquationFromFile() : _parameterCount(0), gravityCount(0) {}

  /// Destructor.
  ~VariationalEquationFromFile() {}

  /** @brief Open a new file.
  * Old open file is closed before. */
  void open(const FileName &fileName, ParametrizationGravityPtr parameterGravity, ParametrizationAccelerationPtr parameterAcceleration,
            const std::vector<Time> &stochasticPulse, EphemeridesPtr ephemerides, UInt integrationDegree);

  /** @brief Close the file. */
  void close();

  /** @brief Number of Arc in file. */
  UInt arcCount() const {return file.arcCount();}

  /** @brief Satellite model. */
  SatelliteModelPtr satellite() const {return file.satellite();}

  /** @brief recompute parameter count.
  * Must only be called if parameters in @a parameterGravity
  * or @a parameterAcceleration are changed.
  * Current arc becomes invalid. */
  void computeIndices();

  /** @brief Total number of unknown parameters.
  * Design matrix (PosDesign, VelDesign) includes
  * 1. gravity field parameters
  * 2. satellite parameters
  * 3. arc related parameters (inclusive 6 satellite state param). */
  UInt parameterCount() const {return _parameterCount;}

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
  * 1. gravity field parameters
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

  /** @brief Replace the star camera from the variational file with this data. */
  void replaceStarCamera(const UInt idEpochStart, const UInt idEpochEnd, std::vector<Rotary3d> rotSat);

  /** @brief Setup observation equations. */
  ObservationEquation integrateArc(Time timeStart, Time timeEnd, Bool computePosition, Bool computeVelocity, std::vector<Rotary3d> rotSat={});

  VariationalEquationArc refineVariationalEquationArc(UInt arcNo, const_MatrixSliceRef x);

private:
  FileVariationalEquation file;
  VariationalEquation     variationalEquation;
  VariationalEquationArc  arc;
  UInt                    arcNo;

  UInt  _parameterCount;
  UInt  idxGravity, gravityCount;
  UInt  idxSat,     satCount;
  UInt  idxSatArc,  satArcCount;

  void getArc(const Time &time);
};

/***********************************************/

#endif
