/***********************************************/
/**
* @file integralEquation.h
*
* @brief Functions for the (short arc) integral equation approach.
*
* @author Torsten Mayer-Guerr
* @date 2005-02-24
*
*/
/***********************************************/

#ifndef __GROOPS_INTEGRALEQUATION__
#define __GROOPS_INTEGRALEQUATION__

#include "base/import.h"
#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Functions for the (short arc) integral equation approach.
* @ingroup misc */
class IntegralEquation
{
  UInt                integrationDegree, interpolationDegree;
  std::vector<Matrix> W;
  mutable Matrix      IntegrationPos;
  mutable Matrix      IntegrationVel;

public:

  /** @brief Linearized observation equations for an arc. */
  class Arc
  {
    public:
    Matrix vPos;         //!< Reference positions computed with a force model.
    Matrix VPos;         //!< Design matrix, partial derivatives of positions with respect to specific forces in TRF.
    Matrix VPosBoundary; //!< Design matrix, partial derivatives of positions with respect to satellite state (boundary values).

    Matrix vVel;         //!< Reference velocities computed with a force model.
    Matrix VVel;         //!< Design matrix, partial derivatives of velocities with respect to specific forces in TRF.
    Matrix VVelBoundary; //!< Design matrix, partial derivatives of velocities with respect to satellite state (boundary values).

    Matrix vAcc;         //!< Reference accelerations computed with a force model.
    Matrix VAcc;         //!< Design matrix, partial derivatives of accelerations with respect to specific forces in TRF.
    Matrix VAccBoundary; //!< Design matrix, partial derivatives of accelerations with respect to satellite state (boundary values).
  };

/** @brief Constructor.
* @param integrationDegree Degree of polynomial used for integration.
* @param interpolationDegree Degree of polynomial used for interpolation. */
IntegralEquation(UInt integrationDegree=7, UInt interpolationDegree=7);

/** @brief Initialization.
* @param integrationDegree Degree of polynomial used for integration.
* @param interpolationDegree Degree of polynomial used for interpolation. */
void init(UInt integrationDegree=7, UInt interpolationDegree=7);

/** @brief Observation equations for positions, velocities, and accelerations.
* @param orbit Approximate positions.
* @param rotEarth Rotation matricies CRF -> TRF at each epoch.
* @param gradientfield Gravity field for indirect effect (gravity gradients)
* @param g Approcimate specific forces along the orbit (x,y,z per position) in CRF; Multiple columns for multiple right hand sides are possible.
* @param accCalibration Observation equations for emprical accealrations (Allowed to be empty).
* @param computeVelocity Compute observation equations for velocities additionally to positions?
* @param computeAcceleration Compute observation equations for accelerations additionally to positions? */
Arc integrateArc(const OrbitArc &orbit, const std::vector<Rotary3d> &rotEarth, GravityfieldPtr gradientfield,
                 const const_MatrixSlice &g, const const_MatrixSlice &accCalibration, Bool computeVelocity=FALSE, Bool computeAcceleration=FALSE) const;

/** @brief Interpolate arc at observation (pod) times.
* @param pod         kinematic positions as pseudo observations (for each right hand side).
* @param orbit       positions that be used to compute the arc with integrateArc.
* @param arc         computed with integrateArc.
* @param[out] l      reduced observation vector.
* @param[out] VPos   design matrix of the unknown forces.
* @param[out] VPosBoundary design matrix of the unknown satellite state. */
void interpolateArc(const std::vector<OrbitArc> &pod, const OrbitArc &orbit, const Arc &arc, Matrix &l, Matrix &VPos, Matrix &VPosBoundary) const;

/** @brief Partial derivatives of positions (in CRF) with respect to specific forces (in CRF).
* Integration of forces with the integral kernel.
* @param T Arc length in seconds.
* @param posCount Number of positions (uniformly distributed). */
Matrix integrationMatrixPosition(Double T, UInt posCount) const;

/** @brief Partial derivatives of velocities (in CRF) with respect to specific forces (in CRF).
* Integration of forces with the integral kernel.
* @param T Arc length in seconds.
* @param posCount Number of positions (uniformly distributed). */
Matrix integrationMatrixVelocity(Double T, UInt posCount) const;
};

/***********************************************/

#endif
