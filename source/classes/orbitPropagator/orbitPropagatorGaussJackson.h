/***********************************************/
/**
* @file orbitPropagatorGaussJackson.h
*
* @brief Propagate a dynamic orbit using Gauss-Jackson method.
* @see orbitPropagator
*
* @author Matthias Ellmer
* @date 2017-02-07
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITPROPAGATORGAUSSJACKSON__
#define __GROOPS_ORBITPROPAGATORGAUSSJACKSON__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagatorGaussJackson = R"(
\subsection{GaussJackson}
This class implements the Gauss-Jackson multi-step predictor-corrector method to
propagate a satellite orbit under the influence of \configClass{Forces}{forcesType}.
Satellite is assumed to be oriented along-track. Implementation is based on [1].
[1] Berry, Matthew M., and Liam M. Healy. 2004. “Implementation of Gauss-Jackson Integration for Orbit Propagation.”
)";
#endif

/***********************************************/

#include "classes/orbitPropagator/orbitPropagator.h"

/***** CLASS ***********************************/

/** @brief Propagate orbit using GaussJackson method.
* @ingroup orbitPropagatorGroup
* @see orbitPropagator
* @see [1] Berry, Matthew M., and Liam M. Healy. 2004. “Implementation of Gauss-Jackson Integration for Orbit Propagation. */
class OrbitPropagatorGaussJackson : public OrbitPropagator
{
  OrbitPropagatorPtr warmup;     // Used for generation of warmup values.
  UInt   order;                  // of Gauss-Jackson method.
  UInt   maxCorrectorIterations; // For Corrector loop in predictor-corrector algorithm
  Double epsilon;                // For warmup and Corrector loop
  Matrix A;                      // Gauss-Jackson coefficients;
  Matrix B;                      // Summed Adams coefficients

  /** @brief Compute summation operator coefficients.
  * @param coefficients vector of propagator coefficients of length  @a order+1.
  * @param flipSign whether to flip the sign of the first element in the second to last row of coefficients.
  * @returns @a Matrix of size @a order+2, @a order+1, with the difference coefficients for the selected propagator.
  * The sign must be flipped to get results consistent with [1], Table 3,
  * computing the beta coefficients. This is not consistent with the
  * formulation given in eq.59, but leads to correct results. One must
  * assume that this is a copy & paste error in the paper. To get results
  * for Table 4 and eq. 63, do not flip the sign.
  * @see [1], eq. 59, eq. 63, Table 3, Table 4  */
  static Matrix differenceCoefficientsTable(const Vector coefficients, const Bool flipSign = FALSE);

  /** @brief Shift difference coefficients to the ordinate form.
  * This allows for the integration to be written as a linear combination of the
  * original function values.
  * @param coefficients Matrix of propagator coefficients.
  * @returns @a Matrix of same size as @a coefficients with the difference coefficients in ordinate form.
  * @see [1], Section Ordinate form, eq. 67  */
  static Matrix differenceCoefficientsToOrdinateForm(const_MatrixSliceRef coefficients);

  /** @brief Compute summed Adams propagator coefficients in ordinate form.
  * @param order of the propagator.
  * @returns @a Matrix of size order+1,order
  * @see [1], Table 5  */
  static Matrix summedAdamsCoefficients(UInt order);

  /** @brief Compute Gauss Jackson propagator coefficients in ordinate form.
  * @param order of the propagator.
  * @returns @a Matrix of size order+1,order
  * @see [1], Table 6  */
  static Matrix gaussJacksonCoefficients(UInt order);

  /** @brief Compute integration constants
  * @param orbit arc for which to compute the sums. Must be of length order+1.
  * @param dt sampling in seconds.
  * @param[out] firstSums first sum integration constants for the epochs in orbit.
  * @param[out] secondSums second sum integration constants for the epochs in orbit. */
  void computeSums(const OrbitArc &orbit, Double dt, std::vector<Vector3d> &firstSums, std::vector<Vector3d> &secondSums) const;

public:
  OrbitPropagatorGaussJackson(Config &config);

  OrbitArc integrateArc(const OrbitEpoch &startEpoch, const Time &sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite,
                        EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const override;
};

/***********************************************/

#endif /* __GROOPS_ORBITPROPAGATORGAUSSJACKSON__ */
