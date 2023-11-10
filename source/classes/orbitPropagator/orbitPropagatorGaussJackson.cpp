/***********************************************/
/**
* @file orbitPropagatorGaussJackson.cpp
*
* @brief Propagate a dynamic orbit using Gauss-Jackson method.
* @see orbitPropagator
*
* @author Matthias Ellmer
* @date 2017-02-07
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "classes/orbitPropagator/orbitPropagatorAdamsBashforthMoulton.h"
#include "classes/orbitPropagator/orbitPropagatorStoermerCowell.h"
#include "classes/orbitPropagator/orbitPropagatorGaussJackson.h"

/***********************************************/

OrbitPropagatorGaussJackson::OrbitPropagatorGaussJackson(Config &config)
{
  try
  {
    readConfig(config, "order",               order,                  Config::MUSTSET,  "8",     "of Gauss-Jackson method.");
    readConfig(config, "warmup",              warmup,                 Config::MUSTSET,  "rungeKutta4", "");
    readConfig(config, "correctorIterations", maxCorrectorIterations, Config::DEFAULT,  "10",    "Maximum number of iterations to run the corrector step for.");
    readConfig(config, "epsilon",             epsilon,                Config::DEFAULT,  "1e-15", "Convergence criteria for position, velocity, and acceleration tests.");
    if (isCreateSchema(config)) return;

    if(order%2)
      throw(Exception("Gauss-Jackson mid-corrector order must be even, is: "+order%"%i"s));

    // The Gauss-Jackson integrator is the summed form of the
    // Stoermer-Cowell integrator, which is in turn based on the Adams-Bashforth coefficients.
    A = gaussJacksonCoefficients(order);
    B = summedAdamsCoefficients(order);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

OrbitArc OrbitPropagatorGaussJackson::integrateArc(const OrbitEpoch &startEpoch, const Time &sampling, UInt posCount, ForcesPtr forces,
                                                   SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const
{
  try
  {
    // Compute warmup
    OrbitArc orbit = flip(warmup->integrateArc(startEpoch, -sampling, order+1, forces, satellite, earthRotation, ephemerides, FALSE));

    // Compute integration constants
    // -----------------------------
    const Double dt = sampling.seconds();
    std::vector<Vector3d> firstSums, secondSums;
    computeSums(orbit, dt, firstSums, secondSums);
    Vector3d firstSum  = firstSums.back();
    Vector3d secondSum = secondSums.back();

    // Propagate rest of orbit arc
    // ---------------------------
    Single::forEach(posCount-1, [&](UInt k)
    {
      secondSum += firstSum + 0.5*orbit.at(k+order).acceleration;

      // Predict
      // -------
      OrbitEpoch epoch;
      epoch.time     = startEpoch.time + (k+1) * sampling;
      epoch.position = dt*dt * secondSum; // eq. 88;
      for(UInt j=0; j<order+1; j++)
        epoch.position += dt*dt * A(order+1, j) * orbit.at(k+j).acceleration;
      epoch.velocity = dt * (firstSum + 0.5 * orbit.at(k+order).acceleration);
      for(UInt j=0; j<order+1; j++)
        epoch.velocity += dt * B(order+1, j) * orbit.at(k+j).acceleration;
      epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);
      orbit.push_back(epoch);

      // Correct iteratively
      // -------------------
      // Precompute the accelerations from epochs up until the predicted one
      Vector3d summedAccelerationsPosition;
      for(UInt j=0; j<order; j++)
        summedAccelerationsPosition += A(order, j) * orbit.at(k+j+1).acceleration;
      Vector3d summedAccelerationsVelocity;
      for(UInt j=0; j<order; j++)
        summedAccelerationsVelocity += B(order, j) * orbit.at(k+j+1).acceleration;

      Vector3d firstSumOld = firstSum;
      for(UInt iter=0; iter<maxCorrectorIterations; iter++)
      {
        firstSum = firstSumOld + 0.5 * (orbit.at(k+order).acceleration + epoch.acceleration);

        // Compute updated state, including accelerations due to current predicted/corrected epoch
        const OrbitEpoch epochOld = epoch;
        epoch.position     = dt * dt * (secondSum + summedAccelerationsPosition + A(order, order) * epoch.acceleration);
        epoch.velocity     = dt      * (firstSum  + summedAccelerationsVelocity + B(order, order) * epoch.acceleration);
        epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);
        if(((epoch.position-epochOld.position).norm()<=epsilon) && ((epoch.velocity-epochOld.velocity).norm()<=epsilon))
          break;
      }

      orbit.at(order+k+1) = epoch;
    }, timing);

    // remove warmup
    orbit.remove(0, order);

    return orbit;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// eq. 59 + 63
// Coefficients of a certain degree are combinations of the coefficients of the next-higher degree.
// Start with the last row of the coefficient matrix and work backwards.
Matrix OrbitPropagatorGaussJackson::differenceCoefficientsTable(const Vector coefficients, const Bool flipSign)
{
  const UInt N = coefficients.rows();
  Matrix mat(N+1, N);
  copy(coefficients.trans(), mat.row(N));

  for(UInt i=N; i-->0;)
  {
    mat(i, 0) = mat(i+1, 0);
    if((i==N-1) && flipSign)
      mat(i, 0) *= -1;

    for(UInt j=1; j<N; j++)
      mat(i, j) = mat(i+1, j) - mat(i+1, j-1);
  }
  return mat;
}

/***********************************************/

// eq. 67
// And flip it. m is the current point with larger positive integers for m being the m-th past point.
// So we want it the other way around.
Matrix OrbitPropagatorGaussJackson::differenceCoefficientsToOrdinateForm(const_MatrixSliceRef coefficients)
{
  Matrix binomial(coefficients.rows(), coefficients.columns(), 1.);
  for(UInt i=2; i<binomial.rows(); i++)
    for(UInt k=1; k<i; k++)
      binomial(i,k) = binomial(i-1,k-1) + binomial(i-1,k);

  Matrix xi(coefficients.rows(), coefficients.columns());
  for(UInt n=0; n<xi.rows(); n++)
    for (UInt m=0; m<xi.columns(); m++)
      for(UInt i=m; i<xi.columns(); i++)
        xi(n, xi.columns()-1-m) += std::pow(-1, m) * coefficients(n, i) * binomial(i, m);
  return xi;
}

/***********************************************/

// Summed Adams coefficients in ordinate form, Table 5
Matrix OrbitPropagatorGaussJackson::summedAdamsCoefficients(UInt order)
{
  Vector summedAdams                 = OrbitPropagatorAdamsBashforthMoulton::factorsBashforth(order+1).row(1, order+1); // Velocity component
  Matrix summedAdamsCoefficientTable = differenceCoefficientsTable(summedAdams);
  Matrix saCoefficientsOrdinate      = differenceCoefficientsToOrdinateForm(summedAdamsCoefficientTable);
  for(UInt i=0; i<saCoefficientsOrdinate.columns(); i++) // removeHalfAccelerationTerm
    saCoefficientsOrdinate(i, i) -= 0.5;
  return saCoefficientsOrdinate; // Array a_jk in [1], Table 5
}

/***********************************************/

// Gauss-Jackson coefficients in ordinate form, Table 6
Matrix OrbitPropagatorGaussJackson::gaussJacksonCoefficients(UInt order)
{
  Vector stoermerCoefficients = OrbitPropagatorStoermerCowell::stoermerCoefficients(order+2).row(2, order+1); // Position component
  return differenceCoefficientsToOrdinateForm(differenceCoefficientsTable(stoermerCoefficients)); // Array b_jk in [1], Table 6
}

/***********************************************/

void OrbitPropagatorGaussJackson::computeSums(const OrbitArc &orbit, Double dt, std::vector<Vector3d> &firstSums, std::vector<Vector3d> &secondSums) const
{
  try
  {
    const UInt N      = order+1;
    const UInt index0 = order/2;

    firstSums  = std::vector<Vector3d>(N); // eq. 75
    secondSums = std::vector<Vector3d>(N); // eq. 86

    // First integration constant, see: eq. 73
    Vector3d C1p = orbit.at(index0).velocity/dt;
    for(UInt k=0; k<N; k++)
      C1p -= B(index0, k) * orbit.at(k).acceleration;
    Vector3d C1 = C1p + orbit.at(index0).acceleration/2;

    // Second integration constant, see: eq. 85
    Vector3d C2 = C1 + orbit.at(index0).position/(dt*dt);
    for(UInt k=0; k<N; k++)
      C2 -= A(index0, k) * orbit.at(k).acceleration;

    // Initial values for first and second sums.
    // -----------------------------------------
    // At midpoint of arc / start epoch
    firstSums.at(index0)  = C1p;
    secondSums.at(index0) = C2 - C1;

    // Epochs around midpoint of arc
    for(UInt n=1; n<=order/2; n++)
    {
      // eq. 75
      firstSums.at(index0+n) = firstSums.at(index0-1+n) + 0.5 * (orbit.at(index0-1+n).acceleration + orbit.at(index0+n).acceleration); // n > 0 case
      firstSums.at(index0-n) = firstSums.at(index0+1-n) - 0.5 * (orbit.at(index0+1-n).acceleration + orbit.at(index0-n).acceleration); // n < 0 case

      // eq. 86
      secondSums.at(index0+n) = secondSums.at(index0-1+n) + firstSums.at(index0-1+n) + 0.5 * orbit.at(index0-1+n).acceleration; // n > 0 case
      secondSums.at(index0-n) = secondSums.at(index0+1-n) - firstSums.at(index0+1-n) + 0.5 * orbit.at(index0+1-n).acceleration; // n < 0 case
    }
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
