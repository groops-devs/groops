/***********************************************/
/**
* @file miscGriddedData.h
*
* @brief Misc functions for values on grid.
*
* @author Torsten Mayer-Guerr
* @date 2008-08-06
*
*/
/***********************************************/

#ifndef __GROOPS_MISCGRIDDEDDATA__
#define __GROOPS_MISCGRIDDEDDATA__

#include "base/vector3d.h"
#include "base/sphericalHarmonics.h"
#include "base/griddedData.h"
#include "parallel/parallel.h"
#include "classes/kernel/kernel.h"

/***********************************************/

/** @brief Functions for values on grid.
* @ingroup miscGroup */
namespace MiscGriddedData
{
  /** @brief Compute some statistics (rms, mean, max...). */
  void statistics(const std::vector<Double> &values, const std::vector<Double> &weights,
                  Double &rms, Double &avg, Double &vmin, Double &vmax, Double &mean);

  /** @brief Prints out some statistics (rms, mean, max...).
  * (only at master node, calls from clients will be ignored). */
  void printStatistics(const GriddedData &grid);

  /** @brief Prints out some statistics (rms, mean, max...).
  * (only at master node, calls from clients will be ignored). */
  void printStatistics(const GriddedDataRectangular &grid);

  /** @brief Generates functionals of a spherical harmonics expansion on a grid.
  * Must be called from every node in parallel computations.
  * @param harmonic spherical harmonics expansion
  * @param points harm is evaluated at these points (fast on rectangular grid)
  * @param kernel define the ouput functional.
  * @param comm   communicator for parallel computation.
  * @param timing start a loop timer for all grid points (only relevant for non-rectangular grids).
  * @return values at @a points (only valid at master). */
  std::vector<Double> synthesisSphericalHarmonics(const SphericalHarmonics &harmonic, const std::vector<Vector3d> &points, KernelPtr kernel, Parallel::CommunicatorPtr comm, Bool timing = TRUE);

  /** @brief Generates a linear functional for the synthesis of spherical harmonics coefficients on a grid.
  * This function generates a matrix A which represents the synthesis of a spherical harmonics vector x by matrix multiplication (y = Ax).
  * The columns of A start a degree zero are sorted degreewise.
  * @param maxDegree maximum expansion degree
  * @param GM geocentric gravitational constant
  * @param R reference radius
  * @param points evaluation points
  * @param kernel define the ouput functional.
  * @param isInterior flag whether interior of the sphere is requested
  * @return Matrix A. */
  Matrix synthesisSphericalHarmonicsMatrix(UInt maxDegree, Double GM, Double R, const std::vector<Vector3d> &points, KernelPtr kernel, Bool isInterior=FALSE);

} // namespace GriddedData

/***********************************************/

#endif /* __GROOPS_GRIDDEDDATA__ */
