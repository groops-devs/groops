/***********************************************/
/**
* @file graceKBandGeometry.h
*
* @brief Compute quantities derived from the specific GRACE K-Band ranging instrument observation geometry.
*
* @author Matthias Ellmer
* @date 2017-11-13
*
*/
/***********************************************/

#ifndef __GROOPS_GRACEKBANDGEOMETRY__
#define __GROOPS_GRACEKBANDGEOMETRY__

#include "base/import.h"
#include "files/fileInstrument.h"

/***********************************************/

/** @brief GRACE K Band Geometry
* Compute quantities derived from the specific GRACE K-Band ranging instrument observation geometry.
* @ingroup miscGroup */
namespace GraceKBandGeometry
{
  /** @brief Computes the antenna offset correction.
  * @param orbit1 Orbit of the first satellite
  * @param orbit2 Orbit of the second satellite
  * @param starCamera1 rotation SRF->CRF
  * @param starCamera2 rotation SRF->CRF
  * @param center1 APC vector in SRF of first satellite
  * @param center2 APC vector in SRF of second satellite
  * @param degree for polynomial differentiation for range rates and range accelerations.
  * @return SatelliteTrackingArc with AOC in ranges, range rates, and range accelerations */
  SatelliteTrackingArc antennaCenterCorrection(const OrbitArc &orbit1, const OrbitArc &orbit2,
                                               const StarCameraArc &starCamera1, const StarCameraArc &starCamera2,
                                               const Vector3d &center1, const Vector3d &center2, UInt degree);

  /** @brief Computes partials of the Antenna Offset Correction with regard to roll/pitch/yaw of the satellites
  * @param pos1 Orbit of the first satellite
  * @param pos2 Orbit of the second satellite
  * @param rotSat1 rotation SRF->CRF
  * @param rotSat2 rotation SRF->CRF
  * @param center1 Taylor point for APC vector in SRF of first satellite
  * @param center2 Taylor point for APC vector in SRF of second satellite
  * @param[out] SparseJacobian1 epochCountx3 Matrix with partials of the aoc with regard to r/p/y of the first satellite
  * @param[out] SparseJacobian2 epochCountx3 Matrix with partials of the aoc with regard to r/p/y of the second satellite
  * @attention The Jacobians really have size [epochCount x 3 * epochCount]. However, for ranges, which we are computing
  * here, all partials outside of the current epoch are 0, so we save some memory by only storing the fields that contain data.  */
  void partialOfAntennaCenterCorrectionWrtRollPitchYaw(const OrbitArc &pos1, const OrbitArc &pos2, const StarCameraArc &rotSat1, const StarCameraArc &rotSat2,
                                                       const Vector3d &center1, const Vector3d &center2, Matrix &SparseJacobian1, Matrix &SparseJacobian2);

  /** @brief Computes the matrix to project SST residuals into residuals of roll/pitch/yaw angles
  * @see partialOfAntennaCenterCorrectionWrtRollPitchYaw() for the jacobian.
  * @param SparseJacobian epochCountx3 Matrix with partials of the aoc with regard to r/p/y of the first satellite
  * @param CovarianceSca  epochCount*3xepochCountx3 SCA covariance matrix in the SRF
  * @param times epochs of the arc
  * @param sstType 0: range, 1: rangeRate, 2: rangeAcceleration
  * @param degree for polynomial differentiation for range rates and range accelerations.
  * @note If arguments following dt are not given, ranges are returned.
  * @returns 3*epochCount x epochCount. To be pre-multiplied with WWe. Returned in range domain.  */
  Matrix sstResidual2RPYresidualProjector(const_MatrixSliceRef SparseJacobian, const_MatrixSliceRef CovarianceSca,
                                          const std::vector<Time> &times, const UInt sstType, const UInt degree);

  /** @brief Computes AOC covariance matrix for one satellite
  * @see partialOfAntennaCenterCorrectionWrtRollPitchYaw() for the jacobian.
  * @param SparseJacobian epochCountx3 Matrix with partials of the aoc with regard to r/p/y of the first satellite
  * @param CovarianceSca  epochCount*3xepochCountx3 SCA covariance matrix in the SRF
  * @param times epochs of the arc
  * @param sstType 0: range, 1: rangeRate, 2: rangeAcceleration
  * @param degree for polynomial differentiation for range rates and range accelerations.
  * @note If arguments after CovarianceSca are not given, ranges are returned.
  * @returns epochCount x epochCount AOC covariance matrix in the requested domain.  */
  Matrix variancePropagationStarCamera2SatelliteRanging(const_MatrixSliceRef SparseJacobian, const_MatrixSliceRef CovarianceSca,
                                                        const std::vector<Time> &times, const UInt sstType, const UInt degree);

  /** @brief Computes TLS observation correction due to orientation residuals
  * @see partialOfAntennaCenterCorrectionWrtRollPitchYaw() for the jacobian.
  * @param SparseJacobian epochCountx3 Matrix with partials of the aoc with regard to r/p/y of the first satellite
  * @param deltaRpy  epochCount*3 vector of updates to orientation
  * @param times epochs of the arc
  * @param sstType 0: range, 1: rangeRate, 2: rangeAcceleration
  * @param degree for polynomial differentiation for range rates and range accelerations.
  * @note If arguments after deltaRpy are not given, ranges are returned.
  * @returns epochCount vector with observation corrections in the requested domain. */
  Vector orientationResidual2rangingCorrection(const_MatrixSliceRef SparseJacobian, const_MatrixSliceRef deltaRpy,
                                               const std::vector<Time> &times, const UInt sstType, const UInt degree);
}

/***********************************************/

#endif // __GROOPS_GRACEKBANDGEOMETRY__
