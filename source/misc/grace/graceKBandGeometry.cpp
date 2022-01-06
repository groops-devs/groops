/***********************************************/
/**
* @file graceKBandGeometry.cpp
*
* @brief Compute quantities derived from the specific GRACE K-Band ranging instrument observation geometry.
*
* @author Matthias Ellmer
* @date 2017-11-13
*
*/
/***********************************************/

#include "base/import.h"
#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "graceKBandGeometry.h"

/***********************************************/

SatelliteTrackingArc GraceKBandGeometry::antennaCenterCorrection(const OrbitArc &orbit1, const OrbitArc &orbit2,
                                                                 const StarCameraArc &starCamera1, const StarCameraArc &starCamera2,
                                                                 const Vector3d &center1, const Vector3d &center2, UInt degree)
{
  try
  {
    Matrix A(orbit1.size(), 4);
    for(UInt i=0; i<A.rows(); i++)
    {
      const Vector3d u = (orbit2.at(i).position - orbit1.at(i).position);                                       // center of mass vector
      const Vector3d v = (starCamera2.at(i).rotary.rotate(center2) - starCamera1.at(i).rotary.rotate(center1)); // combined antenna offset vector
      A(i,1) = u.r() - (u+v).r();  // COM - ANT -> AOC
    }

    const std::vector<Time> times = orbit1.times();
    Polynomial p(times, degree);
    copy(p.derivative(   times, A.column(1)), A.column(2));
    copy(p.derivative2nd(times, A.column(1)), A.column(3));

    return Arc(times, A, Epoch::Type::SATELLITETRACKING);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GraceKBandGeometry::partialOfAntennaCenterCorrectionWrtRollPitchYaw(const OrbitArc &pos1, const OrbitArc &pos2, const StarCameraArc &rotSat1, const StarCameraArc &rotSat2,
                                                                         const Vector3d &center1, const Vector3d &center2,
                                                                         Matrix &SparseJacobian1, Matrix &SparseJacobian2)
{
  try
  {
    Arc::checkSynchronized({pos1, pos2, rotSat1, rotSat2});
    const UInt epochCount = pos1.size();

    // Really the partials are size epochCount x 3 * epochCount
    // However, for ranges, which we are computing here, all partials
    // outside of the current epoch are 0, so we save some memory
    // by only storing the fields that contain data.
    SparseJacobian1 = Matrix(epochCount, 3); ///< size epochCount x 3
    SparseJacobian2 = Matrix(epochCount, 3); ///< size epochCount x 3

    /*
     * Model:
     * S: SRF
     * C: CRF
     * N: Nominal (SRF is such that K-Frame is aligned to Line of Sight)
     * Roll/Pitch/Yaw are differential rotations between N and S
     * aoc = [R(S->C) c(S)] \cdot [u(C)]
     * We write the first rotary as a combination of two rotaries:
     * aoc = [R(N->C) R(S->N) c(S)] \cdot [u(C)]
     * R(N->C) is large rotations, exactly known (from orbit and satellite geometry) and error free
     * R(S->N) is small angles with observation errors (S is observed by SCA) and covariance from fusion
     * Re-write inner product
     * aoc = [R(N->C) R(S->N) c(S)] \cdot [u(C)]
     *     = [R(N->C) R(S->N) c(S)]^T [u(C)]
     *     = c(S)^T R(S->N)^T R(N->C)^T u(C)
     *     = c(S)^T R(S->N)^T [R(C->N) u(C)]
     * Now, the last term [R(N->C)^T u(C)] is the orientation of the baseline in the Nominal system.
     * Per definition, this is exactly the direction of the K-Band Phase center, so the equation is
     * aoc = c(S)^T R(N->S)^T c(S)/norm(c(S))
     * We can then take the derivative of this equation wrt roll/pitch/yaw in R(N->S)^T using the
     * chain rule
     * d aoc/ d rpy = (d aoc / d R) * (d R / d rpy)
     */

    // partial derivatives
    // -------------------
    auto daoc_dR = [](Vector3d c)
    {
      return 1./c.norm() * Matrix({{c.x()*c.x(), c.x()*c.y(), c.x()*c.z(),
                                    c.y()*c.x(), c.y()*c.y(), c.y()*c.z(),
                                    c.z()*c.x(), c.z()*c.y(), c.z()*c.z()}}); ///< size 1x9
    };

    // these are constant throughout the arc
    const auto daoc_dR1 = daoc_dR(center1); ///< size 1x9
    const auto daoc_dR2 = daoc_dR(center2); ///< size 1x9

    // this changes from epoch to epoch
    auto dR_drpy = [](Rotary3d R){
      Angle roll, pitch, yaw;
      R.cardan(roll, pitch, yaw);

      const Double cr = std::cos(roll);
      const Double sr = std::sin(roll);

      const Double cp = std::cos(pitch);
      const Double sp = std::sin(pitch);

      const Double cy = std::cos(yaw);
      const Double sy = std::sin(yaw);

      // roll
      Matrix dR_dr = {{0,  sp*cr*cy - sr*sy,  sp*sr*cy + sy*cr},
                      {0, -sp*sy*cr - sr*cy, -sp*sr*sy + cr*cy},
                      {0, -cp*cr,            -sr*cp           }};

      // pitch
      Matrix dR_dp = {{-sp*cy,  sr*cp*cy, -cp*cr*cy},
                      { sp*sy, -sr*sy*cp,  sy*cp*cr},
                      { cp,     sp*sr,    -sp*cr   }};
      // yaw
      Matrix dR_dy = {{-sy*cp, -sp*sr*sy + cr*cy, sp*sy*cr + sr*cy},
                      {-cp*cy, -sp*sr*cy - sy*cr, sp*cr*cy - sr*sy},
                      { 0,      0,                0               }};

      Matrix dR_drpy(9,3);
      copy(flatten(dR_dr), dR_drpy.column(0));
      copy(flatten(dR_dp), dR_drpy.column(1));
      copy(flatten(dR_dy), dR_drpy.column(2));
      return dR_drpy;
    };

    // Rotation from SRF -> K-Band frame
    const Rotary3d rotKFrame1 = inverse(Rotary3d(center1, Vector3d(0,1,0)));
    const Rotary3d rotKFrame2 = inverse(Rotary3d(center2, Vector3d(0,1,0)));

    // variance propagation
    // --------------------
    for (UInt i = 0; i < epochCount; i++)
    {
      // u is center of mass vector
      Vector3d u = pos2.at(i).position - pos1.at(i).position;
      u.normalize();

      // Nominal rotation CRF -> SRF. Nominal orientation of the SRF is so that the K-Frame is aligned to the line-of sight axis
      const Rotary3d rotSat1Nominal = Rotary3d( u, crossProduct( u, pos1.at(i).position)) * rotKFrame1;
      const Rotary3d rotSat2Nominal = Rotary3d(-u, crossProduct(-u, pos2.at(i).position)) * rotKFrame2;

      // Differential rotation (roll/pitch/yaw) about the nominal orientation (towards other satellite)
      const Rotary3d rotSat1observed = inverse(inverse(rotSat1Nominal) * rotSat1.at(i).rotary);
      const Rotary3d rotSat2observed = inverse(inverse(rotSat2Nominal) * rotSat2.at(i).rotary);

      matMult(1.0, daoc_dR1, dR_drpy(rotSat1observed), SparseJacobian1.row(i));
      matMult(1.0, daoc_dR2, dR_drpy(rotSat2observed), SparseJacobian2.row(i));
    } // for arc.epochCount
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix GraceKBandGeometry::sstResidual2RPYresidualProjector(const_MatrixSliceRef SparseJacobian, const_MatrixSliceRef CovarianceSca,
                                                            const std::vector<Time> &times, const UInt sstType, const UInt degree)
{
  try
  {
    // Each 3x3 block in CovarianceSca is Q_ij
    // Each 1x3 row in partialGrace is P_i
    // Each 3x1 block in I_ij := -1 * Q_ij * P_j.T
    const UInt epochCount = SparseJacobian.rows();
    Matrix A(3*epochCount, epochCount);
    for(UInt i=0; i<epochCount; i++)
      for(UInt j=0; j<epochCount; j++)
        matMult(1.0, CovarianceSca.slice(i*3,j*3,3,3), SparseJacobian.row(j).trans(), A.slice(3*i,j,3,1));

    // Derivative
    // ----------
    if(sstType == 0)
      return A;
    Polynomial p(times, degree);
    if(sstType == 1) // range rate
      return p.derivative(times, A.trans()).trans();
    if(sstType == 2) // range rangeAcceleration
      return p.derivative2nd(times, A.trans()).trans();
    throw(Exception("Unknown sst type"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix GraceKBandGeometry::variancePropagationStarCamera2SatelliteRanging(const_MatrixSliceRef SparseJacobian, const_MatrixSliceRef CovarianceSca,
                                                                          const std::vector<Time> &times, const UInt sstType, const UInt degree)
{
  try
  {
    // Variance propagation but the Jacobian of the range wrt. the quaternion elements is just the sparse representation
    const UInt epochCount = SparseJacobian.rows();

    // Each 3x3 block in CovarianceSca is Q_ij
    // Each 1x3 row in SparseJacobian is P_i
    // Variance propagation is C = P * Q * P.T
    // Each 1x1 block in C_ij := P_i * Q_ij * P_j.T
    Matrix Covariance(epochCount, Matrix::SYMMETRIC, Matrix::UPPER);
    for(UInt i=0; i<epochCount; i++)
      for(UInt j=i; j<epochCount; j++)
        Covariance(i,j) = inner(SparseJacobian.row(i).trans(), CovarianceSca.slice(i*3,j*3,3,3) * SparseJacobian.row(j).trans());

    // Derivative
    // ----------
    if(sstType == 0)
      return Covariance;

    Polynomial p(times, degree);
    fillSymmetric(Covariance);
    if(sstType == 1) // range rate
    {
      Covariance = p.derivative(times, Covariance);
      Covariance = p.derivative(times, Covariance.trans());
    }
    else if(sstType == 2) // range acceleration
    {
      Covariance = p.derivative2nd(times, Covariance);
      Covariance = p.derivative2nd(times, Covariance.trans());
    }
    else
      throw(Exception("Unknown sst type"));
    Covariance.setType(Matrix::SYMMETRIC);
    return Covariance;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GraceKBandGeometry::orientationResidual2rangingCorrection(const_MatrixSliceRef SparseJacobian, const_MatrixSliceRef deltaRpy,
                                                                 const std::vector<Time> &times, const UInt sstType, const UInt degree)
{
  try
  {
    // Observation corrections. The Jacobian of the range wrt. the quaternion elements is just the sparse representation
    // Correction is c = Jacobian * deltaRpy
    Vector correction(SparseJacobian.rows());
    for(UInt i=0; i<correction.rows(); i++)
      correction(i) = inner(SparseJacobian.row(i).trans(), deltaRpy.row(i*3,3));

    // Derivative
    // ----------
    if(sstType == 0)
      return correction;
    Polynomial p(times, degree);
    if(sstType == 1) // range rate
      return p.derivative(times, correction);
    if(sstType == 2) // range rangeAcceleration
      return p.derivative2nd(times, correction);
    throw(Exception("Unknown sst type"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
