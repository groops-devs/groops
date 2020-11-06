/***********************************************/
/**
* @file rotary3d.h
*
* @brief Rotations in 3d space.
*
* @author Torsten Mayer-Guerr
* @date 2001-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_ROTARY3D__
#define __GROOPS_ROTARY3D__

#include "base/importStd.h"
#include "base/matrix.h"

/** @addtogroup vector3dGroup */
/// @{

class Angle;
class Vector3d;
class Tensor3d;
class Transform3d;
class Ellipsoid;
class SphericalHarmonics;

/***** CLASS ***********************************/

/** @brief Rotations in 3d space.
* Rotation of @a Vector3d and @a Tensor3d. */
class Rotary3d
{
  std::array<std::array<Double,3>,3> field;

public:
  Rotary3d() : field{{{1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}} {} //!< Default constructor (Unitary matrix).

  /// Contructor from std::array.
  Rotary3d(const std::array<std::array<Double,3>,3> &x) : field(x) {}

  /// Contructor from quaternions (q0, qx, qy, qz) or rotary matrix (3x3).
  explicit Rotary3d(const_MatrixSliceRef q);

  /** @brief Rotation of a local system.
  * Rotational matrix from local system (defined by the given axis x,y,z)
  * to the global system. The axis will be normalized and orthogonalized before.
  * The z-axis is defined orthogonal to x and y.
  * Zhe y-axis is defined orthogonal to the new z-axis and x. */
  Rotary3d(Vector3d x, Vector3d y);

  /// Quaternions (q0, qx, qy, qz).
  Vector quaternion() const;

  /// Cast to Matrix.
  Matrix matrix() const;

  /** @brief Rotate a vector.
  * \f[ y = D \cdot v \f]. */
  Vector3d rotate(const Vector3d &v) const;
  /** @copydoc rotate */
  Vector3d transform(const Vector3d &v) const;

  /** @brief Inverse rotation of a vector.
  * \f[ y = D^T \cdot v \f]. */
  Vector3d inverseRotate(const Vector3d &v) const;
  /** @copydoc inverseRotate */
  Vector3d inverseTransform(const Vector3d &v) const;

  /** @brief Rotate a tensor.
  * Both sides of the dyadic product are rotated.
  * \f[ y = D \cdot t \cdot D^T \f]. */
  Tensor3d rotate(const Tensor3d &t) const;
  /** @copydoc rotate */
  Tensor3d transform(const Tensor3d &t) const;

  /** @brief Inverse rotation a tensor.
  * Both sides of the dyadic product are rotated.
  * \f[ y = D^T \cdot t \cdot D \f]. */
  Tensor3d inverseRotate(const Tensor3d &t) const;
  /** @copydoc inverseRotate */
  Tensor3d inverseTransform(const Tensor3d &t) const;

  /** @brief Rotate spherical harmonic coefficients. */
  SphericalHarmonics rotate(const SphericalHarmonics &harm) const;

  /** @brief Inverse of spherical harmonic coefficients. */
  SphericalHarmonics inverseRotate(const SphericalHarmonics &harm) const;

  /** @brief Return Euler angles.
  * @code rotary = rotaryZ(gamma)*rotaryX(beta)*rotaryZ(alpha); @endcode */
  void euler(Angle &alpha, Angle &beta, Angle &gamma) const;

  /** @brief Return Euler angles.
  * @code rotary = rotaryZ(yaw)*rotaryY(pitch)*rotaryX(roll); @endcode */
  void cardan(Angle &roll, Angle &pitch, Angle &yaw) const;

  Rotary3d    &operator*=(Rotary3d const &b);
  Rotary3d     operator* (Rotary3d const &b) const;
  Transform3d  operator* (Transform3d const &b) const;

  friend class Transform3d;
  friend Rotary3d rotaryX(Angle angle);
  friend Rotary3d rotaryY(Angle angle);
  friend Rotary3d rotaryZ(Angle angle);
  friend Rotary3d inverse(const Rotary3d &b);
  friend Rotary3d localNorthEastDown(const Vector3d &point);
};

/***********************************************/

/** @brief Rotation about x-axis.
* @param angle in [rad].
* \f[ D=\left(\begin{array}{ccc}
      1 & 0          & 0          \\
      0 & \cos\alpha & \sin\alpha \\
      0 &-\sin\alpha & \cos\alpha
      \end{array}\right) \f] */
Rotary3d rotaryX(Angle angle);

/** @brief Rotation about y-axis.
* @param angle in [rad].
* \f[ D=\left(\begin{array}{ccc}
      \cos\alpha & 0 &-\sin\alpha \\
      0          & 1 & 0          \\
      \sin\alpha & 0 & \cos\alpha
      \end{array}\right) \f] */
Rotary3d rotaryY(Angle angle);

/** @brief Rotation about z-axis.
* @param angle in [rad].
* \f[ D=\left(\begin{array}{ccc}
      \cos\alpha & \sin\alpha & 0 \\
     -\sin\alpha & \cos\alpha & 0 \\
      0          & 0          & 1
      \end{array}\right) \f] */
Rotary3d rotaryZ(Angle angle);

/** @brief Rotary3d of inverse rotation.
* Transposed matrix. */
Rotary3d inverse(const Rotary3d &b);

/** @brief Rotational matrix from local system (north, east, down) to the global system (coordinate system of point). */
//[[deprecated]]
Rotary3d localNorthEastDown(const Vector3d &point);

/** @brief Rotational matrix from local system (north, east, down) to the global system (coordinate system of point).
* Converts @p point to ellipsoidal coordinates beforehand. */
//[[deprecated]]
Rotary3d localNorthEastDown(const Vector3d &point, const Ellipsoid &ellipsoid);

/// @}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

inline Matrix Rotary3d::matrix() const
{
  Matrix R(3,3);
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      R(i,k) = field[i][k];
  return R;
}

/*************************************************/

#endif /* __GROOPS_VECTOR3D__ */


