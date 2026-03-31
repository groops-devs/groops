/***********************************************/
/**
* @file transform3d.h
*
* @brief Orthogonal coordinate transformations in 3d space.
* (rotations and reflections).
*
* @author Torsten Mayer-Guerr
* @date 2019-03-03
*
*/
/***********************************************/

#ifndef __GROOPS_TRANSFORM3D__
#define __GROOPS_TRANSFORM3D__

#include "base/importStd.h"
#include "base/matrix.h"

/** @addtogroup vector3dGroup */
/// @{

class Angle;
class Vector3d;
class Tensor3d;
class Rotary3d;
class Ellipsoid;

/***** CLASS ***********************************/

/** @brief Coordinate transformations in 3d space.
* (rotations and reflections). */
class Transform3d
{
  std::array<std::array<Double,3>,3> field;

public:
  Transform3d() : field{{{1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}} {} //!< Default constructor (Unitary matrix).

  /// Contructor from std::array.
  Transform3d(const std::array<std::array<Double,3>,3> &x) : field(x) {}

  /// Contructor from Rotary3d.
  Transform3d(const Rotary3d &rot);

  /// Contructor from matrix (3x3).
  explicit Transform3d(const_MatrixSliceRef A);

  /** @brief Rotation of a local system.
  * Rotational matrix from local system (defined by the given axis x,y,z)
  * in the global system. The axis will be normalized and orthogonalized before.
  * The z-axis is defined orthogonal to x and y.
  * Zhe y-axis is defined orthogonal to the new z-axis and x. */
  Transform3d(Vector3d x, Vector3d y);

  /// Cast to Matrix.
  Matrix matrix() const;

  /** @brief Transformation of a vector.
  * \f[ y = T \cdot v \f]. */
  Vector3d transform(const Vector3d &v) const;

  /** @brief Inverse transformation of a vector.
  * \f[ y = T^T \cdot v \f]. */
  Vector3d inverseTransform(const Vector3d &v) const;

  /** @brief Transformation of  a tensor.
  * Both sides of the dyadic product are rotated.
  * \f[ y = T \cdot t \cdot T^T \f]. */
  Tensor3d transform(const Tensor3d &t) const;

  /** @brief Inverse transformation of a tensor.
  * Both sides of the dyadic product are rotated.
  * \f[ y = T^T \cdot t \cdot T \f]. */
  Tensor3d inverseTransform(const Tensor3d &t) const;

  Transform3d &operator*=(const Transform3d &b);
  Transform3d &operator*=(const Rotary3d &b);
  Transform3d  operator* (const Transform3d &b) const;
  Transform3d  operator* (const Rotary3d &b) const;

  friend class Rotary3d;
  friend Transform3d inverse(const Transform3d &b);
  friend Transform3d localNorthEastUp(const Vector3d &point);
  friend Transform3d localNorthEastUp(const Vector3d &point, const Ellipsoid &ellipsoid);
};

/***********************************************/

/** @brief flip x-axis. */
inline Transform3d flipX() {return Transform3d(std::array<std::array<Double,3>,3>{{{-1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}});}

/** @brief flip y-axis. */
inline Transform3d flipY() {return Transform3d(std::array<std::array<Double,3>,3>{{{1.,0.,0.}, {0.,-1.,0.}, {0.,0.,1.}}});}

/** @brief flip z-axis. */
inline Transform3d flipZ() {return Transform3d(std::array<std::array<Double,3>,3>{{{1.,0.,0.}, {0.,1.,0.}, {0.,0.,-1.}}});}

/** @brief Transform3d of inverse rotation.
* @ingroup vector3dGroup
* Transposed matrix. */
Transform3d inverse(const Transform3d &b);

/** @brief Rotational matrix from local system (north, east, up) to the global system (coordinate system of point).
* @ingroup vector3dGroup */
Transform3d localNorthEastUp(const Vector3d &point);

/** @brief Rotational matrix from local system (north, east, up) to the global system (coordinate system of point).
* @ingroup vector3dGroup
* Converts @p point to ellipsoidal coordinates beforehand. */
Transform3d localNorthEastUp(const Vector3d &point, const Ellipsoid &ellipsoid);

/// @}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

inline Matrix Transform3d::matrix() const
{
  Matrix R(3, 3, Matrix::NOFILL);
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      R(i,k) = field[i][k];
  return R;
}

/*************************************************/

#endif /* __GROOPS_TRANSFORM3D__ */


