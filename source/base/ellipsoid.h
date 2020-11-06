/***********************************************/
/**
* @file ellipsoid.h
*
* @brief Transformation of ellipsoidial coordinates.
*
* @author Torsten Mayer-Guerr
* @date 2004-10-25
*
*/
/***********************************************/

#ifndef __GROOPS_ELLIPSOID__
#define __GROOPS_ELLIPSOID__

#include "base/importStd.h"
#include "base/angle.h"
#include "base/vector3d.h"

/***** CLASS ***********************************/

/** @brief Transformation of ellipsoidial coordinates.
* @ingroup vector3d */
class Ellipsoid
{
  Double _a,_b;

public:
  /** @brief Ellipsoid from semi major axis and inverse flatenning.
  * if f=0, a sphere is assumed. */
  Ellipsoid(Double a=DEFAULT_GRS80_a, Double f=DEFAULT_GRS80_f) : _a(a), _b((f != 0.) ? (a*(1-1/f)) : a) {}

  /** @brief Computes ellipsoidal coordinates from @a point.
  * @param point point
  * @param[out] L longitude (-PI,PI]
  * @param[out] B latitude [-PI,PI]
  * @param[out] h height [m]
  */
  void operator()(const Vector3d &point, Angle &L, Angle &B, Double &h) const;

  /** @brief Transformation ellipsoidal coordinates in @a Vector3d.
  * @param L longitude (-PI,PI]
  * @param B latitude [-PI,PI]
  * @param h height [m]
  */
  const Vector3d operator()(Angle L, Angle B, Double h) const;

  /// Semi major axis.
  Double a() const {return _a;}

  /// Semi minor axis.
  Double b() const {return _b;}

  /// Numerical  excentricity.
  Double e() const {return sqrt(_a*_a-_b*_b)/_a;}

  /// Flattening.
  Double f() const {return (_a-_b)/_a;}
};

/*************************************************/

#endif
