/***********************************************/
/**
* @file vector3d.h
*
* @brief Vector in 3d space.
*
* @author Torsten Mayer-Guerr
* @date 2001-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_VECTOR3D__
#define __GROOPS_VECTOR3D__

#include "base/importStd.h"
#include "base/angle.h"
#include "base/matrix.h"

/**
* @defgroup vector3dGroup Vector3d
* @brief Coordinates in 3d space.
* @ingroup base */
/// @{

/***** CLASS ***********************************/

/** @brief Vector in 3d space.
* For the representation of e.g. positions, velocities, gravity.
* Can be rotated with @a Rotary3d.
* Can added, subtracted and multiplied with a Double.
* (Internally represented by cartesian coordinates). */
class Vector3d
{
  std::array<Double,3> field;

public:
  Vector3d() : field{0.,0.,0.} {} //!< Default constructor (Zero vector).

  Vector3d(Double x, Double y, Double z) : field{x,y,z} {} //!< Constructor from cartesian coordinates.
  explicit Vector3d(const_MatrixSliceRef x);               //!< Constructor from 3x1 or 1x3 matrix slice

  Double  x() const {return field[0];} //!< Cartesian coordinates: x component.
  Double  y() const {return field[1];} //!< Cartesian coordinates: y component.
  Double  z() const {return field[2];} //!< Cartesian coordinates: z component.
  Double &x()       {return field[0];} //!< Cartesian coordinates: x component.
  Double &y()       {return field[1];} //!< Cartesian coordinates: y component.
  Double &z()       {return field[2];} //!< Cartesian coordinates: z component.

  /** @brief Polar coordinates: longitude (-PI,PI].
  * @f[ \lambda = \arctan2(y,x) @f] */
  Angle lambda()  const;

  /** @brief Polar coordinates: latitude [-PI,PI].
  * @f[ \varphi = \arctan2(z,\sqrt{x^2+y^2}) @f] */
  Angle phi() const;

  /** @brief Polar coordinates: zenit angle [-PI,PI].
  * @f[ \vartheta = \pi/2 - \varphi @f]
  * @see Vector3d::phi */
  Angle theta() const;

  /** @brief Polar coordinates: radius.
  * @f[ r = \sqrt{x^2+y^2+z^2} @f] */
  Double r() const;

  /** @brief Quadratic sum.
  * @f[ x^2+y^2+z^2 @f] */
  Double quadsum() const;

  /** @brief L2-Norm.
  * @f[ \sqrt{x^2+y^2+z^2} @f]*/
  Double norm() const;

  /** @brief Normalize the vector.
  * @f[ \frac{1}{\sqrt{x^2+y^2+z^2}} \cdot (x,y,z)^T @f]
  * @return old length */
  Double normalize();

  /// Cast to Vector.
  Vector vector() const;

  Vector3d &operator+= (Vector3d const &b);
  Vector3d &operator-= (Vector3d const &b);
  Vector3d &operator*= (Double  c);
  Vector3d &operator/= (Double  c);

  friend class Rotary3d;
  friend class Transform3d;
  friend Double inner(const Vector3d &x, const Vector3d &y);
  friend Vector3d polar(Angle lambda, Angle phi, Double r);
  friend Vector3d crossProduct(const Vector3d &x, const Vector3d &y);
};

/***** FUNCTIONS *******************************/

/** @brief Inner Product.
* @f[ c = x^Ty @f] */
inline Double inner(const Vector3d &x, const Vector3d &y);

/** @brief Vector3d from polar coordinates. */
inline Vector3d polar(Angle lambda, Angle phi, Double r);

/** @brief Cross product.
* @f[ z = x \times y @f] */
inline Vector3d crossProduct(const Vector3d &x, const Vector3d &y);

/** @brief Scale vector to unit length. */
inline Vector3d normalize(const Vector3d &x);

inline const Vector3d operator- (const Vector3d &t)                      {return Vector3d(t)  *= -1;}
inline const Vector3d operator+ (const Vector3d &t1, const Vector3d &t2) {return Vector3d(t1) += t2;}
inline const Vector3d operator- (const Vector3d &t1, const Vector3d &t2) {return Vector3d(t1) -= t2;}
inline const Vector3d operator* (Double c, const Vector3d &t)            {return Vector3d(t)  *=c;}
inline const Vector3d operator* (const Vector3d &t, Double c)            {return Vector3d(t)  *=c;}
inline const Vector3d operator/ (const Vector3d &t, Double c)            {return Vector3d(t)  /=c;}

/// @}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

inline Vector3d::Vector3d(const_MatrixSliceRef x) {if(x.size()!=3) throw(Exception("Vector3d constructor dimension error")); for(UInt i=0; i<3; i++) field[i]=(x.rows()>1 ? x(i,0) : x(0,i));}

inline Vector3d &Vector3d::operator+= (Vector3d const &b)  {for(UInt i=0; i<3; i++) field[i]+=b.field[i]; return *this;}
inline Vector3d &Vector3d::operator-= (Vector3d const &b)  {for(UInt i=0; i<3; i++) field[i]-=b.field[i]; return *this;}
inline Vector3d &Vector3d::operator*= (Double  c)          {for(UInt i=0; i<3; i++) field[i]*=c; return *this;}
inline Vector3d &Vector3d::operator/= (Double  c)          {for(UInt i=0; i<3; i++) field[i]/=c; return *this;}

inline Double Vector3d::quadsum()const {return field[0]*field[0]+field[1]*field[1]+field[2]*field[2];}
inline Angle  Vector3d::lambda() const {return Angle(atan2(field[1],field[0]));}
inline Angle  Vector3d::phi()    const {return Angle(atan2(field[2],sqrt(field[0]*field[0]+field[1]*field[1])));}
inline Angle  Vector3d::theta()  const {return Angle(PI/2-atan2(field[2],sqrt(field[0]*field[0]+field[1]*field[1])));}
inline Double Vector3d::r()      const {return norm();}
inline Double Vector3d::norm()   const {return sqrt(quadsum());}
inline Double Vector3d::normalize()    {Double n=norm(); *this*=(1/n); return n;}
inline Vector Vector3d::vector() const {Vector r(3); r(0)=x(); r(1)=y(); r(2)=z(); return r;}

inline Double   inner(const Vector3d &a, const Vector3d &b)        {return a.field[0]*b.field[0]+a.field[1]*b.field[1]+a.field[2]*b.field[2];}
inline Vector3d crossProduct(const Vector3d &a, const Vector3d &b) {return Vector3d(a.field[1]*b.field[2]-a.field[2]*b.field[1],
                                                                                    a.field[2]*b.field[0]-a.field[0]*b.field[2],
                                                                                    a.field[0]*b.field[1]-a.field[1]*b.field[0]);}
inline Vector3d polar(Angle lambda, Angle phi, Double r) {return Vector3d(r*cos(lambda)*cos(phi), r*sin(lambda)*cos(phi), r*sin(phi));}
inline Vector3d normalize(const Vector3d &x) {return x/x.norm();}


/***********************************************/

#endif /* __GROOPS_VECTOR3D__ */
