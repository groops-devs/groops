/***********************************************/
/**
* @file transform3d.cpp
*
* @brief Orthogonal coordinate transformations in 3d space.
* (rotations and reflections).
*
* @author Torsten Mayer-Guerr
* @date 2019-03-03
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/ellipsoid.h"
#include "base/vector3d.h"
#include "base/tensor3d.h"
#include "base/rotary3d.h"
#include "base/transform3d.h"

/***********************************************/

/// Contructor from Rotary3d.
Transform3d::Transform3d(const Rotary3d &rot) : field(rot.field) {}

/***********************************************/

// rotary from 3x3 matrix
Transform3d::Transform3d(const_MatrixSliceRef A)
{
  try
  {
    field[0][0] = A(0,0);
    field[0][1] = A(0,1);
    field[0][2] = A(0,2);
    field[1][0] = A(1,0);
    field[1][1] = A(1,1);
    field[1][2] = A(1,2);
    field[2][0] = A(2,0);
    field[2][1] = A(2,1);
    field[2][2] = A(2,2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Transform3d::Transform3d(Vector3d x, Vector3d y)
{
  // Orthonormalize
  x.normalize();
  Vector3d z = normalize(crossProduct(x,y));
  y = crossProduct(z,x);

  field[0][0] = x.x(); field[0][1] = y.x(); field[0][2] = z.x();
  field[1][0] = x.y(); field[1][1] = y.y(); field[1][2] = z.y();
  field[2][0] = x.z(); field[2][1] = y.z(); field[2][2] = z.z();
}

/***********************************************/

Vector3d Transform3d::transform(const Vector3d &v) const
{
  Vector3d y;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      y.field[i] += field[i][k] * v.field[k];
  return y;
}

/***********************************************/

Vector3d Transform3d::inverseTransform(const Vector3d &v) const
{
  Vector3d y;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      y.field[i] += field[k][i] * v.field[k];
  return y;
}

/***********************************************/

Tensor3d Transform3d::transform(const Tensor3d &t) const
{
  auto rot = [&](UInt i, UInt k)
  {
    return  field[i][0]*field[k][0] * t.xx()
         +  field[i][1]*field[k][1] * t.yy()
         +  field[i][2]*field[k][2] * t.zz()
         + (field[i][0]*field[k][1] + field[i][1]*field[k][0]) * t.xy()
         + (field[i][0]*field[k][2] + field[i][2]*field[k][0]) * t.xz()
         + (field[i][1]*field[k][2] + field[i][2]*field[k][1]) * t.yz();
  };

  Tensor3d t2;
  t2.xx() = rot(0, 0);
  t2.yy() = rot(1, 1);
  t2.zz() = rot(2, 2);
  t2.xy() = rot(0, 1);
  t2.xz() = rot(0, 2);
  t2.yz() = rot(1, 2);
  return t2;
}

/***********************************************/

Tensor3d Transform3d::inverseTransform(const Tensor3d &t) const
{
  auto rot = [&](UInt i, UInt k)
  {
    return  field[0][i]*field[0][k] * t.xx()
         +  field[1][i]*field[1][k] * t.yy()
         +  field[2][i]*field[2][k] * t.zz()
         + (field[0][i]*field[1][k] + field[1][i]*field[0][k]) * t.xy()
         + (field[0][i]*field[2][k] + field[2][i]*field[0][k]) * t.xz()
         + (field[1][i]*field[2][k] + field[2][i]*field[1][k]) * t.yz();
  };

  Tensor3d t2;
  t2.xx() = rot(0, 0);
  t2.yy() = rot(1, 1);
  t2.zz() = rot(2, 2);
  t2.xy() = rot(0, 1);
  t2.xz() = rot(0, 2);
  t2.yz() = rot(1, 2);
  return t2;
}

/***********************************************/

Transform3d &Transform3d::operator*= (const Transform3d &b)
{
  *this = *this * b;
  return *this;
}

/***********************************************/

Transform3d &Transform3d::operator*= (const Rotary3d &b)
{
  *this = *this * b;
  return *this;
}

/***********************************************/

Transform3d Transform3d::operator* (const Transform3d &b) const
{
  Transform3d y;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
    {
      y.field[i][k] = 0.0;
      for(UInt l=0; l<3; l++)
        y.field[i][k] += field[i][l] * b.field[l][k];
    }
  return y;
}

/***********************************************/

Transform3d Transform3d::operator* (const Rotary3d &b) const
{
  Transform3d y;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
    {
      y.field[i][k] = 0.0;
      for(UInt l=0; l<3; l++)
        y.field[i][k] += field[i][l] * b.field[l][k];
    }
  return y;
}

/***********************************************/
/***********************************************/

Transform3d inverse(const Transform3d &b)
{
  Transform3d rot;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      rot.field[i][k] = b.field[k][i];
  return rot;
}

/***********************************************/

Transform3d localNorthEastUp(const Vector3d &point)
{
  const Vector3d z = normalize(point);                            // up
  const Vector3d y = normalize(crossProduct(Vector3d(0,0,1), z)); // east
  const Vector3d x = normalize(crossProduct(z, y));               // north
  Transform3d T;
  T.field[0][0] = x.x(); T.field[0][1] = y.x(); T.field[0][2] = z.x();
  T.field[1][0] = x.y(); T.field[1][1] = y.y(); T.field[1][2] = z.y();
  T.field[2][0] = x.z(); T.field[2][1] = y.z(); T.field[2][2] = z.z();
  return T;
}

/***********************************************/

Transform3d localNorthEastUp(const Vector3d &point, const Ellipsoid &ellipsoid)
{
  Angle  L, B;
  Double h;
  ellipsoid(point, L, B, h);
  return localNorthEastUp(polar(L, B, 1.));
}

/***********************************************/

