/***********************************************/
/**
* @file rotary3d.cpp
*
* @brief Rotations in 3d space.
*
* @author Torsten Mayer-Guerr
* @date 2001-05-31
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/ellipsoid.h"
#include "base/vector3d.h"
#include "base/tensor3d.h"
#include "base/transform3d.h"
#include "base/sphericalHarmonics.h"
#include "base/rotary3d.h"

/***********************************************/

// rotary from 3x3 matrix or quaterions
Rotary3d::Rotary3d(const_MatrixSliceRef q)
{
  try
  {
    if((q.rows()==3) && (q.columns()==3)) // (3x3) rotary matrix
    {
      if( q(0,0)*q(1,1)*q(2,2)+q(0,1)*q(1,2)*q(2,0)+q(0,2)*q(1,0)*q(2,1)
         -q(2,0)*q(1,1)*q(0,2)-q(2,1)*q(1,2)*q(0,0)-q(2,2)*q(1,0)*q(0,1) < 0)
        throw(Exception("rotary matrix contains reflection"));

      field[0][0] = q(0,0);
      field[0][1] = q(0,1);
      field[0][2] = q(0,2);
      field[1][0] = q(1,0);
      field[1][1] = q(1,1);
      field[1][2] = q(1,2);
      field[2][0] = q(2,0);
      field[2][1] = q(2,1);
      field[2][2] = q(2,2);
    }
    else if((q.rows()==4) && (q.columns()==1)) // (4x1) quaternion vector
    {
      field[0][0] = q(1,0)*q(1,0)-q(2,0)*q(2,0)-q(3,0)*q(3,0)+q(0,0)*q(0,0);
      field[0][1] = 2*(q(1,0)*q(2,0)-q(3,0)*q(0,0));
      field[0][2] = 2*(q(1,0)*q(3,0)+q(2,0)*q(0,0));
      field[1][0] = 2*(q(1,0)*q(2,0)+q(3,0)*q(0,0));
      field[1][1] = q(2,0)*q(2,0)-q(1,0)*q(1,0)-q(3,0)*q(3,0)+q(0,0)*q(0,0);
      field[1][2] = 2*(q(2,0)*q(3,0)-q(1,0)*q(0,0));
      field[2][0] = 2*(q(1,0)*q(3,0)-q(2,0)*q(0,0));
      field[2][1] = 2*(q(2,0)*q(3,0)+q(1,0)*q(0,0));
      field[2][2] = q(3,0)*q(3,0)-q(2,0)*q(2,0)-q(1,0)*q(1,0)+q(0,0)*q(0,0);
    }
    else
      throw(Exception("matrix is not quaternion nor rotary matrix"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d::Rotary3d(Vector3d x, Vector3d y)
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

Vector Rotary3d::quaternion() const
{
  Vector q(4);
  if( field[0][0] + field[1][1] + field[2][2] > 0.)
  {
    Double t = 1. + field[0][0] + field[1][1] + field[2][2];
    Double s = 0.5/sqrt(t);
    q(0) = s*t;
    q(3) = (field[1][0] - field[0][1]) * s;
    q(2) = (field[0][2] - field[2][0]) * s;
    q(1) = (field[2][1] - field[1][2]) * s;
  }
  else if((field[0][0] > field[1][1]) && (field[0][0] > field[2][2]))
  {
    Double t = 1. + field[0][0] - field[1][1] - field[2][2];
    Double s = 0.5/sqrt(t);
    q(1) = s*t;
    q(2) = (field[1][0] + field[0][1]) * s;
    q(3) = (field[0][2] + field[2][0]) * s;
    q(0) = (field[2][1] - field[1][2]) * s;
  }
  else if(field[1][1] > field[2][2])
  {
    Double t = 1. - field[0][0] + field[1][1] - field[2][2];
    Double s = 0.5/sqrt(t);
    q(2) = s*t;
    q(1) = (field[1][0] + field[0][1]) * s;
    q(0) = (field[0][2] - field[2][0]) * s;
    q(3) = (field[2][1] + field[1][2]) * s;
  }
  else
  {
    Double t = 1. - field[0][0] - field[1][1] + field[2][2];
    Double s = 0.5/sqrt(t);
    q(3) = s*t;
    q(0) = (field[1][0] - field[0][1]) * s;
    q(1) = (field[0][2] + field[2][0]) * s;
    q(2) = (field[2][1] + field[1][2]) * s;
  }
  return q;
}

/***********************************************/

Vector3d Rotary3d::rotate(const Vector3d &v) const
{
  Vector3d erg;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      erg.field[i] += field[i][k] * v.field[k];
  return erg;
}

Vector3d Rotary3d::transform(const Vector3d &v) const {return rotate(v);}

/***********************************************/

Vector3d Rotary3d::inverseRotate(const Vector3d &v) const
{
  Vector3d erg;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      erg.field[i] += field[k][i] * v.field[k];
  return erg;
}

Vector3d Rotary3d::inverseTransform(const Vector3d &v) const {return inverseRotate(v);}

/***********************************************/

Tensor3d Rotary3d::rotate(const Tensor3d &t) const
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

Tensor3d Rotary3d::transform(const Tensor3d &t) const {return rotate(t);}

/***********************************************/

Tensor3d Rotary3d::inverseRotate(const Tensor3d &t) const
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

Tensor3d Rotary3d::inverseTransform(const Tensor3d &t) const {return inverseRotate(t);}

/***********************************************/

SphericalHarmonics Rotary3d::rotate(const SphericalHarmonics &harm) const
{
  return harm.rotate(*this);
}

/***********************************************/

SphericalHarmonics Rotary3d::inverseRotate(const SphericalHarmonics &harm) const
{
  return harm.rotate(inverse(*this));
}

/***********************************************/

Rotary3d &Rotary3d::operator*= (const Rotary3d &b)
{
  *this = *this * b;
  return *this;
}

/***********************************************/

Rotary3d Rotary3d::operator* (const Rotary3d &b) const
{
  Rotary3d erg;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
    {
      erg.field[i][k] = 0.0;
      for(UInt l=0; l<3; l++)
        erg.field[i][k] += field[i][l] * b.field[l][k];
    }
  return erg;
}

/***********************************************/

Transform3d Rotary3d::operator*(const Transform3d &b) const
{
  Transform3d erg;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
    {
      erg.field[i][k] = 0.0;
      for(UInt l=0; l<3; l++)
        erg.field[i][k] += field[i][l] * b.field[l][k];
    }
  return erg;
}

/***********************************************/

void Rotary3d::euler(Angle &alpha, Angle &beta, Angle &gamma) const
{
  beta  = atan2(sqrt(field[2][0]*field[2][0]+field[2][1]*field[2][1]), field[2][2]);
  if(beta == 0.)
  {
    alpha = atan2(field[0][1], field[0][0]);
    gamma = 0;
    return;
  }
  alpha = atan2(field[2][0], -field[2][1]);
  gamma = atan2(field[0][2],  field[1][2]);
}

/***********************************************/

void Rotary3d::cardan(Angle &roll, Angle &pitch, Angle &yaw) const
{
  roll  = -atan2( field[2][1], field[2][2]);
  pitch = -atan2(-field[2][0], sqrt(field[2][1]*field[2][1]+field[2][2]*field[2][2]));
  yaw   = -atan2( field[1][0], field[0][0]);
}

/***********************************************/
/***********************************************/

Rotary3d rotaryX(Angle angle)
{
  const Double c = cos(angle);
  const Double s = sin(angle);
  Rotary3d rot;
  rot.field[1][1] =  c; rot.field[1][2] = s;
  rot.field[2][1] = -s; rot.field[2][2] = c;
  return rot;
}

/***********************************************/

Rotary3d rotaryY(Angle angle)
{
  const Double c = cos(angle);
  const Double s = sin(angle);
  Rotary3d rot;
  rot.field[0][0] = c; rot.field[0][2] = -s;
  rot.field[2][0] = s; rot.field[2][2] =  c;
  return rot;
}

/***********************************************/

Rotary3d rotaryZ(Angle angle)
{
  const Double c = cos(angle);
  const Double s = sin(angle);
  Rotary3d rot;
  rot.field[0][0] =  c; rot.field[0][1] = s;
  rot.field[1][0] = -s; rot.field[1][1] = c;
  return rot;
}

/***********************************************/

Rotary3d inverse(const Rotary3d &b)
{
  Rotary3d rot;
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      rot.field[i][k] = b.field[k][i];
  return rot;
}

/***********************************************/

Rotary3d localNorthEastDown(const Vector3d &point)
{
  const Vector3d z = normalize(-point);                           // down
  const Vector3d y = normalize(crossProduct(z, Vector3d(0,0,1))); // east
  const Vector3d x = normalize(crossProduct(y, z));               // north
  Rotary3d rot;
  rot.field[0][0] = x.x(); rot.field[0][1] = y.x(); rot.field[0][2] = z.x();
  rot.field[1][0] = x.y(); rot.field[1][1] = y.y(); rot.field[1][2] = z.y();
  rot.field[2][0] = x.z(); rot.field[2][1] = y.z(); rot.field[2][2] = z.z();
  return rot;
}

/***********************************************/

Rotary3d localNorthEastDown(const Vector3d &point, const Ellipsoid &ellipsoid)
{
  Angle  L, B;
  Double h;
  ellipsoid(point, L, B, h);
  return localNorthEastDown(polar(L, B, 1.));
}

/***********************************************/

