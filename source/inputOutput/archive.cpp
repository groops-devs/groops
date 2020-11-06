/***********************************************/
/**
* @file archive.cpp
*
* @brief Interface for streaming of data in different formats.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-01
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/matrix.h"
#include "base/vector3d.h"
#include "base/tensor3d.h"
#include "base/rotary3d.h"
#include "base/transform3d.h"
#include "base/ellipsoid.h"
#include "archive.h"

/***********************************************/

UInt InArchive::versionStr2version(const std::string &str)
{
  if(str.empty())
    return 0;
  if(str == "1.1")
    return 1;
  if(str == "1.2")
    return 2;
  std::stringstream ss(str);
  UInt v;
  ss>>v;
  return v;
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const Vector3d &x)
{
  ar<<nameValue("x", x.x())<<nameValue("y", x.y())<<nameValue("z", x.z());
}

/***********************************************/

template<> void load(InArchive  &ar, Vector3d &x)
{
  ar>>nameValue("x", x.x())>>nameValue("y", x.y())>>nameValue("z", x.z());
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const Tensor3d &x)
{
  ar<<nameValue("xx", x.xx())<<nameValue("yy", x.yy())<<nameValue("zz", x.zz());
  ar<<nameValue("xy", x.xy())<<nameValue("xz", x.xz())<<nameValue("yz", x.yz());
}

/***********************************************/

template<> void load(InArchive  &ar, Tensor3d &x)
{
  if(ar.version() < 20170920)
  {
    ar>>nameValue("xx", x.xx())>>nameValue("xy", x.xy())>>nameValue("xz", x.xz());
    ar>>nameValue("yy", x.yy())>>nameValue("yz", x.yz());
    ar>>nameValue("zz", x.zz());
    return;
  }

  ar>>nameValue("xx", x.xx())>>nameValue("yy", x.yy())>>nameValue("zz", x.zz());
  ar>>nameValue("xy", x.xy())>>nameValue("xz", x.xz())>>nameValue("yz", x.yz());
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const Rotary3d &x)
{
  Vector q = x.quaternion();
  ar<<nameValue("q0", q(0));
  ar<<nameValue("qx", q(1));
  ar<<nameValue("qy", q(2));
  ar<<nameValue("qz", q(3));
}

/***********************************************/

template<> void load(InArchive  &ar, Rotary3d &x)
{
  if(ar.version() < 20170920)
  {
    Matrix R(3,3);
    ar>>nameValue("xx", R(0,0))>>nameValue("xy", R(0,1))>>nameValue("xz", R(0,2));
    ar>>nameValue("yx", R(1,0))>>nameValue("yy", R(1,1))>>nameValue("yz", R(1,2));
    ar>>nameValue("zx", R(2,0))>>nameValue("zy", R(2,1))>>nameValue("zz", R(2,2));
    x = Rotary3d(R);
    return;
  }

  Vector q(4);
  ar>>nameValue("q0", q(0));
  ar>>nameValue("qx", q(1));
  ar>>nameValue("qy", q(2));
  ar>>nameValue("qz", q(3));
  x = Rotary3d(q);
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const Transform3d &x)
{
  Matrix R = x.matrix();
  ar<<nameValue("xx", R(0,0))<<nameValue("xy", R(0,1))<<nameValue("xz", R(0,2));
  ar<<nameValue("yx", R(1,0))<<nameValue("yy", R(1,1))<<nameValue("yz", R(1,2));
  ar<<nameValue("zx", R(2,0))<<nameValue("zy", R(2,1))<<nameValue("zz", R(2,2));
}

/***********************************************/

template<> void load(InArchive  &ar, Transform3d &x)
{
  std::array<std::array<Double,3>,3> T;
  ar>>nameValue("xx", T[0][0])>>nameValue("xy", T[0][1])>>nameValue("xz", T[0][2]);
  ar>>nameValue("yx", T[1][0])>>nameValue("yy", T[1][1])>>nameValue("yz", T[1][2]);
  ar>>nameValue("zx", T[2][0])>>nameValue("zy", T[2][1])>>nameValue("zz", T[2][2]);
  x = Transform3d(T);
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const Ellipsoid &x)
{
  ar<<nameValue("a", x.a())<<nameValue("b", x.b());
}

/***********************************************/

template<> void load(InArchive  &ar, Ellipsoid &x)
{
  Double a,b;
  ar>>nameValue("a", a)>>nameValue("b", b);
  x = Ellipsoid(a, ((a-b)!=0.) ? a/(a-b) : 0.0);
}

/***********************************************/
