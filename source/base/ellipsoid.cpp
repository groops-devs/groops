/***********************************************/
/**
* @file ellipsoid.cpp
*
* @brief Transformation of ellipsoidial coordinates.
*
* @author Torsten Mayer-Guerr
* @date 2004-10-25
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/ellipsoid.h"

/***********************************************/

void Ellipsoid::operator()(const Vector3d &point, Angle &L, Angle &B, Double &h) const
{
  if(point.quadsum() == 0.)
    throw(Exception("In Ellipsoid:\n0 Vector"));

  const Double e2  = sqrt((_a/_b-1)*(_a/_b+1));
  const Double z   = point.z();
  const Double rho = std::sqrt(point.x()*point.x()+point.y()*point.y());

  // start values
  L = (rho != 0.) ? point.lambda() : Angle(0);
  B = Angle(std::atan2(z*(1+e2*e2), rho));
  h = 0.0;

  // Iteration
  Double h1, B1;
  UInt   count = 0;
  do
  {
    h1 = h;
    B1 = B;
    const Double N = _a*(_a/(_b*std::sqrt(1+std::pow(e2*std::cos(B), 2))));
    h = (std::fabs(B) < (60*DEG2RAD)) ? (rho/std::cos(B) - N) : (z/std::sin(B) - N/(1+e2*e2));
    B = Angle(std::atan2(z*(1+e2*e2), rho*(1+e2*e2*h/(N+h))));
  }
  while(((std::fabs(h-h1) > 1e-6) || (std::fabs(B-B1) > 1e-10)) && (count++ < 1000));
}

/***********************************************/

const Vector3d Ellipsoid::operator()(Angle L, Angle B, Double h) const
{
  Double e2 = std::sqrt((_a/_b-1)*(_a/_b+1));
  Double N  = _a*(_a/(_b*std::sqrt(1+std::pow(e2*std::cos(B), 2))));

  return Vector3d((N+h) * std::cos(B) * std::cos(L),
                  (N+h) * std::cos(B) * std::sin(L),
                  (N/(1+e2*e2)+h)* std::sin(B));
}

/***********************************************/
