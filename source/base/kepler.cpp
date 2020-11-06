/***********************************************/
/**
* @file kepler.cpp
*
* @brief Keplerian elements and orbits.
*
* @author Torsten Mayer-Guerr
* @date 2004-09-12
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/vector3d.h"
#include "base/time.h"
#include "base/kepler.h"
#include "base/equinoctial.h"

/***********************************************/

Kepler::Kepler() :
  GM(DEFAULT_GM),
  Omega(0),
  i(0),
  omega(0),
  a(0),
  e(0),
  tau()
{
}

/***********************************************/

Kepler::Kepler(Double _Omega, Double _i, Double _omega, Double _a, Double _e, const Time &_tau, Double _GM) :
  GM(_GM),
  Omega(_Omega),
  i( _i),
  omega(_omega),
  a(_a),
  e(_e),
  tau(_tau)
{
}

/***********************************************/

Kepler::Kepler(const Time &time, Double _Omega, Double _i, Double _omega, Double _a, Double _e, Double _M, Double _GM) :
  GM(_GM),
  Omega(_Omega),
  i( _i),
  omega(_omega),
  a(_a),
  e(_e),
  tau(time - seconds2time(_M*sqrt(a*a*a/GM)))
{
}

/***********************************************/

Kepler::Kepler(const Equinoctial &equinoctial) :
  GM(equinoctial.GM),
  a(equinoctial.a),
  tau(equinoctial.tau)
{
  e = sqrt(equinoctial.h*equinoctial.h+equinoctial.k*equinoctial.k);
  i = 2*atan(sqrt(equinoctial.p*equinoctial.p+equinoctial.q*equinoctial.q));

  const Double sinOmega  = equinoctial.p/sqrt(equinoctial.p*equinoctial.p+equinoctial.q*equinoctial.q);
  const Double cosOmega  = equinoctial.q/sqrt(equinoctial.p*equinoctial.p+equinoctial.q*equinoctial.q);
  Omega = atan2(sinOmega,cosOmega);

  const Double sinXi = equinoctial.h / sqrt(equinoctial.h*equinoctial.h+equinoctial.k*equinoctial.k);
  const Double cosXi = equinoctial.k / sqrt(equinoctial.h*equinoctial.h+equinoctial.k*equinoctial.k);
  const Double xi    = atan2(sinXi,cosXi);
  const Double M     = static_cast<Double>(equinoctial.l0) - xi;

  omega = xi - Omega; //equinoctial.l0 - Omega;
  tau   = tau - seconds2time(M*sqrt(a*a*a/GM));

  if(Omega<0) Omega += 2*PI;
  if(omega<0) omega += 2*PI;
  if(i<0)     i     += 2*PI;
}

/***********************************************/

Kepler::Kepler(const Time &time, const Vector3d &position, const Vector3d &velocity, Double _GM) : GM(_GM)
{
  const Double   r = position.r();
  const Vector3d C = crossProduct(position, velocity);

  i     = atan2(sqrt(C.x()*C.x()+C.y()*C.y()), C.z());
  Omega = atan2(C.x(), -C.y());

  // major semi axis
  a = 1/(2/position.r()-inner(velocity,velocity)/GM);

  // excentricity
  const Double p = inner(C,C)/GM;
  e  = sqrt((a-p)/a);

  // true anomaly
  const Double v = atan2(p*inner(position,velocity)/C.r(), p-r);

  // oribt system
  const Vector3d P = (e+cos(v))/p * position - r*sin(v)/C.r() * velocity;
  const Vector3d Q = sin(v)/p     * position + r*cos(v)/C.r() * velocity;
  const Vector3d K = Vector3d(cos(Omega), sin(Omega), 0);

  // argument of perigee
  omega = atan2(-inner(K,Q), inner(K,P));

  // excentric anomaly
  Double E = atan2(inner(position,Q)/(a*sqrt(1-e*e)), inner(position,P)/a+e);

  // positive angles
  if(i<0)     i     += 2*PI;
  if(Omega<0) Omega += 2*PI;
  if(omega<0) omega += 2*PI;
  if(E<0)     E     += 2*PI;

  // time of perigee
  const Double M = E - e*sin(E);
  tau = time - seconds2time(M*sqrt(a*a*a/GM));
}

/***********************************************/

Double Kepler::meanAnomaly(const Time &time) const
{
  Double M = (time-tau).seconds()*sqrt(GM/(a*a*a));
  if(M<0) M += 2*PI;
  return M;
}

/***********************************************/

Double Kepler::trueAnomaly(const Time &time) const
{
  Double M = meanAnomaly(time);

  // excentric anomaly
  Double E = M;
  Double E_old;
  Double eps = 1e-12;

  UInt iter = 0;
  do
  {
    E_old = E;
    E = M + e*sin(E);
  }
  while((fabs(E-E_old)>eps)&&(iter++<100));

  return atan2(sqrt(1-e*e)*sin(E), cos(E)-e);
}

/***********************************************/

const Vector3d Kepler::position(const Time &time) const
{
  Vector3d position, velocity;
  orbit(time, position, velocity);
  return position;
}

/***********************************************/

void Kepler::orbit(const Time &time, Vector3d &position, Vector3d &velocity) const
{
  Vector3d acceleration;
  orbit(time, position, velocity, acceleration);
}

/***********************************************/

void Kepler::orbit(const Time &time, Vector3d &position, Vector3d &velocity, Vector3d &acceleration) const
{
  Double v = trueAnomaly(time);

  // distance and angular momentum
  Double p = a*(1-e*e);
  Double r = p/(1+e*cos(v));
  Double C = sqrt(GM*p);

  // orbital system
  Vector3d P(cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(i),
             sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(i),
             sin(omega)*sin(i));

  Vector3d Q(-cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(i),
             -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(i),
             cos(omega)*sin(i));


  // position and velocity
  position     = (r*cos(v))    * P + (r*sin(v))       * Q;
  velocity     = (-C/p*sin(v)) * P + (C*(cos(v)+e)/p) * Q;
  acceleration = -GM/pow(position.r(),3) * position;
}

/***********************************************/
/***********************************************/
