/***********************************************/
/**
* @file equinoctial.cpp
*
* @brief Equinoctial elements and orbits.
*
* @author Matthias Ellmer
* @date 2015-06-10
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/vector3d.h"
#include "base/time.h"
#include "base/kepler.h"
#include "base/equinoctial.h"

/***********************************************/

Equinoctial::Equinoctial() :
  GM(DEFAULT_GM),
  a(0),
  h(0),
  k(0),
  p(0),
  q(0),
  l0(0),
  tau()
{
}

/***********************************************/

Equinoctial::Equinoctial(Double a, Double h, Double k, Double p, Double q, Double l0, const Time &tau, Double GM) :
  GM(GM),
  a(a),
  h(h),
  k(k),
  p(p),
  q(q),
  l0(l0),
  tau(tau)
{
}

/***********************************************/

Equinoctial::Equinoctial(const Kepler &kepler) :
  GM(kepler.GM),
  a(kepler.a),
  h(kepler.e*std::sin(kepler.omega + kepler.Omega)),
  k(kepler.e*std::cos(kepler.omega + kepler.Omega)),
  p(tan(kepler.i / 2) *std::sin(kepler.Omega)),
  q(tan(kepler.i / 2) *std::cos(kepler.Omega)),
  l0(kepler.meanAnomaly(kepler.tau) + kepler.omega + kepler.Omega),
  tau(kepler.tau)
{
}

/***********************************************/

Equinoctial::Equinoctial(const Time &epoch, const Vector3d &position, const Vector3d &velocity, Double _GM) :
  Equinoctial(epoch, epoch, position, velocity, _GM)
{
}

/***********************************************/

Equinoctial::Equinoctial(const Time &epoch, const Time &time, const Vector3d &position, const Vector3d &velocity, Double _GM) :
  GM(_GM),
  tau(epoch)
{
  try
  {
    // Semi-major axis
    a = 1 / (2 / position.norm() - velocity.norm() * velocity.norm() / GM);

    // First orthonormal basis vector in equinoctial reference frame
    Vector3d f, g, w;
    w = crossProduct(position, velocity) * (1 / crossProduct(position, velocity).norm());

    // Elements p and q
    p =  w.x() / (1 + w.z());
    q = -w.y() / (1 + w.z());

    // Remaining basis vectors using p and q
    orthonormals(f, g);

    // Ecccentricity vector
    Vector3d e_ = -position * (1 / position.norm()) + crossProduct(velocity, crossProduct(position, velocity)) * (1 / GM);

    // Elements h and k
    h = inner(e_, g);
    k = inner(e_, f);

    // Position in equinoctial frame
    Double X = inner(position, f);
    Double Y = inner(position, g);

    // Eccentric longitude
    Double b = 1 / (1 + std::sqrt(1 - h * h - k * k));
    Double sinF = h + ((1 - h * h * b) * Y - h * k * b * X) / (a * std::sqrt(1 - h * h - k * k));
    Double cosF = k + ((1 - k * k * b) * X - h * k * b * Y) / (a * std::sqrt(1 - h * h - k * k));
    Double    F = atan2(sinF, cosF);

    // Mean longitude
    LongDouble lambda = F + h * std::cos(F) - k * std::sin(F);

    // The easy method with no LongDoubles would be
    //   l0 = lambda - n * (time - epoch).seconds();
    // See Equinoctial::meanLongitude() for reasoning.

    // Refrence mean longitude at epoch
    LongDouble thisGM = GM;
    LongDouble a3 = std::pow(a,3);
    LongDouble n = std::sqrt(thisGM / a3);
    LongDouble l0_Int = (time.mjdInt()-epoch.mjdInt()) * (24*60*60);
    LongDouble l0_Mod = (time.mjdMod()-epoch.mjdMod()) * (24*60*60);

    l0 = lambda - n * l0_Int - n * l0_Mod;

    // correct domain
    const LongDouble lpi = 3.14159265358979323846264338328L;
    l0 = std::fmod(l0, 2*lpi);
    if (l0 < 0)
      l0 += 2 * lpi;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

LongDouble Equinoctial::meanLongitude(const Time &time) const
{
  try
  {
    // See comment in constructor. Again, this would be the easy method:
    //     LongDouble lambda = l0 + (time - tau).seconds() * std::sqrt(GM / (a * a * a));
    //     Double lambda = l0 + ((time - tau).seconds()/a) * std::sqrt(GM / a );
    LongDouble thisGM = GM;
    LongDouble a3 = std::pow(a,3);
    LongDouble n = std::sqrt(thisGM / a3);
    LongDouble l0_Int = (time.mjdInt()-tau.mjdInt()) * (24*60*60);
    LongDouble l0_Mod = (time.mjdMod()-tau.mjdMod()) * (24*60*60);

    LongDouble lambda = l0 + n * l0_Int + n * l0_Mod;

    const LongDouble lpi = 3.14159265358979323846264338328L;
    lambda = std::fmod(lambda, 2*lpi);
    if (lambda < 0)
      lambda += 2 * lpi;

    return lambda;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Equinoctial::eccentricLongitude(const Time &time) const
{
  try
  {
    LongDouble lambda = meanLongitude(time);
    LongDouble F = lambda;
    LongDouble F_old;
    LongDouble eps = 1e-18;

    // Newton-Raphson method
    UInt iter = 0;
    do
    {
      F_old = F;
      F = F - (F - k * std::sin(F) + h * std::cos(F) - lambda) / (1 - k * std::cos(F) - h * std::sin(F));
    }
    while ((std::fabs(F - F_old) > eps) && (iter++ < 100));


    return static_cast<Double>(F);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Equinoctial::trueLongitude(const Time &time) const
{
  try
  {
    Double beta = 1 / (1 + std::sqrt(1 - h * h - k * k)); // Supplementary quantity
    Double F    = eccentricLongitude(time);
    Double sinL = ((1 - k * k * beta) * std::sin(F) + h * k * beta * std::cos(F) - h) / (1 - h * std::sin(F) - k * std::cos(F));
    Double cosL = ((1 - h * h * beta) * std::cos(F) + h * k * beta * std::sin(F) - k) / (1 - h * std::sin(F) - k * std::cos(F));

    return atan2(sinL, cosL);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Equinoctial::meanAnomaly(const Time &time) const
{
  try
  {
    Double sinXi = h / std::sqrt(h * h + k * k);
    Double cosXi = k / std::sqrt(h * h + k * k);
    return static_cast<Double>(meanLongitude(time)) - std::atan2(sinXi, cosXi);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Equinoctial::orthonormals(Vector3d &f, Vector3d &g) const
{
  try
  {
    Vector3d w;
    orthonormals(f, g, w);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

void Equinoctial::orthonormals(Vector3d &f, Vector3d &g, Vector3d &w) const
{
  try
  {
    Double alpha = 1 / (1 + p * p + q * q);
    f = alpha * Vector3d(1 - p * p + q * q, 2 * p * q, -2 * p);
    g = alpha * Vector3d(2 * p * q, 1 + p * p - q * q, 2 * q);
    w = alpha * Vector3d(2 * p, -2 * q, 1 - p * p - q * q);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Equinoctial::stateTransitionMatrix(const Time &time)
{
  try
  {
    // Prerequisites

    Vector3d f, g, w;
    orthonormals(f, g, w);

    Vector3d inPlanePosition, inPlaneVelocity;
    inPlaneOrbit(time, inPlanePosition, inPlaneVelocity);
    Double X  = inPlanePosition.x();
    Double Y  = inPlanePosition.y();
    Double XD = inPlaneVelocity.x();
    Double YD = inPlaneVelocity.y();

    Vector3d position, velocity;
    orbit(time, position, velocity);

    // Supplementary quantities
    Double A = std::sqrt(GM * a);
    Double B = std::sqrt(1 - h * h - k * k);
    Double C = 1 + p * p + q * q;

    Double F =  eccentricLongitude(time);

    Double n = std::sqrt(GM / (a * a * a));
    Double r = a * (1 - h * std::sin(F) - k * std::cos(F)); //inPlanePosition.norm();

    // Partial derivatives
    Double dXdh = -k * XD / (n * (1. + B)) + a * Y * YD / (A * B);
    Double dYdh = -k * YD / (n * (1. + B)) - a * X * YD / (A * B) - a;
    Double dXdk =  h * XD / (n * (1. + B)) + a * Y * XD / (A * B) - a;
    Double dYdk =  h * YD / (n * (1. + B)) - a * X * XD / (A * B);

    Double dXDdh =  a * YD * YD / (A * B) + (A / (r * r * r)) * (a * k * X / (1 + B) - Y * Y / B);
    Double dYDdh = -a * XD * YD / (A * B) + (A / (r * r * r)) * (a * k * Y / (1 + B) + X * Y / B);
    Double dXDdk =  a * XD * YD / (A * B) - (A / (r * r * r)) * (a * h * X / (1 + B) + X * Y / B);
    Double dYDdk = -a * XD * XD / (A * B) - (A / (r * r * r)) * (a * h * Y / (1 + B) - X * X / B);

    // State transition matrix
    // Chapter 2.1.6 -- Danielson, Donald A., Christopher Patrick Sagovac, Beny Neta, and Leo W. Early. 1995. “Semianalytic Satellite Theory.” DTIC Document.
    // Broucke, R. A., and P. J. Cefola. 1972. “On the Equinoctial Orbit Elements.” Celestial Mechanics 5 (3): 303–10. doi:10.1007/BF01228432.
    // TAKE CARE as the derivatives dr/da and drd/da are wrong in Danielson. Check Broucke for these.

    // Order is dr/da dr/dh dr/dk dr/dp dr/dq dr/dlambda
    // Position derivatives
    Matrix Phi(6, 6);

    {
      // dr/da
      Vector3d tmp = (1./a) * (position - velocity * 1.5 * (time -tau).seconds());
      Phi(0, 0) = tmp.x();
      Phi(1, 0) = tmp.y();
      Phi(2, 0) = tmp.z();
    }
    {
      // dr/dh
      Vector3d tmp = dXdh * f + dYdh * g;
      Phi(0, 1) = tmp.x();
      Phi(1, 1) = tmp.y();
      Phi(2, 1) = tmp.z();
    }
    {
      // dr/dk
      Vector3d tmp = dXdk * f + dYdk * g;
      Phi(0, 2) = tmp.x();
      Phi(1, 2) = tmp.y();
      Phi(2, 2) = tmp.z();
    }
    {
      // dr/dp
      Vector3d tmp = 2.* (q * (Y * f - X * g) - X * w) / C;
      Phi(0, 3) = tmp.x();
      Phi(1, 3) = tmp.y();
      Phi(2, 3) = tmp.z();
    }
    {
      // dr/dq
      Vector3d tmp = 2.* (p * (X * g - Y * f) + Y * w) / C;
      Phi(0, 4) = tmp.x();
      Phi(1, 4) = tmp.y();
      Phi(2, 4) = tmp.z();
    }
    {
      // dr/dlambda
      Vector3d tmp = velocity * (1./n);
      Phi(0, 5) = tmp.x();
      Phi(1, 5) = tmp.y();
      Phi(2, 5) = tmp.z();
    }

    // Velocity derivatives
    {
      // drd/da
      Vector3d tmp = - (1./(2*a)) * (velocity - position * 3 * GM * (time -tau).seconds() / (r*r*r) );
      Phi(3, 0) = tmp.x();
      Phi(4, 0) = tmp.y();
      Phi(5, 0) = tmp.z();
    }
    {
      // drd/dh
      Vector3d tmp = dXDdh * f + dYDdh * g;
      Phi(3, 1) = tmp.x();
      Phi(4, 1) = tmp.y();
      Phi(5, 1) = tmp.z();
    }
    {
      // drd/dk
      Vector3d tmp = dXDdk * f + dYDdk * g;
      Phi(3, 2) = tmp.x();
      Phi(4, 2) = tmp.y();
      Phi(5, 2) = tmp.z();
    }
    {
      // drd/dp
      Vector3d tmp = 2.* (q * (YD * f - XD * g) - XD * w) / C;
      Phi(3, 3) = tmp.x();
      Phi(4, 3) = tmp.y();
      Phi(5, 3) = tmp.z();
    }
    {
      // drd/dq
      Vector3d tmp = 2.* (p * (XD * g - YD * f) + YD * w) / C;
      Phi(3, 4) = tmp.x();
      Phi(4, 4) = tmp.y();
      Phi(5, 4) = tmp.z();
    }
    {
      // drd/dlambda
      Vector3d tmp = - n * std::pow(a / r, 3) * position;
      Phi(3, 5) = tmp.x();
      Phi(4, 5) = tmp.y();
      Phi(5, 5) = tmp.z();
    }

    return Phi;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

const Vector3d Equinoctial::position(const Time &time) const
{
  try
  {
    Vector3d position, velocity;
    orbit(time, position, velocity);
    return position;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Equinoctial::orbit(const Time &time, Vector3d &position, Vector3d &velocity) const
{
  try
  {
    Vector3d acceleration;
    orbit(time, position, velocity, acceleration);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Equinoctial::orbit(const Time &time, Vector3d &position, Vector3d &velocity, Vector3d &acceleration) const
{
  try
  {
    Vector3d inPlanePosition, inPlaneVelocity;
    inPlaneOrbit(time, inPlanePosition, inPlaneVelocity);

    Vector3d of, og;
    orthonormals(of, og);

    position = inPlanePosition.x() * of + inPlanePosition.y() * og;
    velocity = inPlaneVelocity.x() * of + inPlaneVelocity.y() * og;
    acceleration = -GM / (std::pow(position.norm(), 3)) * position;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Equinoctial::inPlaneOrbit(const Time &time, Vector3d &position, Vector3d &velocity) const
{
  try
  {
    Vector3d acceleration;
    inPlaneOrbit(time, position, velocity, acceleration);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Equinoctial::inPlaneOrbit(const Time &time, Vector3d &position, Vector3d &velocity, Vector3d &acceleration) const
{
  try
  {
    Double F = eccentricLongitude(time);

    Double beta = 1 / (1 + std::sqrt(1 - h * h - k * k)); // Supplementary quantity
    Double r    = a * (1 - h * std::sin(F) - k * std::cos(F)); // Radius

    // Coordinates and velocities in orbital plane
    position.x() = a * ((1 - h * h * beta) * std::cos(F) + h * k * beta * std::sin(F) - k);
    position.y() = a * ((1 - k * k * beta) * std::sin(F) + h * k * beta * std::cos(F) - h);
    position.z() = 0;

    velocity.x() = std::sqrt(GM * a) / r * ( h * k * beta * std::cos(F) - (1 - h * h * beta) * std::sin(F));
    velocity.y() = std::sqrt(GM * a) / r * (-h * k * beta * std::sin(F) + (1 - k * k * beta) * std::cos(F));
    velocity.z() = 0;

    acceleration = -GM / (r * r * r) * position;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
