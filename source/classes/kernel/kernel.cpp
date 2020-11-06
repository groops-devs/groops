/***********************************************/
/**
* @file kernel.cpp
*
* @brief Isotropic harmonic integral kernels.
*
* @author Torsten Mayer-Guerr
* @date 2003-09-03
*
*/
/***********************************************/

#define DOCSTRING_Kernel

#include "base/import.h"
#include "base/legendrePolynomial.h"
#include "config/configRegister.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/kernel/kernelPoisson.h"
#include "classes/kernel/kernelStokes.h"
#include "classes/kernel/kernelHotine.h"
#include "classes/kernel/kernelSingleLayer.h"
#include "classes/kernel/kernelWaterHeight.h"
#include "classes/kernel/kernelDeformation.h"
#include "classes/kernel/kernelCoefficients.h"
#include "classes/kernel/kernelGeoid.h"
#include "classes/kernel/kernelRadialGradient.h"
#include "classes/kernel/kernelFilterGauss.h"
#include "classes/kernel/kernelSelenoid.h"
#include "classes/kernel/kernelBottomPressure.h"
#include "classes/kernel/kernelBlackmanLowPass.h"
#include "classes/kernel/kernelTruncation.h"
#include "classes/kernel/kernel.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Kernel, "kernelType",
                      KernelGeoid,
                      KernelStokes,
                      KernelHotine,
                      KernelPoisson,
                      KernelSingleLayer,
                      KernelWaterHeight,
                      KernelBottomPressure,
                      KernelDeformation,
                      KernelRadialGradient,
                      KernelCoefficients,
                      KernelFilterGauss,
                      KernelBlackmanLowPass,
                      KernelTruncation,
                      KernelSelenoid)

GROOPS_READCONFIG_CLASS(Kernel, "kernelType")

/***********************************************/

KernelPtr Kernel::create(Config &config, const std::string &name)
{
  try
  {
    KernelPtr   kernel;
    std::string type;
    readConfigChoice(config, name, type, Config::MUSTSET, "", "harmonic and isotropic integral kernels");

    if(readConfigChoiceElement(config, "geoidHeight",    type, "geoid = potential/normalgravity"))
      kernel = KernelPtr(new KernelGeoid(config));
    if(readConfigChoiceElement(config, "anomalies",      type, "gravity anomalies, Stokes kernel"))
      kernel = KernelPtr(new KernelStokes(config));
    if(readConfigChoiceElement(config, "disturbance",    type, "gravity disturbance, Hotine kernel"))
      kernel = KernelPtr(new KernelHotine(config));
    if(readConfigChoiceElement(config, "potential",      type, "Abel-Poisson kernel"))
      kernel = KernelPtr(new KernelPoisson(config));
    if(readConfigChoiceElement(config, "density",        type, "mass on a single layer (1/l kernel)"))
      kernel = KernelPtr(new KernelSingleLayer(config));
    if(readConfigChoiceElement(config, "waterHeight",    type, "equivalent water columns, accounts the loading"))
      kernel = KernelPtr(new KernelWaterHeight(config));
    if(readConfigChoiceElement(config, "bottomPressure", type, "ocean bottom pressure in Pascal, accounts for loading"))
      kernel = KernelPtr(new KernelBottomPressure(config));
    if(readConfigChoiceElement(config, "deformation",    type, "radial deformation by loading"))
      kernel = KernelPtr(new KernelDeformation(config));
    if(readConfigChoiceElement(config, "radialGradient", type, "radial gravity gradient"))
      kernel = KernelPtr(new KernelRadialGradient(config));
    if(readConfigChoiceElement(config, "coefficients",   type, "construct kernel by a legendre-polynomial expansion"))
      kernel = KernelPtr(new KernelCoefficients(config));
    if(readConfigChoiceElement(config, "filterGauss",    type, "smoothing by a gauss filter"))
      kernel = KernelPtr(new KernelFilterGauss(config));
    if(readConfigChoiceElement(config, "blackmanLowPass",type, "Blackman low-pass filter"))
      kernel = KernelPtr(new KernelBlackmanLowPass(config));
    if(readConfigChoiceElement(config, "truncation",     type, "truncate kernel at specific degree"))
      kernel = KernelPtr(new KernelTruncation(config));
    if(readConfigChoiceElement(config, "selenoidHeight", type, "selenoid = potential/normalgravity"))
      kernel = KernelPtr(new KernelSelenoid(config));
    endChoice(config);

    return kernel;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Kernel::kernel(Vector3d const &p, Vector3d const &q) const
{
  try
  {
    return kernel(p, q, coefficients(p, maxDegree()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Kernel::radialDerivative(Vector3d const &p, Vector3d const &q) const
{
  try
  {
    return radialDerivative(p, q, coefficients(p, maxDegree()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d Kernel::gradient(Vector3d const &p, Vector3d const &q) const
{
  try
  {
    return gradient(p, q, coefficients(p, maxDegree()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Tensor3d Kernel::gradientGradient(Vector3d const &p, Vector3d const &q) const
{
  try
  {
    return gradientGradient(p, q, coefficients(p, maxDegree()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Kernel::inverseKernel(Vector3d const &p, Vector3d const &q, const Kernel &kernel2) const
{
  try
  {
    const Double r = p.r();
    const Double R = q.r();
    const Double t = inner(p, q)/r/R; // t = cos(psi)

    const Vector k2     = kernel2.coefficients(p, kernel2.maxDegree());
    const Vector k1     = inverseCoefficients (p, k2.size()-1);
    const UInt   degree = std::min(k1.rows(), k2.rows())-1;
    Double       *p1 = k1.field();
    const Double *p2 = k2.field();
    Double  f1 = R/r;
    const Double  f2 = R/r;

    for(UInt n=0; n<=degree; n++)
    {
      // k1(n) *= (R/r)^(n+1) * sqrt(2n+1) * k2(n));
      *p1++ *= f1 * sqrt(2*n+1.0) * *p2++;
      f1  *= f2;
    }

    return LegendrePolynomial::sum(t, k1, degree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Kernel::inverseKernel(const Time &time, const Vector3d &p, const GravityfieldBase &field) const
{
  try
  {
    SphericalHarmonics harmonics = field.sphericalHarmonics(time, maxDegree());
    Vector kn = inverseCoefficients(p, harmonics.maxDegree(), harmonics.isInterior());
    // Convolution with the kernel
    return inner(kn, harmonics.Yn(p));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Double Kernel::kernel(Vector3d const &p, Vector3d const &q, const Vector &kn) const
{
  const Double r = p.r();
  const Double R = q.r();
  const Double t = inner(p, q)/r/R; // t = cos(psi)
  // Factors: radial = sqrt(2n+1)*(R/r)^(n+1) * k_n
  const Vector radial = computeFactors(r, R, kn);
  // K = sum_n sqrt(2n+1)*(R/r)^(n+1) * k_n * P_n(t)
  return LegendrePolynomial::sum(t, radial, kn.size()-1);
}

/***********************************************/

Double Kernel::radialDerivative(Vector3d const &p, Vector3d const &q, const Vector &kn) const
{
  const Double r = p.r();
  const Double R = q.r();
  const Double t = inner(p, q)/r/R; // t = cos(psi)

  // radial_n = -(n+1)/r*(R/r)^(n+1) * sqrt(2n+1) * k_n
  const Vector radial = computeFactorsRadialDerivative(r, R, kn);
  // K = sum_n -(n+1)/r*(R/r)^(n+1) * sqrt(2n+1) * k_n * P_n(t)
  return LegendrePolynomial::sum(t, radial, kn.size()-1);
}

/***********************************************/

Vector3d Kernel::gradient(Vector3d const &p, Vector3d const &q, const Vector &kn) const
{
  const Double r  = p.r();
  const Double r2 = r*r;
  const Double R  = q.r();
  const Double t  = inner(p, q)/r/R; // t = cos(psi)
  const UInt   degree = kn.size()-1;

  const Vector radial           = computeFactors(r, R, kn);
  const Vector radialDerivative = computeFactorsRadialDerivative(r, R, kn);

  // derivatives of r with respect to x,y,z
  const Double dr_dx = p.x()/r;
  const Double dr_dy = p.y()/r;
  const Double dr_dz = p.z()/r;

  // derivatives of the kernel with respect to t,r
  const Double dK_dr    = LegendrePolynomial::sum(t, radialDerivative, degree);
  const Double dK_dt    = LegendrePolynomial::sumDerivative(t, radial, degree);

  // derivatives of t with respect to x,y,z
  const Double dt_dx = q.x()/r/R-p.x()*t/r2;
  const Double dt_dy = q.y()/r/R-p.y()*t/r2;
  const Double dt_dz = q.z()/r/R-p.z()*t/r2;

  // chain rule
  return Vector3d(dK_dr*dr_dx + dK_dt*dt_dx,
                  dK_dr*dr_dy + dK_dt*dt_dy,
                  dK_dr*dr_dz + dK_dt*dt_dz);
}

/***********************************************/

Tensor3d Kernel::gradientGradient(Vector3d const &p, Vector3d const &q, const Vector &kn) const
{
  const Double r  = p.r();
  const Double r2 = r*r;
  const Double r3 = r2*r;
  const Double r4 = r3*r;
  const Double R  = q.r();
  const Double t  = inner(p, q)/r/R; // t = cos(psi)
  const UInt   degree = kn.size()-1;

  const Vector radial              = computeFactors                   (r, R, kn);
  const Vector radialDerivative    = computeFactorsRadialDerivative   (r, R, kn);
  const Vector radialDerivative2nd = computeFactorsRadialDerivative2nd(r, R, kn);

  // derivatives of r with respect to x,y,z
  const Double dr_dx = p.x()/r;
  const Double dr_dy = p.y()/r;
  const Double dr_dz = p.z()/r;

  // 2nd derivatives of r with respect to x,y,z
  const Double d2r_dx2 = 1/r-p.x()*p.x()/r3;
  const Double d2r_dy2 = 1/r-p.y()*p.y()/r3;
  const Double d2r_dz2 = 1/r-p.z()*p.z()/r3;
  const Double d2r_dxdy =   -p.x()*p.y()/r3;
  const Double d2r_dxdz =   -p.x()*p.z()/r3;
  const Double d2r_dydz =   -p.y()*p.z()/r3;

  // derivatives of t with respect to x,y,z
  const Double dt_dx = q.x()/(r*R)-p.x()*t/r2;
  const Double dt_dy = q.y()/(r*R)-p.y()*t/r2;
  const Double dt_dz = q.z()/(r*R)-p.z()*t/r2;

  // 2nd derivatives of t with respect to x,y,z
  const Double d2t_dx2  = -t/r2 - 2*q.x()*p.x()/r3/R + 3*p.x()*p.x()*t/r4;
  const Double d2t_dy2  = -t/r2 - 2*q.y()*p.y()/r3/R + 3*p.y()*p.y()*t/r4;
  const Double d2t_dz2  = -t/r2 - 2*q.z()*p.z()/r3/R + 3*p.z()*p.z()*t/r4;
  const Double d2t_dxdy = -(q.x()*p.y()+q.y()*p.x())/r3/R+3*p.x()*p.y()*t/r4;
  const Double d2t_dxdz = -(q.x()*p.z()+q.z()*p.x())/r3/R+3*p.x()*p.z()*t/r4;
  const Double d2t_dydz = -(q.y()*p.z()+q.z()*p.y())/r3/R+3*p.y()*p.z()*t/r4;

  // derivatives of the kernel with respect to t,r
  // K = sum B_n*(R/r)^n+1 * P_n
  const Double dK_dt    = LegendrePolynomial::sumDerivative(t, radial, degree);
  const Double d2K_drdt = LegendrePolynomial::sumDerivative(t, radialDerivative, degree);
  const Double dK_dr    = LegendrePolynomial::sum          (t, radialDerivative, degree);
  const Double d2K_dr2  = LegendrePolynomial::sum          (t, radialDerivative2nd, degree);
  const Double d2K_dt2  = LegendrePolynomial::sumDerivative2nd(t, radial, degree);

  // chain rule
  Tensor3d tns;
  tns.xx() = dK_dr*d2r_dx2  + dK_dt*d2t_dx2
           + d2K_dr2*dr_dx*dr_dx + 2*d2K_drdt*dr_dx*dt_dx + d2K_dt2*dt_dx*dt_dx;

  tns.xy() = dK_dr*d2r_dxdy + dK_dt*d2t_dxdy
           + d2K_dr2*dr_dx*dr_dy + d2K_drdt*dr_dx*dt_dy + d2K_drdt*dr_dy*dt_dx + d2K_dt2*dt_dx*dt_dy;

  tns.xz() = dK_dr*d2r_dxdz + dK_dt*d2t_dxdz
           + d2K_dr2*dr_dx*dr_dz + d2K_drdt*dr_dx*dt_dz + d2K_drdt*dr_dz*dt_dx + d2K_dt2*dt_dx*dt_dz;

  tns.yy() = dK_dr*d2r_dy2  + dK_dt*d2t_dy2
           + d2K_dr2*dr_dy*dr_dy + 2*d2K_drdt*dr_dy*dt_dy + d2K_dt2*dt_dy*dt_dy;

  tns.yz() = dK_dr*d2r_dydz + dK_dt*d2t_dydz
           + d2K_dr2*dr_dy*dr_dz + d2K_drdt*dr_dy*dt_dz + d2K_drdt*dr_dz*dt_dy + d2K_dt2*dt_dy*dt_dz;

  tns.zz() = dK_dr*d2r_dz2  + dK_dt*d2t_dz2
           + d2K_dr2*dr_dz*dr_dz + 2*d2K_drdt*dr_dz*dt_dz + d2K_dt2*dt_dz*dt_dz;

  return tns;
}

/***********************************************/
/***********************************************/

// computes sqrt(2n+1)*(R/r)^(n+1) * k_n
Vector Kernel::computeFactors(Double r, Double R, const Vector &kn) const
{
  Vector  radial(kn.size());
  Double  f1   = R/r;
  Double  f2   = R/r;
  Double *rptr = radial.field();
  const Double *kptr = kn.field();

  for(UInt n=0; n<kn.size(); n++)
  {
    *rptr++ = sqrt(2.*n+1.) * f1 * *kptr++;
    f1 *= f2;
  }
  return radial;
}

/***********************************************/

// computes dK/dr = -(n+1)/r*(R/r)^(n+1) * sqrt(2n+1) * k_n
Vector Kernel::computeFactorsRadialDerivative(Double r, Double R, const Vector &kn) const
{
  Vector  radial(kn.size());
  Double  f1   = R/(r*r);
  Double  f2   = R/r;
  Double *rptr = radial.field();
  const Double *kptr = kn.field();

  for(UInt n=0; n<kn.size(); n++)
  {
    *rptr++ = - sqrt(2.*n+1.) * (n+1.) * f1 * *kptr++;
    f1 *= f2;
  }
  return radial;
}

/***********************************************/

// computes d2K/dr2 = (n+1)*(n+2)/r^2*(R/r)^(n+1) * sqrt(2n+1) * k_n
Vector Kernel::computeFactorsRadialDerivative2nd(Double r, Double R, const Vector &kn) const
{
  Vector  radial(kn.size());
  Double  f1   = R/(r*r*r);
  Double  f2   = R/r;
  Double *rptr = radial.field();
  const Double *kptr = kn.field();

  for(UInt n=0; n<kn.size(); n++)
  {
    *rptr++ = sqrt(2.*n+1.) * (n+1.)*(n+2.) * f1 * *kptr++;
    f1 *= f2;
  }
  return radial;
}

/***********************************************/
