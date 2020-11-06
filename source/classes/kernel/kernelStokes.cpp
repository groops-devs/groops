/***********************************************/
/**
* @file kernelStokes.cpp
*
* @brief Stokes kernel (gravity anomalies).
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#include "base/import.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelStokes.h"

/***********************************************/

Vector KernelStokes::coefficients(Vector3d const &q, UInt degree) const
{
  if(degree==INFINITYDEGREE)
    throw(Exception("In KernelStokes::coefficients: INFINITYDEGREE requested"));

  Double  R = q.r();
  Vector  k(degree+1);
  Double *kn = k.field();

  for(UInt n=0; n<=degree; n++)
    *kn++ = R * ((n==1) ? (0.0) : (1.0/(n-1.0)) );

  return k;
}

/***********************************************/

Vector KernelStokes::inverseCoefficients(Vector3d const &p, UInt degree, Bool interior) const
{
  if(degree==INFINITYDEGREE)
    throw(Exception("In KernelStokes::inverseCoefficients: INFINITYDEGREE requested"));

  Double  r = p.r();
  Vector  k(degree+1);
  Double *kn = k.field();

  if(interior)
    for(UInt n=0; n<=degree; n++)
      *kn++ = -(n+2.)/r;
  else
    for(UInt n=0; n<=degree; n++)
      *kn++ = (n-1.)/r;

  return k;
}

/***********************************************/

Double KernelStokes::kernel(Vector3d const &p, Vector3d const &q) const
{
  Vector3d diff  = p-q;
  Double r       = p.r();
  Double R       = q.r();
  Double l       = diff.r();
  Double cos_psi = (R*R+r*r-l*l)/(2*R*r);

  return R*(2*R/l-3*R*l/(r*r)-R*R/(r*r)*cos_psi*(5+3*log((l+r-R*cos_psi)/(2*r))));
}

/***********************************************/

Double KernelStokes::radialDerivative(Vector3d const &p, Vector3d const &q) const
{
  Vector3d diff  = p-q;
  Double r       = p.r();
  Double R       = q.r();
  Double l       = diff.r();
  Double cos_psi = (R*R+r*r-l*l)/(2*R*r);

  Double term0   = (R*R+r*r-l*l);
  Double term    = (r+l-R)*(r+l+R);
  Double ln_term = log(term/(4*r*r));
  Double r2      = r*r;
  Double r3      = r2*r;
  Double r4      = r3*r;

  // Ableitungen Stokes nach r,l
  // Kettenregel: ds/dr = dS/du*du/dr + dS/dl*dl/dr  mit u = r
  Double dS_du = R*(6*R*l/r3+1.5*R*term0*(5+3*ln_term)/r4
                 -R*(5+3*ln_term)/r2
                 -3*R*term0*(0.25*term0/r3-0.5*(r-0.5*term0/r+l)/r2)/(r2*(r-0.5*term0/r+l)));

  Double dS_dl = R*(R/r3*(-2*r3/(l*l)-3*r+5*l+3*l*ln_term
                -3*term0*(l+r)/term));

  Double dl_dr = (r-R*cos_psi)/l;

  Double dS_dr = dS_du+dS_dl*dl_dr;

  return dS_dr;
}

/***********************************************/

Vector3d KernelStokes::gradient(Vector3d const &p, Vector3d const &q) const
{
  // Vorausberechnungen
  Vector3d diff  = p-q;
  Double   r     = p.r();
  Double   R     = q.r();
  Double   l     = diff.r();
  Double term0   = (R*R+r*r-l*l);
  Double term    = (r+l-R)*(r+l+R);
  Double ln_term = log(term/(4*r*r));
  Double r2      = r*r;
  Double r3      = r2*r;
  Double r4      = r3*r;

  // Ableitungen r,l nach x,y,z
  Vector3d dr = 1/r * p;
  Vector3d dl = 1/l * diff;

  // Ableitungen Stokes nach r,l
  Double dS_dr = -R/r2+6*R*l/r3+1.5*R/r4*term0*(5+3*ln_term)
                 -R/r2*(5+3*ln_term)
                 -3*R/r3*term0*((r+l)/term-1/r);

  Double dS_dl = R/r3*(-2*r3/(l*l)-3*r+5*l+3*l*ln_term
                -3*term0*(l+r)/term);

  // Kettenregel
  return R*(dS_dr*dr + dS_dl*dl);
}

/***********************************************/

Tensor3d KernelStokes::gradientGradient(Vector3d const &p, Vector3d const &q) const
{
  Tensor3d tns;

  // Vorausberechnungen
  Double r       = p.r();
  Double R       = q.r();
  Double l       = (p-q).r();
  Double term0   = (R*R+r*r-l*l);
  Double term    = (r+l-R)*(r+l+R);
  Double ln_term = log(term/(4*r*r));
  Double r2      = r*r;
  Double r3      = r2*r;
  Double r4      = r3*r;
  Double r5      = r4*r;

  // Ableitungen r,l nach x,y,z
  Double dr_dx = p.x()/r;   Double dl_dx = (p.x()-q.x())/l;
  Double dr_dy = p.y()/r;   Double dl_dy = (p.y()-q.y())/l;
  Double dr_dz = p.z()/r;   Double dl_dz = (p.z()-q.z())/l;

  // Zweite Ableitungen r,l nach x,y,z
  Double d2r_dx2 = 1/r-p.x()*p.x()/r3;
  Double d2r_dy2 = 1/r-p.y()*p.y()/r3;
  Double d2r_dz2 = 1/r-p.z()*p.z()/r3;
  Double d2r_dxdy =   -p.x()*p.y()/r3;
  Double d2r_dxdz =   -p.x()*p.z()/r3;
  Double d2r_dydz =   -p.y()*p.z()/r3;

  Double d2l_dx2 = 1/l-(p.x()-q.x())*(p.x()-q.x())/(l*l*l);
  Double d2l_dy2 = 1/l-(p.y()-q.y())*(p.y()-q.y())/(l*l*l);
  Double d2l_dz2 = 1/l-(p.z()-q.z())*(p.z()-q.z())/(l*l*l);
  Double d2l_dxdy =   -(p.x()-q.x())*(p.y()-q.y())/(l*l*l);
  Double d2l_dxdz =   -(p.x()-q.x())*(p.z()-q.z())/(l*l*l);
  Double d2l_dydz =   -(p.y()-q.y())*(p.z()-q.z())/(l*l*l);

  // Ableitungen Stokes nach r,l

  Double dS_dr = -R/r2+6*R*l/r3+1.5*R*term0*(5+3*ln_term)/r4
                 -R*(5+3*ln_term)/r2
                 -3*R*term0*(0.25*term0/r3
                 -0.5*(r-0.5*term0/r+l)/r2)/(r2*(r-0.5*term0/r+l));

  Double dS_dl = R/r3*(-2*r3/(l*l)-3*r+5*l+3*l*ln_term
                -3*term0*(l+r)/term);

  // Zweite Ableitungen Stokes nach r,l
  Double d2S_dr2 = R*(-3/r3-18*l/r4-30*(R*R-l*l)/r5
                     +(15*r2-18*term0)*ln_term/r5
                     +12*(R*R-r*l-l*l)*(R*R-l*l)/(r5*term)
                     +6*term0*(R*R-r*l-l*l)*(r+l)/(r4*term*term)
                     -3*term0*(3*l*l+2*r*l-3*R*R)/(r5*term));

  Double d2S_drdl= R*(6.0/r3 - 15*l/r4 - 9*l*ln_term/r4
                     + 6*(l*(R*R-l*r-l*l)-r2*(l+r))/(r4*term)
                     + 3*term0*(3*l+2*r)/(r4*term)
                     + 6*(r+l)*(r+l)*term0/(r3*term*term));

  Double d2S_dl2 = R/(r3)*(4*r3/(l*l*l)+5+3*ln_term
                     + (12*l*(r+l)-3*term0)/term
                     + 6*(r+l)*(r+l)*term0/(term*term));

  // Kettenregel
  tns.xx() = dS_dr*d2r_dx2  + dS_dl*d2l_dx2
           + d2S_dr2*dr_dx*dr_dx + 2*d2S_drdl*dr_dx*dl_dx + d2S_dl2*dl_dx*dl_dx;

  tns.yy() = dS_dr*d2r_dy2  + dS_dl*d2l_dy2
           + d2S_dr2*dr_dy*dr_dy + 2*d2S_drdl*dr_dy*dl_dy + d2S_dl2*dl_dy*dl_dy;

  tns.zz() = dS_dr*d2r_dz2  + dS_dl*d2l_dz2
           + d2S_dr2*dr_dz*dr_dz + 2*d2S_drdl*dr_dz*dl_dz + d2S_dl2*dl_dz*dl_dz;

  tns.xy() = dS_dr*d2r_dxdy + dS_dl*d2l_dxdy
           + d2S_dr2*dr_dx*dr_dy + d2S_drdl*dr_dx*dl_dy + d2S_drdl*dr_dy*dl_dx + d2S_dl2*dl_dx*dl_dy;

  tns.xz() = dS_dr*d2r_dxdz + dS_dl*d2l_dxdz
           + d2S_dr2*dr_dx*dr_dz + d2S_drdl*dr_dx*dl_dz + d2S_drdl*dr_dz*dl_dx + d2S_dl2*dl_dx*dl_dz;

  tns.yz() = dS_dr*d2r_dydz + dS_dl*d2l_dydz
           + d2S_dr2*dr_dy*dr_dz + d2S_drdl*dr_dy*dl_dz + d2S_drdl*dr_dz*dl_dy + d2S_dl2*dl_dy*dl_dz;

  return R*tns;
}

/***********************************************/

Double KernelStokes::inverseKernel(Vector3d const &p, Vector3d const &q, const Kernel &kernel) const
{
  // anomalies = -dK/dr - 2K/r
  return -kernel.radialDerivative(p,q) - 2*kernel.kernel(p,q)/p.r();;
}

/***********************************************/

Double KernelStokes::inverseKernel(const Time &time, const Vector3d &p, const GravityfieldBase &field) const
{
  // anomalies = -dK/dr - 2K/r
  return -field.radialGradient(time, p) - 2*field.potential(time, p)/p.r();
}

/***********************************************/
