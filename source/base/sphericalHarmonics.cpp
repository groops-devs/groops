/***********************************************/
/**
* @file sphericalHarmonics.cpp
*
* @brief Spherical harmonics and functions of the gravity field.
*
* (4Pi normalized).
*
* @author Torsten Mayer-Guerr
* @date 2001-05-31
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/tensor3d.h"
#include "base/rotary3d.h"
#include "base/sphericalHarmonics.h"

/***********************************************/

Matrix SphericalHarmonics::factor1;
Matrix SphericalHarmonics::factor2;

/***********************************************/

void SphericalHarmonics::computeFactors(UInt degree)
{
  // Enough factors?
  if(factor1.rows()>degree)
    return;

  factor1 = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
  factor2 = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);

  // factors for the recursion P[n-1][n-1] -> P[n][n]
  if(degree>0) factor1(1,1) = std::sqrt(3.);
  for(UInt n=2; n<=degree; n++)
    factor1(n,n) =  std::sqrt((2.*n+1.)/(2.*n));

  // factors for the recursion P[m][n-1] and P[m][n-2] -> P[m][n]
  for(UInt m=0; m<degree; m++)
    for(UInt n=m+1; n<=degree; n++)
    {
      Double f = (2.*n+1.)/static_cast<Double>((n+m)*(n-m));
      factor1(n,m) =  std::sqrt(f*(2.*n-1.));
      factor2(n,m) = -std::sqrt(f*(n-m-1.)*(n+m-1.)/(2.*n-3.));
    }
}

/***********************************************/

// Basis functions Ynm (Cnm and Snm)
void SphericalHarmonics::CnmSnm(const Vector3d &point, UInt degree, Matrix &Cnm, Matrix &Snm, Bool interior)
{
  computeFactors(degree);

  Cnm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
  Snm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);

  Double rr, x, y, z;
  if(!interior)
  {
    rr = pow(1/point.r(),2);
    x  = point.x() * rr;
    y  = point.y() * rr;
    z  = point.z() * rr;
    Cnm(0,0) = 1e280/point.r(); // dirty trick: to account for small numbers in very high degrees.
  }
  else
  {
    rr = pow(point.r(),2);
    x  = point.x();
    y  = point.y();
    z  = point.z();
    Cnm(0,0) = 1e280; // dirty trick: to account for small numbers in very high degrees.
  }

  // Recursion diagonal
  // C(n-1,n-1) -> C(n,n)
  for(UInt n=1; n<=degree; n++)
  {
    Cnm(n,n) = factor1(n,n) * (x * Cnm(n-1,n-1) - y * Snm(n-1,n-1));
    Snm(n,n) = factor1(n,n) * (y * Cnm(n-1,n-1) + x * Snm(n-1,n-1));
  }

  // Recursion secondary diagonal
  // C(n-1,n-1) -> C(n,n-1)
  for(UInt n=1; n<=degree; n++)
  {
    Cnm(n,n-1) = factor1(n,n-1) * z * Cnm(n-1,n-1);
    Snm(n,n-1) = factor1(n,n-1) * z * Snm(n-1,n-1);
  }

  // Recursion others
  // C(n-1,m),C(n-1,m) -> C(n,m)
  if(degree>1)
  {
    for(UInt m=0; m<(degree-1); m++)
    {
      Double *c  = &Cnm(m+2,m),     *s  = &Snm(m+2,m);
      Double *f1 = &factor1(m+2,m), *f2 = &factor2(m+2,m);

      for(UInt n=m+2; n<=degree; n++)
      {
        *c = *f1 * z * c[-1] + *f2 * rr * c[-2];
        *s = *f1 * z * s[-1] + *f2 * rr * s[-2];
        c++; s++; f1++; f2++;
      }
    }
  }

  Cnm *= 1e-280; // dirty trick: to account for small numbers in very high degrees.
  Snm *= 1e-280;
}

/***********************************************/

Matrix SphericalHarmonics::Pnm(Angle theta, Double _r, UInt degree, Bool interior)
{
  computeFactors(degree);

  Matrix Pnm(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);

  Double r;
  if(!interior)
  {
    r = 1/_r;
    Pnm(0,0) = 1e280*r; // dirty trick: to account for small numbers in very high degrees.
  }
  else
  {
    r = _r;
    Pnm(0,0) = 1e280; // dirty trick: to account for small numbers in very high degrees.
  }
  const Double rr = pow(r,2);
  const Double x  = sin(theta) * r;
  const Double z  = cos(theta) * r;

  // Recursion diagonal
  // C(n-1,n-1) -> C(n,n)
  for(UInt n=1; n<=degree; n++)
    Pnm(n,n) = factor1(n,n) * (x * Pnm(n-1,n-1));

  // Recursion secondary diagonal
  // C(n-1,n-1) -> C(n,n-1)
  for(UInt n=1; n<=degree; n++)
    Pnm(n,n-1) = factor1(n,n-1) * z * Pnm(n-1,n-1);

  // Recursion others
  // C(n-1,m),C(n-1,m) -> C(n,m)
  if(degree>1)
  {
    for(UInt m=0; m<(degree-1); m++)
    {
      Double *c  = &Pnm(m+2,m);
      Double *f1 = &factor1(m+2,m);
      Double *f2 = &factor2(m+2,m);

      for(UInt n=m+2; n<=degree; n++)
      {
        *c = *f1 * z * c[-1] + *f2 * rr * c[-2];
        c++; f1++; f2++;
      }
    }
  }

  Pnm *= 1e-280; // dirty trick: to account for small numbers in very high degrees.
  return Pnm;
}

/***********************************************/

SphericalHarmonics::SphericalHarmonics(Bool interior) :
  _GM(DEFAULT_GM),
  _R(DEFAULT_R),
  _maxDegree(0),
  _cnm(_maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER),
  _snm(_maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER),
  _interior(interior)
{
}

/***********************************************/

SphericalHarmonics::SphericalHarmonics(Double GM, Double R, const const_MatrixSlice &cnm, const const_MatrixSlice &snm, Bool interior) :
    _GM(GM),
    _R(R),
    _maxDegree(cnm.rows()-1),
    _cnm(cnm),
    _snm(snm),
    _interior(interior)
{
  try
  {
    if(cnm.rows()!=snm.rows())
      throw(Exception("Dimensions of potential coefficients do not agree"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics::SphericalHarmonics(Double GM, Double R, const const_MatrixSlice &cnm, const const_MatrixSlice &snm, const const_MatrixSlice &sigma2cnm, const const_MatrixSlice &sigma2snm, Bool interior) :
    _GM(GM),
    _R(R),
    _maxDegree(cnm.rows()-1),
    _cnm(cnm),
    _snm(snm),
    _sigma2cnm(sigma2cnm),
    _sigma2snm(sigma2snm),
    _interior(interior)
{
  try
  {
    if((cnm.rows() != snm.rows()) ||
       ((sigma2cnm.size() || sigma2snm.size()) && ((cnm.rows() != sigma2cnm.rows()) || (snm.rows() != sigma2snm.rows()))))
      throw(Exception("Dimensions of potential coefficients do not agree"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics SphericalHarmonics::get(UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  if(GM <= 0) GM = _GM;
  if(R  <= 0) R  = _R;
  if(maxDegree == INFINITYDEGREE)
    maxDegree = _maxDegree;

  // Quick return possible?
  if((GM == _GM) && (R == _R) && (minDegree == 0) && (maxDegree == _maxDegree))
    return SphericalHarmonics(GM, R, _cnm, _snm, _sigma2cnm, _sigma2snm, _interior);

  const UInt maxd = std::min(maxDegree, _maxDegree);
  if(minDegree > maxd)
    return SphericalHarmonics(GM, R, Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER), Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER), _interior);

  Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
  Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
  copy(_cnm.slice(minDegree, 0, maxd-minDegree+1, maxd+1), cnm.slice(minDegree, 0, maxd-minDegree+1, maxd+1));
  copy(_snm.slice(minDegree, 0, maxd-minDegree+1, maxd+1), snm.slice(minDegree, 0, maxd-minDegree+1, maxd+1));

  Matrix sigma2cnm, sigma2snm;
  if(_sigma2cnm.size())
  {
    sigma2cnm = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    sigma2snm = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    copy(_sigma2cnm.slice(minDegree,0,maxd-minDegree+1,maxd+1), sigma2cnm.slice(minDegree,0,maxd-minDegree+1,maxd+1));
    copy(_sigma2snm.slice(minDegree,0,maxd-minDegree+1,maxd+1), sigma2snm.slice(minDegree,0,maxd-minDegree+1,maxd+1));
  }

  // adjust factors
  if(R != _R)
  {
    Double factor = (_interior) ? (_GM/GM*pow(R/_R, minDegree+1)) : (_GM/GM*pow(_R/R, minDegree));
    for(UInt n=minDegree; n<=maxd; n++)
    {
      cnm.slice(n, 0, 1, n+1) *= factor;
      snm.slice(n, 0, 1, n+1) *= factor;
      if(sigma2cnm.size()) sigma2cnm.slice(n, 0, 1, n+1) *= factor*factor;
      if(sigma2snm.size()) sigma2snm.slice(n, 0, 1, n+1) *= factor*factor;
      factor *= ((_interior) ? R/_R : _R/R);
    }
  }
  else if(GM != _GM)
  {
    const Double factor = _GM/GM;
    cnm *= factor;
    snm *= factor;
    if(sigma2cnm.size()) sigma2cnm *= factor*factor;
    if(sigma2snm.size()) sigma2snm *= factor*factor;
  }

  return SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm, _interior);
}

/***********************************************/

// Based on
// "On a Fortran procedure for rotating spherical-harmonic coefficients"
// R. H. Gooding, C. A. Wagner, DOI 10.1007/s10569-010-9293-3
SphericalHarmonics SphericalHarmonics::rotate(const Rotary3d &rotary) const
{
  try
  {
    Angle alpha, beta, gamma;
    rotary.euler(alpha, beta, gamma); // rotary = rotaryZ(gamma)*rotaryX(beta)*rotaryZ(alpha);

    Double s = sin(0.5*beta);
    Double c = cos(0.5*beta);
    if(s*s<0.5) // consistency for I/2 is set up
      s = std::sqrt((1.-c)*(1.+c)); // xImod from 2*asin(s)*radian if wanted
    else
      c = std::sqrt((1.-s)*(1.+s)); // xImod from 2*acos(c)*radian if wanted

    // initialize d-matrix degree 1/2 (Risbo comment)
    Matrix d(2,2);
    d(0,0) =  c;
    d(0,1) =  s;
    d(1,0) = -s;
    d(1,1) =  c;

    // MAIN LOOP: combination of Risbo (for key d-file) & Gooding (for rotation)
    Matrix cnmRot(_maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snmRot(_maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    cnmRot(0,0) = _cnm(0,0);
    for(UInt j=2; j<=2*_maxDegree; j++)
    {
      // recursion in two directions from d to dd (Risbo comment)
      Matrix dd(j+1, j+1);
      for(UInt i=0; i<j; i++)
        for(UInt k=0; k<j; k++)
        {
          dd(i,  k)   += std::sqrt((j-i)*(j-k))/j * c * d(i,k);
          dd(i+1,k)   -= std::sqrt((i+1)*(j-k))/j * s * d(i,k);
          dd(i,  k+1) += std::sqrt((j-i)*(k+1))/j * s * d(i,k);
          dd(i+1,k+1) += std::sqrt((i+1)*(k+1))/j * c * d(i,k);
        }
      d = dd;

      // do the rotation (Gooding), ignoring half-integer degrees
      if((j%2)==0)
      {
        const UInt n = j/2;
        for(UInt k=0; k<=n; k++)
          for(UInt m=0; m<=n; m++)
          {
            const Double dplus  = d(n+k, n+m) * ((m==0) ? std::sqrt(0.5) : 1.) * ((k==0) ? std::sqrt(0.5) : 1.);
            const Double dminus = d(n+k, n-m) * ((m==0) ? std::sqrt(0.5) : 1.) * ((k==0) ? std::sqrt(0.5) : 1.) * (((m%2)==1) ? -1. : 1.);
            const Double cm     = std::cos(m*(alpha-PI/2));
            const Double sm     = std::sin(m*(alpha-PI/2));
            const Double ck     = std::cos(k*(-gamma-PI/2));
            const Double sk     = std::sin(k*(-gamma-PI/2));
            cnmRot(n,k) += (ck*cm*(dplus+dminus)+sk*sm*(dplus-dminus)) * _cnm(n,m) + (ck*sm*(dplus+dminus)-sk*cm*(dplus-dminus)) * _snm(n,m);
            snmRot(n,k) += (sk*cm*(dplus+dminus)-ck*sm*(dplus-dminus)) * _cnm(n,m) + (sk*sm*(dplus+dminus)+ck*cm*(dplus-dminus)) * _snm(n,m);
          }
      }
    } // for(j)

    return SphericalHarmonics(_GM, _R, cnmRot, snmRot);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics &SphericalHarmonics::operator+= (const SphericalHarmonics &harm)
{
  try
  {
    if(_interior != harm._interior)
      throw(Exception("Cannot add inner and outer space spherical harmonics"));

    // swap: add smaller to larger matrix and adjust GM and R
    const SphericalHarmonics harm2 = (maxDegree()<harm.maxDegree()) ? get(INFINITYDEGREE, 0, harm.GM(), harm.R()) : harm.get(INFINITYDEGREE, 0, GM(), R());
    if(maxDegree()<harm.maxDegree())
      *this = harm;

    _cnm.slice(0, 0, harm2.maxDegree()+1, harm2.maxDegree()+1) += harm2.cnm();
    _snm.slice(0, 0, harm2.maxDegree()+1, harm2.maxDegree()+1) += harm2.snm();
    if(harm2.sigma2cnm().size() && !_sigma2cnm.size()) _sigma2cnm = Matrix(maxDegree()+1, Matrix::SYMMETRIC, Matrix::LOWER);
    if(harm2.sigma2snm().size() && !_sigma2snm.size()) _sigma2snm = Matrix(maxDegree()+1, Matrix::SYMMETRIC, Matrix::LOWER);
    if(harm2.sigma2cnm().size()) _sigma2cnm.slice(0, 0, harm2.maxDegree()+1, harm2.maxDegree()+1) += harm2.sigma2cnm();
    if(harm2.sigma2snm().size()) _sigma2snm.slice(0, 0, harm2.maxDegree()+1, harm2.maxDegree()+1) += harm2.sigma2snm();

    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics &SphericalHarmonics::operator-= (const SphericalHarmonics &harm)
{
  try
  {
    if(_interior != harm._interior)
      throw(Exception("Cannot subtract inner and outer space spherical harmonics"));

    // swap: add smaller to larger matrix and adjust GM and R
    const SphericalHarmonics harm2 = (maxDegree()<harm.maxDegree()) ? get(INFINITYDEGREE, 0, harm.GM(), harm.R()) : harm.get(INFINITYDEGREE, 0, GM(), R());
    if(maxDegree()<harm.maxDegree())
      *this = harm;

    _cnm.slice(0, 0, harm2.maxDegree()+1, harm2.maxDegree()+1) -= harm2.cnm();
    _snm.slice(0, 0, harm2.maxDegree()+1, harm2.maxDegree()+1) -= harm2.snm();
    if(harm2.sigma2cnm().size() && !_sigma2cnm.size()) _sigma2cnm = Matrix(maxDegree()+1, Matrix::SYMMETRIC, Matrix::LOWER);
    if(harm2.sigma2snm().size() && !_sigma2snm.size()) _sigma2snm = Matrix(maxDegree()+1, Matrix::SYMMETRIC, Matrix::LOWER);
    if(harm2.sigma2cnm().size()) _sigma2cnm.slice(0, 0, harm2.maxDegree()+1, harm2.maxDegree()+1) += harm2.sigma2cnm();
    if(harm2.sigma2snm().size()) _sigma2snm.slice(0, 0, harm2.maxDegree()+1, harm2.maxDegree()+1) += harm2.sigma2snm();

    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics &SphericalHarmonics::operator*= (Double c)
{
  _cnm *= c;
  _snm *= c;
  if(_sigma2cnm.size()) _sigma2cnm *= c*c;
  if(_sigma2snm.size()) _sigma2snm *= c*c;
  return *this;
}

/***********************************************/

// Yn = GM/R * sum (cnm*Cnm + snm*Snm)
Vector SphericalHarmonics::Yn(const Vector3d &point, UInt maxDegree) const
{
  if(maxDegree==INFINITYDEGREE)
   maxDegree = this->maxDegree();

  Matrix Cnm, Snm;
  CnmSnm(1/R() * point, std::min(maxDegree, this->maxDegree()), Cnm, Snm, _interior);

  Vector Y(maxDegree+1);
  Double *Yn = &Y(0);
  Double fak = GM()/R();
  for(UInt n=0; n<=std::min(maxDegree,_maxDegree); n++)
    *(Yn++) = fak * (inner(_cnm.slice(n,0,1,n+1), Cnm.slice(n,0,1,n+1))
                    +inner(_snm.slice(n,1,1,n),   Snm.slice(n,1,1,n)));
  return Y;
}

/***********************************************/

Double SphericalHarmonics::potential(const Vector3d &point, UInt maxDegree, UInt minDegree) const
{
  try
  {
    maxDegree = std::min(maxDegree, this->maxDegree());

    const Vector Y = Yn(point, maxDegree);

    Double sum = 0.0;
    for(UInt n=minDegree; n<=maxDegree; n++)
      sum += Y(n);

    return sum;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SphericalHarmonics::radialGradient(const Vector3d &point, UInt maxDegree, UInt minDegree) const
{
  try
  {
    maxDegree  = std::min(maxDegree, this->maxDegree());
    const Double r = point.r();
    const Vector Y = Yn(point, maxDegree);
    Double sum = 0.0;
    if(!_interior)
      for(UInt n=minDegree; n<=maxDegree; n++)
        sum += -(n+1.)/r * Y(n);
    else
      for(UInt n=minDegree; n<=maxDegree; n++)
        sum += n/r * Y(n);
    return sum;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d SphericalHarmonics::gravity(const Vector3d &point, UInt maxDegree, UInt minDegree) const
{
  try
  {
    if(_interior)
      throw(Exception("not implemented yet for inner space"));

    maxDegree = std::min(maxDegree, this->maxDegree());

    Matrix Cnm, Snm;
    CnmSnm(1/R() * point, maxDegree+1, Cnm, Snm);

    Vector3d g, g_ges;

    // all degrees
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      // 0. Order
      Double wm0 = std::sqrt(static_cast<Double>(n+1)*(n+1));
      Double wp1 = std::sqrt(static_cast<Double>(n+1)*(n+2)) / std::sqrt(2.0);

      Double Cm0 = wm0*Cnm(n+1,0);
      Double Cp1 = wp1*Cnm(n+1,1); Double Sp1 = wp1*Snm(n+1,1);

      g.x() = _cnm(n,0) * (-2*Cp1);
      g.y() = _cnm(n,0) * (-2*Sp1);
      g.z() = _cnm(n,0) * (-2*Cm0);

      // all other orders
      for(UInt m=1; m<=n; m++)
      {
        Double wm1 = std::sqrt(static_cast<Double>(n-m+1)*(n-m+2)) * ((m==1) ? std::sqrt(2.0) : 1.0);
        Double wm0 = std::sqrt(static_cast<Double>(n-m+1)*(n+m+1));
        Double wp1 = std::sqrt(static_cast<Double>(n+m+1)*(n+m+2));

        Double Cm1 = wm1*Cnm(n+1,m-1);  Double Sm1 = wm1*Snm(n+1,m-1);
        Double Cm0 = wm0*Cnm(n+1,m  );  Double Sm0 = wm0*Snm(n+1,m  );
        Double Cp1 = wp1*Cnm(n+1,m+1);  Double Sp1 = wp1*Snm(n+1,m+1);

        g.x() += _cnm(n,m) * ( Cm1 - Cp1) + _snm(n,m) * (Sm1 - Sp1);
        g.y() += _cnm(n,m) * (-Sm1 - Sp1) + _snm(n,m) * (Cm1 + Cp1);
        g.z() += _cnm(n,m) * (-2*Cm0)     + _snm(n,m) * (-2*Sm0);
      } // for(m)

      g_ges  += std::sqrt((2.*n+1.)/(2.*n+3.))*g;
    } // for(n)

    return GM()/(2*R()*R()) * g_ges;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Tensor3d SphericalHarmonics::gravityGradient(const Vector3d &point, UInt maxDegree, UInt minDegree) const
{
  try
  {
    if(_interior)
      throw(Exception("not implemented yet for inner space"));

    maxDegree = std::min(maxDegree, this->maxDegree());

    Matrix Cnm, Snm;
    CnmSnm(1/R() * point, maxDegree+2, Cnm, Snm);

    Tensor3d K, K_ges;

    // all degrees
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      // 0. Order
      Double wm0 = std::sqrt(static_cast<Double>(n+1)*(n+2)*(n+1)*(n+2));
      Double wp1 = std::sqrt(static_cast<Double>(n+1)*(n+1)*(n+2)*(n+3)) / std::sqrt(2.0);
      Double wp2 = std::sqrt(static_cast<Double>(n+1)*(n+2)*(n+3)*(n+4)) / std::sqrt(2.0);

      Double Cm0 = wm0*Cnm(n+2,0);
      Double Cp1 = wp1*Cnm(n+2,1);  Double Sp1 = wp1*Snm(n+2,1);
      Double Cp2 = wp2*Cnm(n+2,2);  Double Sp2 = wp2*Snm(n+2,2);

      K.xx() = _cnm(n,0) * (-2*Cm0 + 2*Cp2);
      K.xy() = _cnm(n,0) * ( 2*Sp2);
      K.xz() = _cnm(n,0) * ( 4*Cp1);
      K.yy() = _cnm(n,0) * (-2*Cm0 - 2*Cp2);
      K.yz() = _cnm(n,0) * ( 4*Sp1);
      K.zz() = _cnm(n,0) * ( 4*Cm0);

      // 1. order
      if(n>0)
      {
        UInt m=1;
        Double wm1 = std::sqrt(static_cast<Double>(n-m+1)*(n-m+2)*(n-m+3)*(n+m+1)) * std::sqrt(2.0);
        Double wm0 = std::sqrt(static_cast<Double>(n-m+1)*(n-m+2)*(n+m+1)*(n+m+2));
        Double wp1 = std::sqrt(static_cast<Double>(n-m+1)*(n+m+1)*(n+m+2)*(n+m+3));
        Double wp2 = std::sqrt(static_cast<Double>(n+m+1)*(n+m+2)*(n+m+3)*(n+m+4));

        Double Cm1 = wm1*Cnm(n+2,m-1);  Double Sm1 = wm1*Snm(n+2,m-1);
        Double Cm0 = wm0*Cnm(n+2,m  );  Double Sm0 = wm0*Snm(n+2,m  );
        Double Cp1 = wp1*Cnm(n+2,m+1);  Double Sp1 = wp1*Snm(n+2,m+1);
        Double Cp2 = wp2*Cnm(n+2,m+2);  Double Sp2 = wp2*Snm(n+2,m+2);

        K.xx() += _cnm(n,m) * (- 3*Cm0 + Cp2)  + _snm(n,m) * (- Sm0 + Sp2);
        K.xy() += _cnm(n,m) * (-   Sm0 + Sp2)  + _snm(n,m) * (- Cm0 - Cp2);
        K.xz() += _cnm(n,m) * (-2*Cm1 + 2*Cp1) + _snm(n,m) * (-2*Sm1 + 2*Sp1);
        K.yy() += _cnm(n,m) * (-   Cm0 - Cp2)  + _snm(n,m) * (- 3*Sm0 - Sp2);
        K.yz() += _cnm(n,m) * (2*Sp1)          + _snm(n,m) * (-2*Cm1 - 2*Cp1);
        K.zz() += _cnm(n,m) * (4*Cm0)          + _snm(n,m) * (4*Sm0);
      } // end 1. order

      // all other orders
      for(UInt m=2; m<=n; m++)
      {
        Double wm2 = std::sqrt(static_cast<Double>(n-m+1)*(n-m+2)*(n-m+3)*(n-m+4)) * ((m==2) ? std::sqrt(2.0) : 1.0);
        Double wm1 = std::sqrt(static_cast<Double>(n-m+1)*(n-m+2)*(n-m+3)*(n+m+1));
        Double wm0 = std::sqrt(static_cast<Double>(n-m+1)*(n-m+2)*(n+m+1)*(n+m+2));
        Double wp1 = std::sqrt(static_cast<Double>(n-m+1)*(n+m+1)*(n+m+2)*(n+m+3));
        Double wp2 = std::sqrt(static_cast<Double>(n+m+1)*(n+m+2)*(n+m+3)*(n+m+4));

        Double Cm2 = wm2*Cnm(n+2,m-2);  Double Sm2 = wm2*Snm(n+2,m-2);
        Double Cm1 = wm1*Cnm(n+2,m-1);  Double Sm1 = wm1*Snm(n+2,m-1);
        Double Cm0 = wm0*Cnm(n+2,m  );  Double Sm0 = wm0*Snm(n+2,m  );
        Double Cp1 = wp1*Cnm(n+2,m+1);  Double Sp1 = wp1*Snm(n+2,m+1);
        Double Cp2 = wp2*Cnm(n+2,m+2);  Double Sp2 = wp2*Snm(n+2,m+2);

        K.xx() += _cnm(n,m) * ( Cm2 - 2*Cm0 + Cp2) + _snm(n,m) * ( Sm2 - 2*Sm0 + Sp2);
        K.xy() += _cnm(n,m) * (-Sm2         + Sp2) + _snm(n,m) * ( Cm2         - Cp2);
        K.xz() += _cnm(n,m) * (-2*Cm1 + 2*Cp1)     + _snm(n,m) * (-2*Sm1 + 2*Sp1);
        K.yy() += _cnm(n,m) * (-Cm2 - 2*Cm0 - Cp2) + _snm(n,m) * (-Sm2 - 2*Sm0 - Sp2);
        K.yz() += _cnm(n,m) * ( 2*Sm1 + 2*Sp1)     + _snm(n,m) * (-2*Cm1 - 2*Cp1);
        K.zz() += _cnm(n,m) * (4*Cm0)              + _snm(n,m) * (4*Sm0);
      }  // for(m)

      K_ges  += std::sqrt((2.*n+1.)/(2.*n+5.))*K;
    } // for(n)

    return GM()/(4*R()*R()*R()) * K_ges;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d SphericalHarmonics::deformation(const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, UInt maxDegree, UInt minDegree) const
{
  try
  {
    if(_interior)
      throw(Exception("not implemented yet for inner space"));

    maxDegree = std::min(maxDegree, this->maxDegree());
    Vector3d up = normalize(point);

    Matrix Cnm, Snm;
    CnmSnm(1/R()*point, maxDegree+1, Cnm, Snm);

    Vector3d disp;
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      // 0. order
      Double wm0 = std::sqrt((n+1.)*(n+1.));
      Double wp1 = std::sqrt((n+1.)*(n+2.)) / std::sqrt(2.0);
      Double Cm0 = wm0*Cnm(n+1,0);
      Double Cp1 = wp1*Cnm(n+1,1);
      Double Sp1 = wp1*Snm(n+1,1);

      Double   Vn  =  _cnm(n,0) * Cnm(n,0);
      Vector3d gradVn(_cnm(n,0) * (-2*Cp1),
                      _cnm(n,0) * (-2*Sp1),
                      _cnm(n,0) * (-2*Cm0));

      // other orders
      for(UInt m=1; m<=n; m++)
      {
        Double wm1 = std::sqrt((n-m+1.)*(n-m+2.)) * ((m==1) ? std::sqrt(2.0) : 1.0);
        Double wm0 = std::sqrt((n-m+1.)*(n+m+1.));
        Double wp1 = std::sqrt((n+m+1.)*(n+m+2.));
        Double Cm1 = wm1*Cnm(n+1,m-1);  Double Sm1 = wm1*Snm(n+1,m-1);
        Double Cm0 = wm0*Cnm(n+1,m  );  Double Sm0 = wm0*Snm(n+1,m  );
        Double Cp1 = wp1*Cnm(n+1,m+1);  Double Sp1 = wp1*Snm(n+1,m+1);

        Vn         += _cnm(n,m) *  Cnm(n,m)    + _snm(n,m) * Snm(n,m);
        gradVn.x() += _cnm(n,m) * ( Cm1 - Cp1) + _snm(n,m) * (Sm1 - Sp1);
        gradVn.y() += _cnm(n,m) * (-Sm1 - Sp1) + _snm(n,m) * (Cm1 + Cp1);
        gradVn.z() += _cnm(n,m) * (-2*Cm0)     + _snm(n,m) * (-2*Sm0);
      } // for(m)

      Vn     *= GM()/R();
      gradVn *= GM()/(2*R()) * std::sqrt((2*n+1.)/(2*n+3.));

      disp += (hn(n)/gravity*Vn) * up + (ln(n)/gravity) * (gradVn-inner(gradVn,up)*up);
    } // for(n)

    return disp;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Vector SphericalHarmonics::x() const
{
  Vector x((_maxDegree+1)*(_maxDegree+1));
  UInt idx = 0;
  for(UInt n=0; n<=_maxDegree; n++)
  {
    x(idx++) = _cnm(n,0);
    for(UInt m=1; m<=n; m++)
    {
      x(idx++) = _cnm(n,m);
      x(idx++) = _snm(n,m);
    }
  }
  return x;
}

/***********************************************/

Vector SphericalHarmonics::sigma2x() const
{
  Vector x((_maxDegree+1)*(_maxDegree+1));
  if(_sigma2cnm.size() == 0)
    return x;

  UInt idx = 0;
  for(UInt n=0; n<=_maxDegree; n++)
  {
    x(idx++) = _sigma2cnm(n,0);
    for(UInt m=1; m<=n; m++)
    {
      x(idx++) = _sigma2cnm(n,m);
      x(idx++) = _sigma2snm(n,m);
    }
  }
  return x;
}

/***********************************************/
/***********************************************/
