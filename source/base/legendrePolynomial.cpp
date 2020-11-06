/***********************************************/
/**
* @file legendrePolynomial.cpp
*
* @brief Legendre polynomials.
*
* fully normalized.
*
* @author Torsten Mayer-Guerr
* @author Annette Eicker
* @date 2001-05-31
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/matrix.h"
#include "base/legendrePolynomial.h"

/***** VARIABLES ********************************/

// factors for the recursion
Vector LegendrePolynomial::factor1,            LegendrePolynomial::factor2;
Vector LegendrePolynomial::factor1Derivate,    LegendrePolynomial::factor2Derivate;
Vector LegendrePolynomial::factor1Derivate2nd, LegendrePolynomial::factor2Derivate2nd;
Vector LegendrePolynomial::factor1Integral,    LegendrePolynomial::factor2Integral;

/***********************************************/

void LegendrePolynomial::computeFactors(UInt degree)
{
  // Enough factors already computed?
  if(factor1.rows()>degree)
    return;

  factor1 = Vector(degree+1);
  factor2 = Vector(degree+1);

  factor1(0) = 1.0;
  factor2(0) = sqrt(3.0);

  for(UInt n=2; n<=degree; n++)
  {
    factor1(n) = sqrt(2.*n+1)*sqrt(2.*n-1.)/n;
    factor2(n) = -(n-1.)*sqrt(2.*n+1.)/sqrt(2.*n-3.)/n;
  }
}

/***********************************************/

void LegendrePolynomial::computeFactorsDerivate(UInt degree)
{
  // Enough factors already computed?
  if(factor1Derivate.rows()>degree)
    return;

  factor1Derivate = Vector(degree+1);
  factor2Derivate = Vector(degree+1);

  factor1Derivate(0) = sqrt(3.);
  factor2Derivate(0) = 3.*sqrt(5.);

  for(UInt n=2; n<=degree; n++)
  {
    factor1Derivate(n) = sqrt(2.*n+1.)*sqrt(2.*n-1.)/(n-1.);
    factor2Derivate(n) = -sqrt(2.*n+1.)*n/sqrt(2.*n-3.)/(n-1.);
  }
}

/***********************************************/

void LegendrePolynomial::computeFactorsDerivate2nd(UInt degree)
{
  // Enough factors already computed?
  if(factor1Derivate2nd.rows()>degree)
    return;

  factor1Derivate2nd = Vector(degree+1);
  factor2Derivate2nd = Vector(degree+1);

  factor1Derivate2nd(0) =  3. * sqrt(5.);
  factor2Derivate2nd(0) = 15. * sqrt(7.);

  for(UInt n=3; n<=degree; n++)
  {
    factor1Derivate2nd(n) = sqrt(2.*n+1.)*sqrt(2.*n-1.)/(n-2.);
    factor2Derivate2nd(n) = -(n+1.)*sqrt(2.*n+1.)/sqrt(2.*n-3.)/(n-2.);
  }
}

/***********************************************/

void LegendrePolynomial::computeFactorsIntegral(UInt degree)
{
  // Enough factors already computed?
  if(factor1Integral.rows()>degree)
    return;

  factor1Integral = Vector(degree+1);
  factor2Integral = Vector(degree+1);

  for(UInt n=1; n<=degree; n++)
  {
    factor1Integral(n) = -1./sqrt((2.*n+1.)*(2.*n+3.));
    factor2Integral(n) =  1./sqrt((2.*n+1.)*(2.*n-1.));
  }
}

/***********************************************/

const Vector LegendrePolynomial::compute(Double t, UInt degree)
{
  computeFactors(degree);

  Vector P(degree+1);
  P(0) = 1.0;

  if(degree>=1) P(1) = sqrt(3.)*t;

  for(UInt n=2; n<=degree; n++)
    P(n) = factor1(n)*t*P(n-1) + factor2(n)*P(n-2);

  return P;
}

/***********************************************/

const Vector LegendrePolynomial::derivative(Double t, UInt degree)
{
  computeFactorsDerivate(degree);

  Vector P(degree+1);
  if(degree>=1) P(1) = sqrt(3.);

  for(UInt n=2; n<=degree; n++)
    P(n) = factor1Derivate(n)*t*P(n-1) + factor2Derivate(n)*P(n-2);

  return P;
}

/***********************************************/

const Vector LegendrePolynomial::derivative2nd(Double t, UInt degree)
{
  computeFactorsDerivate2nd(degree);

  Vector P(degree+1);
  if(degree>=2) P(2) = 3.*sqrt(5.);

  for(UInt n=3; n<=degree; n++)
    P(n) = factor1Derivate2nd(n)*t*P(n-1) + factor2Derivate2nd(n)*P(n-2);

  return P;
}

/***********************************************/

const Vector LegendrePolynomial::integral(Double t, UInt degree)
{
  computeFactorsIntegral(degree);

  Vector R(degree+1);
  Vector P = compute(t,degree);

  R(0)=1-t;
  for(UInt n=1; n<=degree-1; n++)
    R(n) = factor1Integral(n)*P(n+1) + factor2Integral(n)*P(n-1);

  return R;
}

/***********************************************/

Double LegendrePolynomial::sum(Double t, const Vector &koeff, UInt degree)
{
  computeFactors(degree);

  // pointer arithemtic to be as fast as possible
  const Double *aptr = factor1.field()+degree;
  const Double *bptr = factor2.field()+degree;
  const Double *kptr = koeff.field()+degree;
  const Double *stop = koeff.field();

  Double u;
  Double u2 = *(kptr--);
  Double u1 = *(aptr--) * t * u2 + *(kptr--);

  while(kptr!=stop)
  {
    u  = *(aptr--) * t * u1 + *(bptr--) * u2 + *(kptr--);
    u2 = u1;
    u1 = u;
  }

  return *kptr  + *bptr * u2 + u1 * sqrt(3.) * t;
}

/***********************************************/

Double LegendrePolynomial::sumDerivative(Double t, const Vector &koeff, UInt degree)
{
  computeFactorsDerivate(degree);

  // pointer arithemtic to be as fast as possible
  const Double *aptr = factor1Derivate.field()+degree;
  const Double *bptr = factor2Derivate.field()+degree;
  const Double *kptr = koeff.field()+degree;
  const Double *stop = koeff.field()+1;

  Double u;
  Double u2 = *(kptr--);
  Double u1 = *(aptr--) * t * u2 + *(kptr--);

  while(kptr!=stop)
  {
    u  = *(aptr--) * t * u1 + *(bptr--) * u2 + *(kptr--);
    u2 = u1;
    u1 = u;
  }

  return (*kptr  + *bptr * u2) * sqrt(3.) + u1 * (3.*sqrt(5.)) * t;
}

/***********************************************/

Double LegendrePolynomial::sumDerivative2nd(Double t, const Vector &koeff, UInt degree)
{
  computeFactorsDerivate2nd(degree);

  // pointer arithemtic to be as fast as possible
  const Double *aptr = factor1Derivate2nd.field()+degree;
  const Double *bptr = factor2Derivate2nd.field()+degree;
  const Double *kptr = koeff.field()+degree;
  const Double *stop = koeff.field()+2;

  Double u;
  Double u2 = *(kptr--);
  Double u1 = *(aptr--) * t * u2 + *(kptr--);

  while(kptr!=stop)
  {
    u  = *(aptr--) * t * u1 + *(bptr--) * u2 + *(kptr--);
    u2 = u1;
    u1 = u;
  }

  return (*kptr  + *bptr * u2) * (3.*sqrt(5.)) + u1 * (15.*sqrt(7.)) * t;
}

/***********************************************/

void LegendrePolynomial::zeros(UInt degree, Vector &zeros, Vector &weights)
{
  zeros   = Vector(degree);
  weights = Vector(degree);

  Double pf=0, pd=1;
  for(UInt k=0; k<(degree+1)/2; k++)
  {
    Double x = std::cos(PI*(4*(k+1)-1)/(4*degree+2));
    Double x0;
    do
    {
      x0 = x;
      Double f0 = 1.;
      Double f1 = x;
      for(UInt n=2; n<=degree; n++)
      {
        pf = (2-1./n)*x*f1-(1-1./n)*f0;
        pd = n*(x*pf-f1)/(x*x-1);
        f0 = f1;
        f1 = pf;
      }
      x -= pf/pd;
    }
    while(std::fabs(x-x0) > std::fabs(x)*1e-15);

    zeros(k)            =  x;
    zeros(degree-1-k)   = -x;
    weights(k)          = 2./((1.-x*x)*pd*pd);
    weights(degree-1-k) = weights(k);
  }
}

/***********************************************/
