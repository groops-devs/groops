/***********************************************/
/**
* @file legendreFunction.cpp
*
* @brief Associated Legendre functions.
* (fully normalized).
*
* @author Torsten Mayer-Guerr
* @author Annette Eicker
* @date 2001-05-31
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/matrix.h"
#include "legendreFunction.h"

/***** VARIABLES *******************************/

Matrix LegendreFunction::factor1;
Matrix LegendreFunction::factor2;
Matrix LegendreFunction::factor1Integral;
Matrix LegendreFunction::factor2Integral;
Vector LegendreFunction::factorSmall; //integration for small thetas

/***********************************************/

void LegendreFunction::computeFactors(UInt degree)
{
  // Enough factors already computed?
  if(factor1.rows()>degree)
    return;

  factor1 = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
  factor2 = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);

  // factors for recursion P[n-1][n-1] -> P[n][n]
  if(degree>0) factor1(1,1) = std::sqrt(3.0);
  for(UInt n=2; n<=degree; n++)
    factor1(n,n) = std::sqrt((2.*n+1)/(2.*n));

  // factors for recursion P[m][n-1] and P[m][n-2] -> P[m][n]
  for(UInt m=0; m<degree; m++)
    for(UInt n=m+1; n<=degree; n++)
    {
      Double f = (2.*n+1)/((n+m)*(n-m));
      factor1(n,m) =  std::sqrt(f*(2.*n-1));
      factor2(n,m) = -std::sqrt(f*(n-m-1.)*(n+m-1.)/(2.*n-3));
    }
}

/***********************************************/

void LegendreFunction::computeFactorsIntegral(UInt degree)
{
  // Enough factors already computed?
  if(factor1Integral.rows()>degree)
    return;

  factor1Integral = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
  factor2Integral = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);

  // factors for recursion P[n-1][n-2] and int_P[n-2][n-2] -> P[n][n]
  for(UInt n=2; n<=degree; n++)
  {
    factor1Integral(n,n) =  1./(2.*n+2.)*std::sqrt((2.*n+1.)/(n*(n-1.)));
    factor2Integral(n,n) =  1./(2.*n+2.)*std::sqrt(n*(2.*n+1.)*(2.*n-1.)/(n-1.));
  }

  // factors for recursion P[n-1][m] and int_P[n-2][m] -> P[n][m]
  for(UInt m=0; m<degree; m++)
    for(UInt n=m+1; n<=degree; n++)
    {
      factor1Integral(n,m) = -1.0/(n+1.)*std::sqrt((2.*n+1.)*(2.*n-1.)/((n-m)*(n+m)));
      factor2Integral(n,m) = (n-2.)/(n+1.)*std::sqrt((2.*n+1.)*(n+m-1.)*(n-m-1.)/((2.*n-3.)*(n+m)*(n-m)));
    }

  // factors for Hmain diagonal for small thetas
  factorSmall = Vector(degree+1);

  for(UInt n=3; n<=degree; n++)
  {
    factorSmall(n)=1.0;
    for(UInt k=5; k<=(2*n-1); k+=2)
      factorSmall(n) *= k/(k+1.);
    factorSmall(n) = std::sqrt(factorSmall(n)*(2*n+1));
  }
}

/***********************************************/

const Matrix LegendreFunction::compute(Double t, UInt degree)
{
  computeFactors(degree);

  Matrix Fkt(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);

  Fkt(0,0) = 1e280; // dirty trick: to account for small numbers in very high degrees.

  // recursion P[n-1][n-1] -> P[n][n] (main diagonals)
  for(UInt n=1; n<=degree; n++)
    Fkt(n,n) = factor1(n,n) *std::sqrt(1-t*t)* Fkt(n-1, n-1);

  // recursion P[m][n-1] and P[m][n-2] -> P[m][n]
  // secondary diagonal m=n-1
  for(UInt n=1; n<=degree; n++)
    Fkt(n,n-1) = factor1(n,n-1) * t * Fkt(n-1,n-1);

  // all other functions
  for(UInt m=0; m<(degree-1); m++)
    for(UInt n=m+2; n<=degree; n++)
      Fkt(n,m) = factor1(n,m)*t*Fkt(n-1,m) + factor2(n,m)*Fkt(n-2,m);

  return 1e-280 * Fkt; // dirty trick: to account for small numbers in very high degrees.
}

/***********************************************/

const Matrix LegendreFunction::integral(Double t1, Double t2, UInt degree)
{
  computeFactorsIntegral(degree);

  Matrix intP(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);

  Double qt1    = t1*t1;
  Double qt2    = t2*t2;
  Double y1     = std::sqrt(1-qt1);
  Double y2     = std::sqrt(1-qt2);
  Double qy1    = 1-qt1;
  Double qy2    = 1-qt2;
  Double theta1 = acos(t1);
  Double theta2 = acos(t2);

  Matrix P1 = LegendreFunction::compute(t1, degree);
  Matrix P2 = LegendreFunction::compute(t2, degree);

  intP(0,0) = t2-t1;
  intP(1,1) = std::sqrt(3.0)/2.0  * (t2*y2-theta2-t1*y1+theta1);
  intP(2,2) = std::sqrt(5.0/12.0) * (3*t2-std::pow(t2, 3)-3*t1+std::pow(t1, 3));

  // recursion P[n-1][n-2] and int_P[n-2][n-2] -> int_P[n][n]

  // recursions are instable for small theta
  if(theta1<=10.0*DEG2RAD && theta2<=10.0*DEG2RAD)
  {
    for(UInt n=3; n<=degree; n++)
    {
      Double x1  = 0.0;
      Double x2  = 0.0;
      Double fak = 1.0;

      for(UInt k=0; k<=20; k++)
      {
        for(UInt i=1; i<=k; i++)
          fak *= (2.*i-1.)/(2.*i);

        x1 += fak*std::pow(y1, 2*k)/(n+2+2*k);
        x2 += fak*std::pow(y2, 2*k)/(n+2+2*k);
      }

      intP(n,n) = -factorSmall(n)*(std::pow(y2, n+2)*x2 - std::pow(y1, n+2)*x1);
    }
  }
  else // for large thetas
  {
    for(UInt n=3; n<=degree; n++)
      intP(n,n) = factor1Integral(n,n) * (qy2*P2(n-1,n-2)-qy1*P1(n-1,n-2)) + factor2Integral(n,n) * intP(n-2,n-2);
  }

  // recursion P[n-1][m] and int_P[n-2][m] -> int_P[n][m]
  for(UInt m=0; m<degree; m++)
  {
    intP(m+1,m) = factor1Integral(m+1,m) * (qy2*P2(m,m)-qy1*P1(m,m));
    for(UInt n=m+2; n<=degree; n++)
      intP(n,m) = factor1Integral(n,m) * (qy2*P2(n-1,m)-qy1*P1(n-1,m))
                + factor2Integral(n,m) * intP(n-2,m);
  }

  return intP;
}

/***********************************************/
