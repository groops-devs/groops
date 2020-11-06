/***********************************************/
/**
* @file parametrizationGravitySphericalHarmonics.cpp
*
* @brief Spherical Harmonics.
* @see ParametrizationGravity
*
* @author Torsten Mayer-Guerr
* @date 2003-03-10
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationGravity/parametrizationGravitySphericalHarmonics.h"

/***********************************************/

ParametrizationGravitySphericalHarmonics::ParametrizationGravitySphericalHarmonics(Config &config)
{
  try
  {
    readConfig(config, "minDegree", minDegree, Config::MUSTSET, "2", "");
    readConfig(config, "maxDegree", maxDegree, Config::MUSTSET, "",  "");
    readConfig(config, "GM",        GM,        Config::DEFAULT, STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",         R,         Config::DEFAULT, STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "numbering", numbering, Config::MUSTSET, "",  "numbering scheme");
    if(isCreateSchema(config)) return;

    numbering->numbering(maxDegree, minDegree, idxC, idxS);
    _parameterCount = numbering->parameterCount(maxDegree, minDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravitySphericalHarmonics::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<UInt> degree, order, cs;
    numbering->numbering(maxDegree, minDegree, degree, order, cs);
    for(UInt i=0; i<cs.size(); i++)
      name.push_back(ParameterName("", "sphericalHarmonics."+std::string((cs.at(i)==0) ? "c_" : "s_")+degree.at(i)%"%i"s+"_"+order.at(i)%"%i"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravitySphericalHarmonics::field(const Time &/*time*/, const Vector3d &point, const Kernel &kernel, MatrixSliceRef A) const
{
  try
  {
    Matrix Cnm, Snm;
    SphericalHarmonics::CnmSnm(1/R * point, maxDegree, Cnm, Snm);
    Vector coeff = GM/R * kernel.inverseCoefficients(point, maxDegree);

    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      if(idxC[n][0]!=NULLINDEX) A(0, idxC[n][0]) = coeff(n) * Cnm(n,0);
      for(UInt m=1; m<=n; m++)
      {
        if(idxC[n][m]!=NULLINDEX) A(0, idxC[n][m]) = coeff(n) * Cnm(n,m);
        if(idxS[n][m]!=NULLINDEX) A(0, idxS[n][m]) = coeff(n) * Snm(n,m);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravitySphericalHarmonics::potential(const Time &/*time*/, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix Cnm, Snm;
    SphericalHarmonics::CnmSnm(1/R * point, maxDegree, Cnm, Snm);

    Double factor = GM/R;
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      if(idxC[n][0]!=NULLINDEX) A(0, idxC[n][0]) = factor * Cnm(n,0);
      for(UInt m=1; m<=n; m++)
      {
        if(idxC[n][m]!=NULLINDEX) A(0, idxC[n][m]) = factor * Cnm(n,m);
        if(idxS[n][m]!=NULLINDEX) A(0, idxS[n][m]) = factor * Snm(n,m);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravitySphericalHarmonics::radialGradient(const Time &/*time*/, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix Cnm, Snm;
    SphericalHarmonics::CnmSnm(1/R * point, maxDegree, Cnm, Snm);

    Double factor = -GM/R/point.r();
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      if(idxC[n][0]!=NULLINDEX) A(0, idxC[n][0]) = factor * (n+1) * Cnm(n,0);
      for(UInt m=1; m<=n; m++)
      {
        if(idxC[n][m]!=NULLINDEX) A(0, idxC[n][m]) = factor * (n+1) * Cnm(n,m);
        if(idxS[n][m]!=NULLINDEX) A(0, idxS[n][m]) = factor * (n+1) * Snm(n,m);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravitySphericalHarmonics::gravity(const Time &/*time*/, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix Cnm, Snm;
    SphericalHarmonics::CnmSnm(1/R * point, maxDegree+1, Cnm, Snm);

    // 0. Order
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      Double factor = sqrt((2.*n+1.)/(2.*n+3.))*GM/(2.*R*R);

      Double wm0 = sqrt((n+1.)*(n+1.));
      Double wp1 = sqrt((n+1.)*(n+2.)) / sqrt(2.0);

      Double Cm0 = wm0*Cnm(n+1,0);
      Double Cp1 = wp1*Cnm(n+1,1); Double Sp1 = wp1*Snm(n+1,1);

      if(idxC[n][0]!=NULLINDEX)
      {
        A(0, idxC[n][0]) = factor*(-2*Cp1);
        A(1, idxC[n][0]) = factor*(-2*Sp1);
        A(2, idxC[n][0]) = factor*(-2*Cm0);
      }
    }

    // all other orders
    for(UInt m=1; m<=maxDegree; m++)
      for(UInt n=std::max(minDegree,m); n<=maxDegree; n++)
      {
        Double factor = sqrt((2.*n+1.)/(2.*n+3.))*GM/(2.*R*R);

        Double wm1 = sqrt((n-m+1.)*(n-m+2.)) * ((m==1) ? sqrt(2.0) : 1.0);
        Double wm0 = sqrt((n-m+1.)*(n+m+1.));
        Double wp1 = sqrt((n+m+1.)*(n+m+2.));

        Double Cm1 = wm1*Cnm(n+1,m-1);  Double Sm1 = wm1*Snm(n+1,m-1);
        Double Cm0 = wm0*Cnm(n+1,m  );  Double Sm0 = wm0*Snm(n+1,m  );
        Double Cp1 = wp1*Cnm(n+1,m+1);  Double Sp1 = wp1*Snm(n+1,m+1);

        if(idxC[n][m]!=NULLINDEX)
        {
          A(0, idxC[n][m]) = factor*( Cm1 - Cp1);
          A(1, idxC[n][m]) = factor*(-Sm1 - Sp1);
          A(2, idxC[n][m]) = factor*(-2*Cm0);
        }

        if(idxS[n][m]!=NULLINDEX)
        {
          A(0, idxS[n][m]) = factor*(Sm1 - Sp1);
          A(1, idxS[n][m]) = factor*(Cm1 + Cp1);
          A(2, idxS[n][m]) = factor*(-2*Sm0);
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravitySphericalHarmonics::gravityGradient(const Time &/*time*/, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix Cnm, Snm;
    SphericalHarmonics::CnmSnm(1/R * point, maxDegree+2, Cnm, Snm);

    // 0. Order
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      Double factor = sqrt((2.*n+1.)/(2.*n+5.))*GM/(4.*R*R*R);

      Double wm0 = sqrt((n+1.)*(n+2.)*(n+1.)*(n+2.));
      Double wp1 = sqrt((n+1.)*(n+1.)*(n+2.)*(n+3.)) / sqrt(2.0);
      Double wp2 = sqrt((n+1.)*(n+2.)*(n+3.)*(n+4.)) / sqrt(2.0);

      Double Cm0 = wm0*Cnm(n+2,0);
      Double Cp1 = wp1*Cnm(n+2,1);  Double Sp1 = wp1*Snm(n+2,1);
      Double Cp2 = wp2*Cnm(n+2,2);  Double Sp2 = wp2*Snm(n+2,2);

      if(idxC[n][0]!=NULLINDEX)
      {
        A(0, idxC[n][0]) = factor * (-2*Cm0 + 2*Cp2);
        A(1, idxC[n][0]) = factor * ( 2*Sp2);
        A(2, idxC[n][0]) = factor * ( 4*Cp1);
        A(3, idxC[n][0]) = factor * (-2*Cm0 - 2*Cp2);
        A(4, idxC[n][0]) = factor * ( 4*Sp1);
        A(5, idxC[n][0]) = factor * ( 4*Cm0);
      }
    }

    // 1. order
    UInt m=1;
    for(UInt n=std::max(minDegree, static_cast<UInt>(1)); n<=maxDegree; n++)
    {
      Double factor = sqrt((2.*n+1.)/(2.*n+5.))*GM/(4.*R*R*R);

      Double wm1 = sqrt((n-m+1.)*(n-m+2.)*(n-m+3.)*(n+m+1.)) * sqrt(2.0);
      Double wm0 = sqrt((n-m+1.)*(n-m+2.)*(n+m+1.)*(n+m+2.));
      Double wp1 = sqrt((n-m+1.)*(n+m+1.)*(n+m+2.)*(n+m+3.));
      Double wp2 = sqrt((n+m+1.)*(n+m+2.)*(n+m+3.)*(n+m+4.));

      Double Cm1 = wm1*Cnm(n+2,m-1);  Double Sm1 = wm1*Snm(n+2,m-1);
      Double Cm0 = wm0*Cnm(n+2,m  );  Double Sm0 = wm0*Snm(n+2,m  );
      Double Cp1 = wp1*Cnm(n+2,m+1);  Double Sp1 = wp1*Snm(n+2,m+1);
      Double Cp2 = wp2*Cnm(n+2,m+2);  Double Sp2 = wp2*Snm(n+2,m+2);

      if(idxC[n][m]!=NULLINDEX)
      {
        A(0, idxC[n][m]) = factor * (- 3*Cm0 + Cp2);
        A(1, idxC[n][m]) = factor * (-   Sm0 + Sp2);
        A(2, idxC[n][m]) = factor * (-2*Cm1 + 2*Cp1);
        A(3, idxC[n][m]) = factor * (-   Cm0 - Cp2);
        A(4, idxC[n][m]) = factor * (2*Sp1);
        A(5, idxC[n][m]) = factor * (4*Cm0);
      }

      if(idxS[n][m]!=NULLINDEX)
      {
        A(0, idxS[n][m]) = factor * (- Sm0 + Sp2);
        A(1, idxS[n][m]) = factor * (- Cm0 - Cp2);
        A(2, idxS[n][m]) = factor * (-2*Sm1 + 2*Sp1);
        A(3, idxS[n][m]) = factor * (- 3*Sm0 - Sp2);
        A(4, idxS[n][m]) = factor * (-2*Cm1 - 2*Cp1);
        A(5, idxS[n][m]) = factor * (4*Sm0);
      }
    } // end 1. order

    // all other orders
    for(UInt m=2; m<=maxDegree; m++)
      for(UInt n=std::max(minDegree,m); n<=maxDegree; n++)
      {
        Double factor = sqrt((2.*n+1.)/(2.*n+5.))*GM/(4.*R*R*R);

        Double wm2 = sqrt((n-m+1.)*(n-m+2.)*(n-m+3.)*(n-m+4.)) * ((m==2) ? sqrt(2.0) : 1.0);
        Double wm1 = sqrt((n-m+1.)*(n-m+2.)*(n-m+3.)*(n+m+1.));
        Double wm0 = sqrt((n-m+1.)*(n-m+2.)*(n+m+1.)*(n+m+2.));
        Double wp1 = sqrt((n-m+1.)*(n+m+1.)*(n+m+2.)*(n+m+3.));
        Double wp2 = sqrt((n+m+1.)*(n+m+2.)*(n+m+3.)*(n+m+4.));

        Double Cm2 = wm2*Cnm(n+2,m-2);  Double Sm2 = wm2*Snm(n+2,m-2);
        Double Cm1 = wm1*Cnm(n+2,m-1);  Double Sm1 = wm1*Snm(n+2,m-1);
        Double Cm0 = wm0*Cnm(n+2,m  );  Double Sm0 = wm0*Snm(n+2,m  );
        Double Cp1 = wp1*Cnm(n+2,m+1);  Double Sp1 = wp1*Snm(n+2,m+1);
        Double Cp2 = wp2*Cnm(n+2,m+2);  Double Sp2 = wp2*Snm(n+2,m+2);

        if(idxC[n][m]!=NULLINDEX)
        {
          A(0, idxC[n][m]) = factor * ( Cm2 - 2*Cm0 + Cp2);
          A(1, idxC[n][m]) = factor * (-Sm2         + Sp2);
          A(2, idxC[n][m]) = factor * (-2*Cm1 + 2*Cp1);
          A(3, idxC[n][m]) = factor * (-Cm2 - 2*Cm0 - Cp2);
          A(4, idxC[n][m]) = factor * ( 2*Sm1 + 2*Sp1);
          A(5, idxC[n][m]) = factor * (4*Cm0);
        }

        if(idxS[n][m]!=NULLINDEX)
        {
          A(0, idxS[n][m]) = factor * ( Sm2 - 2*Sm0 + Sp2);
          A(1, idxS[n][m]) = factor * ( Cm2         - Cp2);
          A(2, idxS[n][m]) = factor * (-2*Sm1 + 2*Sp1);
          A(3, idxS[n][m]) = factor * (-Sm2 - 2*Sm0 - Sp2);
          A(4, idxS[n][m]) = factor * (-2*Cm1 - 2*Cp1);
          A(5, idxS[n][m]) = factor * (4*Sm0);
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravitySphericalHarmonics::deformation(const Time &/*time*/, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, MatrixSliceRef A) const
{
  try
  {
    Vector3d up = normalize(point);
    Matrix Cnm, Snm;
    SphericalHarmonics::CnmSnm(1/R * point, maxDegree+1, Cnm, Snm);

    // 0. order
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      Double wm0 = sqrt((n+1.)*(n+1.));
      Double wp1 = sqrt((n+1.)*(n+2.)) / sqrt(2.0);
      Double Cm0 = wm0*Cnm(n+1,0);
      Double Cp1 = wp1*Cnm(n+1,1); Double Sp1 = wp1*Snm(n+1,1);

      if(idxC[n][0]!=NULLINDEX)
      {
        Double   Vn     = GM/R * Cnm(n,0);
        Vector3d gradVn = GM/(2*R) * sqrt((2*n+1.)/(2*n+3.)) * Vector3d(-2*Cp1, -2*Sp1, -2*Cm0);


        Vector3d disp = (hn(n)/gravity*Vn) * up // vertical
                      + (ln(n)/gravity) * (gradVn-inner(gradVn,up)*up); // horizontal

        A(0, idxC[n][0]) = disp.x();
        A(1, idxC[n][0]) = disp.y();
        A(2, idxC[n][0]) = disp.z();
      }
    }

    // other orders
    for(UInt m=1; m<=maxDegree; m++)
      for(UInt n=std::max(minDegree,m); n<=maxDegree; n++)
      {
        Double wm1 = sqrt((n-m+1.)*(n-m+2.)) * ((m==1) ? sqrt(2.0) : 1.0);
        Double wm0 = sqrt((n-m+1.)*(n+m+1.));
        Double wp1 = sqrt((n+m+1.)*(n+m+2.));
        Double Cm1 = wm1*Cnm(n+1,m-1);  Double Sm1 = wm1*Snm(n+1,m-1);
        Double Cm0 = wm0*Cnm(n+1,m  );  Double Sm0 = wm0*Snm(n+1,m  );
        Double Cp1 = wp1*Cnm(n+1,m+1);  Double Sp1 = wp1*Snm(n+1,m+1);

        if(idxC[n][m]!=NULLINDEX)
        {
          Double   Vn     = GM/R * Cnm(n,m);
          Vector3d gradVn = GM/(2*R) * sqrt((2*n+1.)/(2*n+3.)) * Vector3d(Cm1-Cp1, -Sm1-Sp1, -2*Cm0);

          Vector3d disp = (hn(n)/gravity*Vn) * up // vertical
                        + (ln(n)/gravity) * (gradVn-inner(gradVn,up)*up); // horizontal

          A(0, idxC[n][m]) = disp.x();
          A(1, idxC[n][m]) = disp.y();
          A(2, idxC[n][m]) = disp.z();
        }

        if(idxS[n][m]!=NULLINDEX)
        {
          Double   Vn     = GM/R * Snm(n,m);
          Vector3d gradVn = GM/(2*R) * sqrt((2*n+1.)/(2*n+3.)) * Vector3d(Sm1-Sp1, Cm1+Cp1, -2*Sm0);

          Vector3d disp = (hn(n)/gravity*Vn) * up // vertical
                        + (ln(n)/gravity) * (gradVn-inner(gradVn,up)*up); // horizontal

          A(0, idxS[n][m]) = disp.x();
          A(1, idxS[n][m]) = disp.y();
          A(2, idxS[n][m]) = disp.z();
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

SphericalHarmonics ParametrizationGravitySphericalHarmonics::sphericalHarmonics(const Time &/*time*/, const Vector &x, UInt _maxDegree) const
{
  try
  {
    Matrix cnm(this->maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(this->maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      if(idxC[n][0]!=NULLINDEX) cnm(n,0) = x(idxC[n][0]);
      for(UInt m=1; m<=n; m++)
      {
        if(idxC[n][m]!=NULLINDEX) cnm(n,m) = x(idxC[n][m]);
        if(idxS[n][m]!=NULLINDEX) snm(n,m) = x(idxS[n][m]);
      }
    }

    if(_maxDegree==INFINITYDEGREE)
      _maxDegree = this->maxDegree;
    return SphericalHarmonics(GM,R, cnm,snm).get(maxDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics ParametrizationGravitySphericalHarmonics::sphericalHarmonics(const Time &/*time*/, const Vector &x, const Vector &sigma2x, UInt _maxDegree) const
{
  try
  {
    Matrix cnm      (this->maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm      (this->maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2cnm(this->maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2snm(this->maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    for(UInt n=minDegree; n<=this->maxDegree; n++)
    {
      if(idxC[n][0]!=NULLINDEX) cnm(n,0)       = x(idxC[n][0]);
      if(idxC[n][0]!=NULLINDEX) sigma2cnm(n,0) = sigma2x(idxC[n][0]);
      for(UInt m=1; m<=n; m++)
      {
        if(idxC[n][0]!=NULLINDEX) cnm(n,m)       = x(idxC[n][m]);
        if(idxS[n][m]!=NULLINDEX) snm(n,m)       = x(idxS[n][m]);
        if(idxC[n][0]!=NULLINDEX) sigma2cnm(n,m) = sigma2x(idxC[n][m]);
        if(idxS[n][m]!=NULLINDEX) sigma2snm(n,m) = sigma2x(idxS[n][m]);
      }
    }

    if(_maxDegree==INFINITYDEGREE)
      _maxDegree = this->maxDegree;
    return SphericalHarmonics(GM,R, cnm,snm, sigma2cnm, sigma2snm).get(_maxDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
