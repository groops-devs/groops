/***********************************************/
/**
* @file gravityfield.cpp
*
* @brief Function values of (time variable) gravity fields.
*
* @author Torsten Mayer-Guerr
* @date 2001-08-15
*
*/
/***********************************************/

#define DOCSTRING_Gravityfield

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/configRegister.h"
#include "classes/gravityfield/gravityfieldPotentialCoefficients.h"
#include "classes/gravityfield/gravityfieldPotentialCoefficientsInterior.h"
#include "classes/gravityfield/gravityfieldInInterval.h"
#include "classes/gravityfield/gravityfieldFromParametrization.h"
#include "classes/gravityfield/gravityfieldTimeSplines.h"
#include "classes/gravityfield/gravityfieldTrend.h"
#include "classes/gravityfield/gravityfieldOscillation.h"
#include "classes/gravityfield/gravityfieldTides.h"
#include "classes/gravityfield/gravityfieldTopography.h"
#include "classes/gravityfield/gravityfieldEarthquakeOscillation.h"
#include "classes/gravityfield/gravityfieldFilter.h"
#include "classes/gravityfield/gravityfield.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Gravityfield, "gravityfieldType",
                      GravityfieldPotentialCoefficients,
                      GravityfieldPotentialCoefficientsInterior,
                      GravityfieldFromParametrization,
                      GravityfieldTimeSplines,
                      GravityfieldTrend,
                      GravityfieldOscillation,
                      GravityfieldInInterval,
                      GravityfieldTides,
                      GravityfieldTopography,
                      GravityfieldEarthquakeOscillation,
                      GravityfieldFilter)

GROOPS_READCONFIG_UNBOUNDED_CLASS(Gravityfield, "gravityfieldType")

/***********************************************/

Gravityfield::Gravityfield(Config &config, const std::string &name)
{
  try
  {
    GravityfieldPotentialCoefficients         *coeff         = nullptr;
    GravityfieldPotentialCoefficientsInterior *coeffInterior = nullptr;

    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "functions of the gravity field"))
    {
      renameDeprecatedChoice(config, type, "fromRepresentation", "fromParametrization", date2time(2020, 6, 3));

      if(readConfigChoiceElement(config, "potentialCoefficients", type, "file with potential coefficients"))
      {
        if(!coeff)
          gravityfield.push_back(coeff = new GravityfieldPotentialCoefficients(config));
        else
          coeff->addPotentialCoefficients(GravityfieldPotentialCoefficients(config).sphericalHarmonics(Time()));
      }
      if(readConfigChoiceElement(config, "potentialCoefficientsInterior", type, "file with potential coefficients"))
      {
        if(!coeffInterior)
          gravityfield.push_back(coeffInterior = new GravityfieldPotentialCoefficientsInterior(config));
        else
          coeffInterior->addPotentialCoefficients(GravityfieldPotentialCoefficientsInterior(config).sphericalHarmonics(Time()));
      }
      if(readConfigChoiceElement(config, "fromParametrization",   type, "from a solution vector with given parametrization"))
        gravityfield.push_back(new GravityfieldFromParametrization(config));
      if(readConfigChoiceElement(config, "timeSplines",           type, "file with splines in time domain"))
        gravityfield.push_back(new GravityfieldTimeSplines(config));
      if(readConfigChoiceElement(config, "trend",                 type, "gravityfield as trend"))
        gravityfield.push_back(new GravityfieldTrend(config));
      if(readConfigChoiceElement(config, "oscillation",           type, "gravityfield as oscillation"))
        gravityfield.push_back(new GravityfieldOscillation(config));
      if(readConfigChoiceElement(config, "inInterval",            type, "gravityfields only valid in specific time interval"))
        gravityfield.push_back(new GravityfieldInInterval(config));
      if(readConfigChoiceElement(config, "tides",                 type, "gravityfield from tide models"))
        gravityfield.push_back(new GravityfieldTides(config));
      if(readConfigChoiceElement(config, "topography",            type, "gravityfield from topographic masses"))
        gravityfield.push_back(new GravityfieldTopography(config));
      if(readConfigChoiceElement(config, "earthquakeOscillation", type, "gravityfield from earthquake oscillation"))
        gravityfield.push_back(new GravityfieldEarthquakeOscillation(config));
      if(readConfigChoiceElement(config, "filter",                type, "filtered spherical harmonics"))
        gravityfield.push_back(new GravityfieldFilter(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    };
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Gravityfield::~Gravityfield()
{
  for(UInt i=0; i<gravityfield.size(); i++)
    delete gravityfield.at(i);
}

/***********************************************/

Double Gravityfield::potential(const Time &time, const Vector3d &point) const
{
  Double sum = 0.0;
  for(UInt i=0; i<gravityfield.size(); i++)
    sum += gravityfield.at(i)->potential(time, point);
  return sum;
}

/***********************************************/

Double Gravityfield::radialGradient(const Time &time, const Vector3d &point) const
{
  Double sum = 0.0;
  for(UInt i=0; i<gravityfield.size(); i++)
    sum += gravityfield.at(i)->radialGradient(time, point);
  return sum;
}

/***********************************************/

Double Gravityfield::field(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  Double sum = 0.0;
  for(UInt i=0; i<gravityfield.size(); i++)
    sum += gravityfield.at(i)->field(time, point, kernel);
  return sum;
}

/***********************************************/

Vector3d Gravityfield::gravity(const Time &time, const Vector3d &point) const
{
  Vector3d sum;
  for(UInt i=0; i<gravityfield.size(); i++)
    sum += gravityfield.at(i)->gravity(time, point);
  return sum;
}

/***********************************************/

Tensor3d Gravityfield::gravityGradient(const Time &time, const Vector3d &point) const
{
  Tensor3d sum;
  for(UInt i=0; i<gravityfield.size(); i++)
    sum += gravityfield.at(i)->gravityGradient(time, point);
  return sum;
}

/***********************************************/

Vector3d Gravityfield::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  Vector3d sum;
  for(UInt i=0; i<gravityfield.size(); i++)
    sum += gravityfield.at(i)->deformation(time, point, gravity, hn, ln);
  return sum;
}

/***********************************************/

void Gravityfield::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity, const Vector &hn, const Vector &ln, std::vector< std::vector<Vector3d> > &disp) const
{
  for(UInt i=0; i<gravityfield.size(); i++)
    gravityfield.at(i)->deformation(time, point, gravity, hn, ln, disp);
}

/***********************************************/

SphericalHarmonics Gravityfield::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    if(gravityfield.size()==0)
      return SphericalHarmonics().get(maxDegree, minDegree, GM, R);
    SphericalHarmonics harmonics = gravityfield.at(0)->sphericalHarmonics(time, maxDegree, minDegree, GM, R);
    for(UInt i=1; i<gravityfield.size(); i++)
      harmonics += gravityfield.at(i)->sphericalHarmonics(time, maxDegree, minDegree, GM, R);
    return harmonics;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Gravityfield::sphericalHarmonicsCovariance(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    Matrix Cov;
    for(UInt i=0; i<gravityfield.size(); i++)
    {
      Matrix Cov2 = gravityfield.at(i)->sphericalHarmonicsCovariance(time, maxDegree, minDegree, GM, R);

      if(Cov2.size() == 0)
        continue;

      if(Cov.size() == 0)
      {
        Cov = Cov2;
        continue;
      }

      if(Cov.rows()<Cov2.columns())
      {
        Matrix tmp = Cov;
        Cov = Cov2;
        if((Cov.getType()==Matrix::SYMMETRIC) && (Cov2.getType()==Matrix::SYMMETRIC))
        {
          axpy(1., tmp, Cov.slice(0,0,tmp.rows(),tmp.columns()));
          continue;
        }
        else if((Cov.getType()!=Matrix::SYMMETRIC) && (Cov2.getType()!=Matrix::SYMMETRIC))
        {
          axpy(1., tmp, Cov.row(0,tmp.rows()));
          continue;
        }
        else if((Cov.getType()!=Matrix::SYMMETRIC) && (Cov2.getType()==Matrix::SYMMETRIC))
        {
          for(UInt i=0; i<tmp.rows(); i++)
            Cov(i,i) += tmp(i,0);
          continue;
        }
        else if((Cov.getType()==Matrix::SYMMETRIC) && (Cov2.getType()!=Matrix::SYMMETRIC))
        {
          Cov = Matrix(Cov2.rows(), Matrix::SYMMETRIC);
          copy(tmp, Cov.slice(0,0,tmp.rows(),tmp.columns()));
        }
        else
          throw(Exception("something strange"));
      }

      if((Cov.getType()==Matrix::SYMMETRIC) && (Cov2.getType()==Matrix::SYMMETRIC))
        Cov.slice(0,0,Cov2.rows(),Cov2.columns()) += Cov2;
      else if((Cov.getType()!=Matrix::SYMMETRIC) && (Cov2.getType()!=Matrix::SYMMETRIC))
        Cov.row(0,Cov2.rows()) += Cov2;
      else if((Cov.getType()==Matrix::SYMMETRIC) && (Cov2.getType()!=Matrix::SYMMETRIC))
        for(UInt i=0; i<Cov2.rows(); i++)
          Cov(i,i) += Cov2(i,0);
      else
        throw(Exception("something strange"));
    }

    return Cov;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Gravityfield::variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel) const
{
  Matrix D(point.size(), Matrix::SYMMETRIC);
  for(UInt i=0; i<gravityfield.size(); i++)
    gravityfield.at(i)->variance(time, point, kernel, D);
  return D;
}

/***********************************************/

Double Gravityfield::variance(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  Double sigma2 = 0;
  for(UInt i=0; i<gravityfield.size(); i++)
    sigma2 += gravityfield.at(i)->variance(time, point, kernel);
  return sigma2;
}

/***********************************************/

Double Gravityfield::covariance(const Time &time, const Vector3d &point1, const Vector3d &point2, const Kernel &kernel) const
{
  Double sigma2 = 0;
  for(UInt i=0; i<gravityfield.size(); i++)
    sigma2 += gravityfield.at(i)->covariance(time, point1, point2, kernel);
  return sigma2;
}

/***********************************************/
/***********************************************/

Matrix GravityfieldBase::deformationMatrix(const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                           const Vector &hn, const Vector &ln, Double GM, Double R, UInt maxDegree)
{
  try
  {
    Matrix A(3*point.size(), (maxDegree+1)*(maxDegree+1));

    for(UInt k=0; k<point.size(); k++)
    {
      Vector3d up = normalize(point.at(k));

      Matrix Cnm, Snm;
      SphericalHarmonics::CnmSnm(1/R * point.at(k), maxDegree+1, Cnm, Snm);

      // 0. order
      for(UInt n=0; n<=maxDegree; n++)
      {
        Double wm0 = sqrt((n+1.)*(n+1.));
        Double wp1 = sqrt((n+1.)*(n+2.)) / sqrt(2.0);
        Double Cm0 = wm0*Cnm(n+1,0);
        Double Cp1 = wp1*Cnm(n+1,1); Double Sp1 = wp1*Snm(n+1,1);

        Double   Vn     = GM/R * Cnm(n,0);
        Vector3d gradVn = GM/(2*R) * sqrt((2*n+1.)/(2*n+3.)) * Vector3d(-2*Cp1, -2*Sp1, -2*Cm0);


        Vector3d disp = (hn(n)/gravity.at(k)*Vn) * up // vertical
                      + (ln(n)/gravity.at(k)) * (gradVn-inner(gradVn,up)*up); // horizontal

        A(3*k+0, n*n) = disp.x();
        A(3*k+1, n*n) = disp.y();
        A(3*k+2, n*n) = disp.z();
      }

      // other orders
      for(UInt m=1; m<=maxDegree; m++)
        for(UInt n=m; n<=maxDegree; n++)
        {
          Double wm1 = sqrt((n-m+1.)*(n-m+2.)) * ((m==1) ? sqrt(2.0) : 1.0);
          Double wm0 = sqrt((n-m+1.)*(n+m+1.));
          Double wp1 = sqrt((n+m+1.)*(n+m+2.));
          Double Cm1 = wm1*Cnm(n+1,m-1);  Double Sm1 = wm1*Snm(n+1,m-1);
          Double Cm0 = wm0*Cnm(n+1,m  );  Double Sm0 = wm0*Snm(n+1,m  );
          Double Cp1 = wp1*Cnm(n+1,m+1);  Double Sp1 = wp1*Snm(n+1,m+1);

          Double   Vn     = GM/R * Cnm(n,m);
          Vector3d gradVn = GM/(2*R) * sqrt((2*n+1.)/(2*n+3.)) * Vector3d(Cm1-Cp1, -Sm1-Sp1, -2*Cm0);

          Vector3d disp = (hn(n)/gravity.at(k)*Vn) * up // vertical
                        + (ln(n)/gravity.at(k)) * (gradVn-inner(gradVn,up)*up); // horizontal

          A(3*k+0, n*n+2*m-1) = disp.x();
          A(3*k+1, n*n+2*m-1) = disp.y();
          A(3*k+2, n*n+2*m-1) = disp.z();

          Vn     = GM/R * Snm(n,m);
          gradVn = GM/(2*R) * sqrt((2*n+1.)/(2*n+3.)) * Vector3d(Sm1-Sp1, Cm1+Cp1, -2*Sm0);

          disp = (hn(n)/gravity.at(k)*Vn) * up // vertical
               + (ln(n)/gravity.at(k)) * (gradVn-inner(gradVn,up)*up); // horizontal

          A(3*k+0, n*n+2*m) = disp.x();
          A(3*k+1, n*n+2*m) = disp.y();
          A(3*k+2, n*n+2*m) = disp.z();
        }
    } // for(k=point)

    return A;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Default implementation
Double GravityfieldBase::field(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  try
  {
    return kernel.inverseKernel(time, point, *this);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Default implementation
Matrix GravityfieldBase::sphericalHarmonicsCovariance(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  // only diagonal matrix
  return sphericalHarmonics(time, maxDegree, minDegree, GM, R).sigma2x();
}

/***********************************************/

// Default implementation
Double GravityfieldBase::variance(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  std::vector<Vector3d> p(1);
  p[0] = point;
  Matrix D(p.size(), Matrix::SYMMETRIC);
  variance(time, p, kernel, D);
  return D(0,0);
}

/***********************************************/

// Default implementation
Double GravityfieldBase::covariance(const Time &time, const Vector3d &point1, const Vector3d &point2, const Kernel &kernel) const
{
  std::vector<Vector3d> p(2);
  p[0] = point1;
  p[1] = point2;
  Matrix D(p.size(), Matrix::SYMMETRIC, Matrix::UPPER);
  variance(time, p, kernel, D);
  return D(0,1);
}

/***********************************************/
