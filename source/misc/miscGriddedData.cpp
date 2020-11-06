/***********************************************/
/**
* @file miscGriddedData.cpp
*
* @brief Misc functions for values on grid.
*
* @author Torsten Mayer-Guerr
* @author Christian Pock
* @date 2008-08-06
*
*/
/***********************************************/

#include "base/import.h"
#include "inputOutput/logging.h"
#include "miscGriddedData.h"

/***********************************************/

namespace MiscGriddedData
{

/***********************************************/

void statistics(const std::vector<Double> &values, const std::vector<Double> &weights,
                Double &rms, Double &avg, Double &vmin, Double &vmax, Double &mean)
{
  try
  {
    mean = 0.;
    avg  = 0.;
    rms  = 0.;
    vmin = NAN_EXPR;
    vmax = NAN_EXPR;
    if(!values.size())
      return;

    const Double sumWeight = std::accumulate(weights.begin(), weights.end(), Double(0.));
    for(UInt i=0; i<values.size(); i++)
    {
      const Double w = (sumWeight) ? (weights.at(i)/sumWeight) : (1./values.size());
      const Double v = values.at(i);
      mean += w * v;
      avg  += w * std::fabs(v);
      rms  += w * v*v;
    }
    rms  = std::sqrt(rms);
    vmin = *std::min_element(values.begin(), values.end());
    vmax = *std::max_element(values.begin(), values.end());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static void statistics(const std::vector<std::vector<Double>> &values, const std::vector<Double> &weights,
                       Vector &rms, Vector &avg, Vector &vmin, Vector &vmax, Vector &mean)
{
  try
  {
    rms = avg = vmin = vmax = mean = Vector(values.size());
    for(UInt i=0; i<values.size(); i++)
      statistics(values.at(i), weights, rms(i), avg(i), vmin(i), vmax(i), mean(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static void printStatistics(const Vector &rms, const Vector &avg, const Vector &vmin, const Vector &vmax, const Vector &mean)
{
  try
  {
    if(rms.rows() == 0)
      return;

    const UInt width = 13;
    auto line = [&](const std::string &str, const Vector &x)
    {
      std::stringstream ss;
      for(UInt i=0; i<rms.rows(); i++)
        ss<<std::setw(width)<<std::left<<x(i);
      logInfo<<"  "<<str<<ss.str()<<Log::endl;
    };

    logInfo<<"data statistics"<<Log::endl;
    if(rms.rows()>1)
    {
      std::stringstream ss;
      for(UInt i=0; i<rms.rows(); i++)
        ss<<"data"<<std::setw(width-4)<<std::left<<(i);
      logInfo<<"        "<<ss.str()<<Log::endl;
    }

    line("rms:  ", rms);
    line("avg:  ", avg);
    line("min:  ", vmin);
    line("max:  ", vmax);
    line("mean: ", mean);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void printStatistics(const GriddedData &grid)
{
  try
  {
    if(!Parallel::isMaster()) return;

    const Double totalArea = std::accumulate(grid.areas.begin(), grid.areas.end(), Double(0.));

    logInfo<<"grid statistics"<<Log::endl;
    std::vector<Angle>  lambda, phi;
    std::vector<Double> radius;
    if(grid.isRectangle(lambda, phi, radius))
    {
      Angle  lon1, lon2, lat1, lat2;
      Double h1, h2;
      grid.ellipsoid(grid.points.front(), lon1, lat1, h1);
      grid.ellipsoid(grid.points.back(),  lon2, lat2, h2);
      logInfo<<"  regular grid ("<<phi.size()<<" x "<<lambda.size()<<") = "<<grid.points.size()<<Log::endl;
      logInfo<<"  longitude: "<<lon1*RAD2DEG<<"° -- "<<lon2*RAD2DEG<<"°"<<Log::endl;
      logInfo<<"  latitude:  "<<lat1*RAD2DEG<<"° -- "<<lat2*RAD2DEG<<"°"<<Log::endl;
      logInfo<<"  area:      "<<totalArea/(4*PI)*100<<"% of Earth's surface ("<<totalArea*std::pow(DEFAULT_R/1000,2)<<" km^2)"<<Log::endl;
    }
    else
    {
      logInfo<<"  count: "<<grid.points.size()<<Log::endl;
      logInfo<<"  area:  "<<totalArea/(4*PI)*100<<"% of Earth's surface ("<<totalArea*pow(DEFAULT_R/1000,2)<<" km^2)"<<Log::endl;
    }

    if(grid.values.size() && grid.points.size())
    {
      Vector rms, avg, vmin, vmax, mean;
      statistics(grid.values, grid.areas, rms, avg, vmin, vmax, mean);
      printStatistics(rms, avg, vmin, vmax, mean);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void printStatistics(const GriddedDataRectangular &grid)
{
  try
  {
    if(!Parallel::isMaster()) return;

    if(!grid.isValid())
    {
      logInfo<<"grid is not valid"<<Log::endl;
      return;
    }

    const UInt rows = grid.latitudes.size();
    const UInt cols = grid.longitudes.size();

    std::vector<Double> radius, dLambda, dPhi;
    std::vector<Angle>  lambda, phi;
    grid.geocentric(lambda, phi, radius, dLambda, dPhi);

    Double totalArea = 0.0;
    Double dL = 0;
    for(UInt s=0; s<cols; s++)
      dL += dLambda[s];
    for(UInt z=0; z<rows; z++)
      totalArea += std::cos(phi.at(z))*dPhi.at(z)*dL;

    logInfo<<"grid statistics"<<Log::endl;
    logInfo<<"  regular grid ("<<rows<<" x "<<cols<<") = "<<rows*cols<<Log::endl;
    logInfo<<"  longitude: "<<grid.longitudes.at(0)*RAD2DEG<<"° -- "<<grid.longitudes.back()*RAD2DEG<<"°"<<Log::endl;
    logInfo<<"  latitude:  "<<grid.latitudes.at(0) *RAD2DEG<<"° -- "<<grid.latitudes.back() *RAD2DEG<<"°"<<Log::endl;
    logInfo<<"  area:      "<<totalArea/(4*PI)*100<<"% of Earth's surface ("<<totalArea*std::pow(DEFAULT_R/1000,2)<<" km^2)"<<Log::endl;

    if(!grid.values.size())
      return;

    Vector mean(grid.values.size());
    Vector vmin(grid.values.size()), vmax(grid.values.size());
    Vector avg(grid.values.size()),  rms(grid.values.size());
    for(UInt i=0; i<grid.values.size(); i++)
    {
      for(UInt z=0; z<rows; z++)
      {
        const Double weight = dPhi.at(z)*std::cos(phi.at(z))/totalArea;
        for(UInt s=0; s<cols; s++)
        {
          const Double w = dLambda[s]*weight;
          const Double v = grid.values[i](z,s);
          mean(i) += w * v;
          avg(i)  += w * std::fabs(v);
          rms(i)  += w * v * v;
        }
      }
      vmin(i) = min(grid.values[i]);
      vmax(i) = max(grid.values[i]);
      rms(i)  = std::sqrt(rms(i));
    }

    printStatistics(rms, avg, vmin, vmax, mean);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::vector<Double> synthesisSphericalHarmonics(const SphericalHarmonics &harm, const std::vector<Vector3d> &points, KernelPtr kernel, Parallel::CommunicatorPtr comm, Bool timing)
{
  try
  {
    std::vector<Double>   field(points.size(), 0.);
    std::vector<Angle>    lambda, phi;
    std::vector<Double>   r;
    if(GriddedData(Ellipsoid(), points, std::vector<Double>(), std::vector<std::vector<Double>>()).isRectangle(lambda, phi, r))
    {
      if(!Parallel::isMaster(comm))
        return field;

      // spherical harmonics with recatangular grid
      Matrix cossinm(lambda.size(), 2*harm.maxDegree()+1);
      for(UInt i=0; i<lambda.size(); i++)
      {
        cossinm(i,0) = 1.;
        for(UInt m=1; m<=harm.maxDegree(); m++)
        {
          cossinm(i,2*m-1) = cos(m*static_cast<Double>(lambda.at(i)));
          cossinm(i,2*m+0) = sin(m*static_cast<Double>(lambda.at(i)));
        }
      }

      // Compute Legendre functions for each phi (row)
      // Parallel::forEach(phi.size(), [&](UInt i)
      for(UInt i=0; i<phi.size(); i++)
      {
        const Vector3d p  = polar(lambda.at(0), phi.at(i), r.at(i));
        const Vector   kn = kernel->inverseCoefficients(p, harm.maxDegree(), harm.isInterior());

        Matrix Pnm = SphericalHarmonics::Pnm(Angle(PI/2-phi.at(i)), r.at(i)/harm.R(), harm.maxDegree(), harm.isInterior());
        for(UInt n=0; n<=harm.maxDegree(); n++)
          Pnm.slice(n,0,1,n+1) *= harm.GM()/harm.R()*kn(n);

        Vector sum(2*harm.maxDegree()+1);
        sum(0) = inner(harm.cnm().column(0), Pnm.column(0));
        for(UInt m=1; m<=harm.maxDegree(); m++)
        {
          sum(2*m-1) = inner(harm.cnm().slice(m,m,harm.maxDegree()-m+1,1), Pnm.slice(m,m,harm.maxDegree()-m+1,1));
          sum(2*m+0) = inner(harm.snm().slice(m,m,harm.maxDegree()-m+1,1), Pnm.slice(m,m,harm.maxDegree()-m+1,1));
        }
        Vector row = cossinm * sum;
        for(UInt k=0; k<lambda.size(); k++)
          field.at(i*lambda.size()+k) = row(k);
      } //);
      // Parallel::reduceSum(field);
      return field;
    } // if(isRectangle)

    // spherical harmonics with arbitrary point distribution
    Parallel::forEach(field, [&](UInt i)
                      {return inner(kernel->inverseCoefficients(points.at(i), harm.maxDegree(), harm.isInterior()),
                                    harm.Yn(points.at(i), harm.maxDegree()));}, comm, timing);
    return field;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix synthesisSphericalHarmonicsMatrix(UInt maxDegree, Double GM, Double R, const std::vector<Vector3d> &points, KernelPtr kernel, Bool isInterior)
{
  try
  {
    Matrix A = Matrix(points.size(), (maxDegree+1)*(maxDegree+1));
    for(UInt k=0; k<points.size(); k++)
    {
      Matrix Cnm, Snm;
      SphericalHarmonics::CnmSnm(1./R * points.at(k), maxDegree, Cnm, Snm, isInterior);
      Vector kn = kernel->inverseCoefficients(points.at(k), maxDegree, isInterior);
      UInt idx = 0;
      for(UInt n=0; n<=maxDegree; n++)
      {
        A(k,idx++) =  kn(n) * GM/R * Cnm(n,0);
        for(UInt m=1; m<=n; m++)
        {
          A(k,idx++) = kn(n) * GM/R * Cnm(n,m);
          A(k,idx++) = kn(n) * GM/R * Snm(n,m);
        }
      }
    }

    return A;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

} // end namespace GriddedData
