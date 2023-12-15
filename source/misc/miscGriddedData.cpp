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
    if(!grid.isValid())
    {
      logInfo<<"grid is not valid"<<Log::endl;
      return;
    }

    std::vector<Double> dLambda, dPhi;
    const Double totalArea = grid.areaElements(dLambda, dPhi);
    const UInt   rows      = grid.latitudes.size();
    const UInt   cols      = grid.longitudes.size();

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
    for(UInt idx=0; idx<grid.values.size(); idx++)
    {
      for(UInt i=0; i<rows; i++)
        for(UInt k=0; k<cols; k++)
        {
          const Double w = dLambda[k]*dPhi[i]/totalArea;
          const Double v = grid.values[idx](i,k);
          mean(idx) += w * v;
          avg(idx)  += w * std::fabs(v);
          rms(idx)  += w * v * v;
        }
      vmin(idx) = min(grid.values[idx]);
      vmax(idx) = max(grid.values[idx]);
      rms(idx)  = std::sqrt(rms(idx));
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

    // spherical harmonics with recatangular grid
    if(GriddedData(Ellipsoid(), points, std::vector<Double>(), std::vector<std::vector<Double>>()).isRectangle(lambda, phi, r))
    {
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
      Parallel::forEach(phi.size(), [&](UInt i)
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
      }, comm, timing);
      Parallel::reduceSum(field, 0, comm);
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
/***********************************************/

std::vector<SphericalHarmonics> analysisSphericalHarmonics(const GriddedData &grid, KernelPtr kernel, UInt minDegree, UInt maxDegree, Double GM, Double R,
                                                           Bool useLeastSquares, Parallel::CommunicatorPtr comm, Bool timing)
{
  try
  {
    // test rectangular grid
    // ---------------------
    std::vector<Angle>  lambda, phi;
    std::vector<Double> radius;
    const Bool isRectangle = grid.isRectangle(lambda, phi, radius);
    // precompute cos(m*lambda), sin(m*lambda)
    Matrix cosml, sinml;
    if(isRectangle)
    {
      cosml = sinml = Matrix(lambda.size(), maxDegree+1);
      for(UInt m=0; m<=maxDegree; m++)
        for(UInt k=0; k<lambda.size(); k++)
        {
          cosml(k, m) = std::cos(m*static_cast<Double>(lambda.at(k)));
          sinml(k, m) = std::sin(m*static_cast<Double>(lambda.at(k)));
        }
    }

    // quadrature formular
    // -------------------
    if(!useLeastSquares)
    {
      if(timing) logStatus<<"computing quadrature formular"<<Log::endl;
      std::vector<Matrix> cnm(grid.values.size(), Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
      std::vector<Matrix> snm(grid.values.size(), Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));

      if(isRectangle)
      {
        Parallel::forEach(phi.size(), [&](UInt i)
        {
          // legendre functions with kernel coefficients
          Matrix Pnm = SphericalHarmonics::Pnm(Angle(0.5*PI-phi.at(i)), radius.at(i)/R, maxDegree, TRUE);
          Vector kn  = kernel->coefficients(polar(Angle(0.), phi.at(i), radius.at(i)), maxDegree);
          for(UInt n=0; n<=maxDegree; n++)
            Pnm.slice(n, 0, 1, n+1) *= kn(n) * R/(4*PI*GM);

          for(UInt idx=0; idx<grid.values.size(); idx++)
            for(UInt k=0; k<lambda.size(); k++)
            {
              const Double f = grid.values.at(idx).at(i*lambda.size()+k) * grid.areas.at(i*lambda.size()+k);
              for(UInt m=0; m<=maxDegree; m++)
              {
                axpy(f * cosml(k, m), Pnm.column(m), cnm.at(idx).column(m));
                axpy(f * sinml(k, m), Pnm.column(m), snm.at(idx).column(m));
              }
            }
        }, comm, timing);
      }
      else
      {
        Parallel::forEach(grid.points.size(), [&](UInt i)
        {
          const Vector kn = kernel->coefficients(grid.points.at(i), maxDegree);
          Matrix Cnm, Snm;
          SphericalHarmonics::CnmSnm(normalize(grid.points.at(i)), maxDegree, Cnm, Snm);
          for(UInt idx=0; idx<grid.values.size(); idx++)
          {
            Double rR = std::pow(grid.points.at(i).r()/R, minDegree+1);
            for(UInt n=minDegree; n<=maxDegree; n++)
            {
              axpy((kn(n)* R/(4*PI*GM) * rR * grid.values.at(idx).at(i) * grid.areas.at(i)), Cnm.row(n), cnm.at(idx).row(n));
              axpy((kn(n)* R/(4*PI*GM) * rR * grid.values.at(idx).at(i) * grid.areas.at(i)), Snm.row(n), snm.at(idx).row(n));
              rR *= grid.points.at(i).r()/R;
            }
          }
        }, comm, timing);
      }

      std::vector<SphericalHarmonics> harm(grid.values.size());
      for(UInt idx=0; idx<grid.values.size(); idx++)
      {
        Parallel::reduceSum(cnm.at(idx), 0, comm);
        Parallel::reduceSum(snm.at(idx), 0, comm);
        harm.at(idx) = SphericalHarmonics(GM, R, cnm.at(idx), snm.at(idx)).get(maxDegree, minDegree);
      }

      return harm;
    }

    // least squares adjustment order by order
    // ---------------------------------------
    if(timing) logStatus<<"least squares adjustment (order by order)"<<Log::endl;
    if(!isRectangle)
      throw(Exception("GriddedData must be a rectangle grid"));

    // system of normal equations (order by order)
    std::vector<Matrix> N, n;
    N.push_back(Matrix(maxDegree+1, Matrix::SYMMETRIC));
    n.push_back(Matrix(maxDegree+1, grid.values.size()));
    for(UInt m=1; m<=maxDegree; m++)
    {
      N.push_back(Matrix(2*(maxDegree+1-m), Matrix::SYMMETRIC));
      n.push_back(Matrix(2*(maxDegree+1-m), grid.values.size()));
    }

    Vector lPl(grid.values.size());
    for(UInt idx=0; idx<grid.values.size(); idx++)
      for(UInt i=0; i<grid.points.size(); i++)
        lPl(idx) += grid.values.at(idx).at(i) * grid.areas.at(i)/(4*PI) * grid.values.at(idx).at(i);

    if(timing) logStatus<<"accumulate normal equations"<<Log::endl;
    Parallel::forEach(phi.size(), [&](UInt i)
    {
      Matrix l(lambda.size(), grid.values.size());
      for(UInt idx=0; idx<grid.values.size(); idx++)
        for(UInt k=0; k<lambda.size(); k++)
          l(k, idx) = grid.values.at(idx).at(i*lambda.size()+k);

      // legendre functions with kernel coefficients
      Matrix Pnm = SphericalHarmonics::Pnm(Angle(0.5*PI-phi.at(i)), radius.at(i)/R, maxDegree);
      Vector kn  = kernel->inverseCoefficients(polar(Angle(0.), phi.at(i), radius.at(i)), maxDegree);
      for(UInt n=0; n<=maxDegree; n++)
        Pnm.slice(n, 0, 1, n+1) *= GM/R * kn(n);

      const Double weight = grid.areas.at(i*lambda.size())/(4*PI); // assume same area for all longitudes
      Matrix A = cosml.column(0) * Pnm.column(0).trans();
      rankKUpdate(weight, A, N.at(0));
      matMult(weight, A.trans(), l, n.at(0));
      for(UInt m=1; m<=maxDegree; m++)
      {
        Matrix A(lambda.size(), 2*(maxDegree+1-m));
        matMult(1.0, cosml.column(m), Pnm.slice(m, m, maxDegree+1-m, 1).trans(), A.column(0, maxDegree+1-m));
        matMult(1.0, sinml.column(m), Pnm.slice(m, m, maxDegree+1-m, 1).trans(), A.column(maxDegree+1-m, maxDegree+1-m));
        rankKUpdate(weight, A, N.at(m));
        matMult(weight, A.trans(), l, n.at(m));
      }
    }, comm, timing);

    for(UInt m=0; m<=maxDegree; m++)
    {
      Parallel::reduceSum(N.at(m), 0, comm);
      Parallel::reduceSum(n.at(m), 0, comm);
    }

    // solve normals
    // -------------
    std::vector<SphericalHarmonics> harm(grid.values.size());
    if(Parallel::isMaster(comm))
    {
      if(timing) logStatus<<"solve the system of equations"<<Log::endl;
      std::vector<Matrix> x(maxDegree+1);
      std::vector<Vector> sigma2x(maxDegree+1);
      UInt parameterCount = 0;
      for(UInt m=0; m<=maxDegree; m++)
      {
        parameterCount += n.at(m).rows();
        for(UInt k=0; k<N.at(m).rows(); k++)
          if(N.at(m)(k, k) == 0.)
          {
            N.at(m)(k, k) = 1;
            parameterCount--;
            // logWarning<<k<<". parameter has zero diagonal element -> set to one"<<Log::endl;
          }
        x.at(m) = solve(N.at(m), n.at(m));
        inverse(N.at(m));               // inverse of the cholesky matrix
        sigma2x.at(m) = Vector(x.at(m).rows());
        for(UInt k=0; k<N.at(m).rows(); k++)
          sigma2x.at(m)(k) = quadsum(N.at(m).slice(k, k, 1, N.at(m).columns()-k));
      }

      // aposteriori sigma
      Vector sigma2(grid.values.size());
      for(UInt idx=0; idx<grid.values.size(); idx++)
      {
        Double ePe = lPl(idx);
        for(UInt m=0; m<=maxDegree; m++)
          ePe -= inner(x.at(m).column(idx), n.at(m).column(idx));
        sigma2(idx) = std::max(ePe/(grid.points.size()-parameterCount), 0.);
      }

      // potential coefficients
      Matrix cnm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix snm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix sigma2cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix sigma2snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      for(UInt idx=0; idx<grid.values.size(); idx++)
      {
        copy(x.at(0).column(idx), cnm.column(0));
        copy(sigma2x.at(0), sigma2cnm.column(0));
        for(UInt m=1; m<=maxDegree; m++)
        {
          copy(x.at(m).slice(0,             idx, maxDegree+1-m, 1),    cnm.slice(m, m, maxDegree+1-m, 1));
          copy(x.at(m).slice(maxDegree+1-m, idx, maxDegree+1-m, 1),    snm.slice(m, m, maxDegree+1-m, 1));
          copy(sigma2x.at(m).row(0,              maxDegree+1-m), sigma2cnm.slice(m, m, maxDegree+1-m, 1));
          copy(sigma2x.at(m).row(maxDegree+1-m,  maxDegree+1-m), sigma2snm.slice(m, m, maxDegree+1-m, 1));
        }
        harm.at(idx) = SphericalHarmonics(GM, R, cnm, snm, sigma2(idx)*sigma2cnm, sigma2(idx)*sigma2snm).get(maxDegree, minDegree);
      }
    } // if(Parallel::isMaster(comm))

    return harm;
 }
 catch(std::exception &e)
 {
   GROOPS_RETHROW(e)
 }
}

/***********************************************/

} // end namespace GriddedData
