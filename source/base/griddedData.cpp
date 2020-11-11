/***********************************************/
/**
* @file griddedData.cpp
*
* @brief Gridded values.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#include "base/importStd.h"
#include "griddedData.h"

/***********************************************/

void GriddedData::sort()
{
  try
  {
    std::vector<UInt> index(points.size());
    std::iota(index.begin(), index.end(), 0);
    std::stable_sort(index.begin(), index.end(), [&](UInt i, UInt k)
                     {return (std::fabs(points.at(i).theta()-points.at(k).theta()) < 1e-10) ? // same latitude?
                             (points.at(i).lambda() < points.at(k).lambda()) :                // sort longitudes
                             (points.at(i).theta()  < points.at(k).theta());});               // sort latitudes

    auto tmpPoints = points;
    for(UInt i=0; i<tmpPoints.size(); i++)
      points.at(i) = tmpPoints.at(index.at(i));

    auto tmpAreas = areas;
    for(UInt i=0; i<tmpAreas.size(); i++)
      areas.at(i) = tmpAreas.at(index.at(i));

    for(UInt k=0; k<values.size(); k++)
    {
      auto tmpValues = values.at(k);
      for(UInt i=0; i<tmpValues.size(); i++)
        values.at(k).at(i) = tmpValues.at(index.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GriddedData::isRectangle(std::vector<Angle> &lambda, std::vector<Angle> &phi, std::vector<Double> &radius) const
{
  try
  {
    if(!points.size())
      return FALSE;

    const Angle  phi1 = points.at(0).phi();
    const Double r1   = points.at(0).r();
    const Double eps  = 1e-10 * r1;

    // first row
    lambda.clear();
    for(UInt i=0; (i<points.size()) && (std::fabs(points.at(i).r()-r1) < eps) && (std::fabs(points.at(i).phi()-phi1) < 1e-10); i++)
      lambda.push_back(points.at(i).lambda());

    if(points.size() % lambda.size())
      return FALSE;

    // check points
    phi.clear();
    radius.clear();
    for(UInt i=0; i*lambda.size()<points.size(); i++)
    {
      phi.push_back(points.at(i*lambda.size()).phi());
      radius.push_back(points.at(i*lambda.size()).r());
      for(UInt k=0; k<lambda.size(); k++)
        if((points.at(i*lambda.size()+k)-polar(lambda.at(k), phi.at(i), radius.at(i))).r() > eps)
          return FALSE;
    }

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GriddedData::computeArea()
{
  try
  {
    std::vector<Angle>  lambda, phi;
    std::vector<Double> radius;
    if(!isRectangle(lambda, phi, radius)) // || !std::is_sorted(phi.begin(), phi.end()) || !std::is_sorted(lambda.begin(), lambda.end()))
      return FALSE;

    Vector dx(lambda.size(), 2*PI);
    if(lambda.size() > 1)
    {
      dx(0) = std::fabs(lambda.at(1)-lambda.at(0));
      for(UInt i=1; i<lambda.size()-1; i++)
        dx(i) = 0.5*std::fabs(lambda.at(i+1)-lambda.at(i-1));
      dx(dx.rows()-1) = std::fabs(lambda.at(lambda.size()-1)-lambda.at(lambda.size()-2));
    }

    Vector dy(phi.size(), 2.);
    if(phi.size() > 1)
    {
      dy(0) = std::fabs(std::cos(phi.at(0)) * 2*std::sin(0.5*(phi.at(0)-phi.at(1))));
      for(UInt i=1; i<phi.size()-1; i++)
        dy(i) = std::fabs(std::cos(phi.at(i)) * (std::sin(0.5*(phi.at(i)-phi.at(i+1)))+std::sin(0.5*(phi.at(i-1)-phi.at(i)))));
      dy(dy.rows()-1) = std::fabs(std::cos(phi.at(phi.size()-1)) * 2*std::sin(0.5*(phi.at(phi.size()-2)-phi.at(phi.size()-1))));
    }

    areas.resize(points.size());
    for(UInt z=0; z<phi.size(); z++)
      for(UInt s=0; s<lambda.size(); s++)
        areas.at(z*lambda.size()+s) = dx(s)*dy(z);

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GriddedData::isValid() const
{
  try
  {
    const UInt count = points.size();

    if((areas.size()!=0) && (areas.size()!=count))
      return FALSE;

    for(UInt i=0; i<values.size(); i++)
      if(values.at(i).size() != count)
        return FALSE;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void GriddedDataRectangular::geocentric(std::vector<Angle> &lambda, std::vector<Angle> &phi, std::vector<Double> &radius,
                                        std::vector<Double> &dLambda, std::vector<Double> &dPhi) const
{
  try
  {
    lambda = longitudes;
    phi.resize(latitudes.size());
    radius.resize(heights.size());
    for(UInt i=0; i<phi.size(); i++)
    {
      const Vector3d points = ellipsoid(Angle(0), latitudes.at(i), heights.at(i));
      phi.at(i)      = points.phi();
      radius.at(i)   = points.r();
    }

    // areas elements
    // -------------
    const UInt cols = lambda.size();
    dLambda.clear();
    dLambda.resize(cols, 2*PI);
    if(cols > 1)
    {
      dLambda.at(0) = std::fabs(lambda.at(1)-lambda.at(0));
      for(UInt s=1; s<cols-1; s++)
        dLambda.at(s) = std::fabs(0.5*(lambda.at(s+1)-lambda.at(s-1)));
      dLambda.at(cols-1) = std::fabs(lambda.at(cols-1)-lambda.at(cols-2));
    }

    // \int_{B0-dB/2}^{B0+dB/2} cosB dB = cosB0 * 2*sin(dB/2)
    const UInt rows = phi.size();
    dPhi.clear();
    dPhi.resize(rows, 2.);
    if(rows > 1)
    {
      dPhi.at(0) = std::fabs(2*std::sin((phi.at(0)-phi.at(1))/2));
      for(UInt i=1; i<rows-1; i++)
        dPhi.at(i) = std::fabs(2*std::sin((phi.at(i-1)-phi.at(i+1))/4));
      dPhi.at(rows-1) = std::fabs(2*std::sin((phi.at(rows-2)-phi.at(rows-1))/2));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GriddedDataRectangular::convert(GriddedData &grid) const
{
  try
  {
    if(!isValid())
      throw(Exception("GriddedDataRectangular is not valid"));

    std::vector<Double> radius, dLambda, dPhi;
    std::vector<Angle>  lambda;
    std::vector<Angle>  phi;
    geocentric(lambda, phi, radius, dLambda, dPhi);
    const UInt rows = phi.size();
    const UInt cols = lambda.size();

    std::vector<Double> cosL(cols), sinL(cols);
    for(UInt s=0; s<cols; s++)
    {
      cosL[s] = std::cos(lambda[s]);
      sinL[s] = std::sin(lambda[s]);
    }

    grid.ellipsoid = ellipsoid;
    grid.points.resize(rows*cols);
    grid.areas.resize(rows*cols);
    for(UInt z=0; z<rows; z++)
    {
      const Double cosB = std::cos(phi[z]);
      const Double sinB = std::sin(phi[z]);
      for(UInt s=0; s<cols; s++)
      {
        grid.points[z*cols+s] = Vector3d(radius[z]*cosB*cosL[s], radius[z]*cosB*sinL[s], radius[z]*sinB);
        grid.areas[z*cols+s]  = dLambda[s]*dPhi[z]*cosB;
      }
    }

    // values
    grid.values.resize(values.size());
    for(UInt i=0; i<grid.values.size(); i++)
    {
      grid.values.at(i).resize(rows*cols);
      for(UInt z=0; z<rows; z++)
        for(UInt s=0; s<cols; s++)
          grid.values[i][z*cols+s] = values[i](z,s);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GriddedDataRectangular::init(const GriddedData &grid)
{
  try
  {
    if(!grid.isValid())
      return FALSE;

    std::vector<Double> radius;
    std::vector<Angle>  lambda;
    std::vector<Angle>  phi;
    if(!grid.isRectangle(lambda, phi, radius))
      return FALSE;

    const UInt rows = phi.size();
    const UInt cols = lambda.size();
    ellipsoid = grid.ellipsoid;

    longitudes = lambda;
    latitudes.resize(rows);
    heights.resize(rows);
    Angle L0;
    for(UInt z=0; z<rows; z++)
      ellipsoid(polar(Angle(0), phi.at(z), radius.at(z)), L0, latitudes.at(z), heights.at(z));

    values.resize( grid.values.size() );
    for(UInt i=0; i<values.size(); i++)
    {
      values.at(i) = Matrix(rows, cols);
      for(UInt z=0; z<rows; z++)
        for(UInt s=0; s<cols; s++)
          values.at(i)(z,s) = grid.values.at(i).at(z*cols+s);
    }

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GriddedDataRectangular::isValid() const
{
  try
  {
    const UInt rows = latitudes.size();
    const UInt cols = longitudes.size();

    if(heights.size() != rows)
      return FALSE;

    for(UInt i=0; i<values.size(); i++)
      if((values.at(i).rows() != rows) || (values.at(i).columns() != cols))
        return FALSE;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
