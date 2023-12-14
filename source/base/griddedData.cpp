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

void GriddedData::init(const GriddedDataRectangular &grid)
{
  try
  {
    if(!grid.isValid())
      throw(Exception("GriddedDataRectangular is not valid"));

    std::vector<Double> radius, dLambda, dPhi;
    std::vector<Angle>  lambda, phi;
    grid.geocentric(lambda, phi, radius);
    grid.areaElements(dLambda, dPhi);

    std::vector<Double> cosL(lambda.size()), sinL(lambda.size());
    for(UInt k=0; k<lambda.size(); k++)
    {
      cosL[k] = std::cos(lambda[k]);
      sinL[k] = std::sin(lambda[k]);
    }

    ellipsoid = grid.ellipsoid;
    points.resize(phi.size()*lambda.size());
    areas.resize(phi.size()*lambda.size());
    for(UInt i=0; i<phi.size(); i++)
    {
      const Double cosPhi = std::cos(phi[i]);
      const Double sinPhi = std::sin(phi[i]);
      for(UInt k=0; k<lambda.size(); k++)
      {
        points[i*lambda.size()+k] = Vector3d(radius[i]*cosPhi*cosL[k], radius[i]*cosPhi*sinL[k], radius[i]*sinPhi);
        areas [i*lambda.size()+k] = dLambda[k]*dPhi[i];
      }
    }

    // values
    values.resize(grid.values.size());
    for(UInt idx=0; idx<values.size(); idx++)
    {
      values.at(idx).resize(phi.size()*lambda.size());
      for(UInt i=0; i<phi.size(); i++)
        for(UInt k=0; k<lambda.size(); k++)
          values[idx][i*lambda.size()+k] = grid.values[idx](i,k);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

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

    // TODO: check regional over date boundary

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

    const Angle  phi1 = points.front().phi();
    const Double r1   = points.front().r();
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
    if(!isRectangle(lambda, phi, radius))
      return FALSE;

    std::vector<Double> dx(lambda.size(), 2*PI);
    if(lambda.size() > 1)
    {
      dx.front() = std::fabs(std::remainder(lambda.at(1)-lambda.at(0), 2*PI));
      for(UInt k=1; k<lambda.size()-1; k++)
        dx.at(k) = 0.5*std::fabs(std::remainder(lambda.at(k+1)-lambda.at(k-1), 2*PI));
      dx.back() = std::fabs(std::remainder(lambda.at(lambda.size()-1)-lambda.at(lambda.size()-2), 2*PI));
    }

    // TODO: boundaries should be between ellipsoidal latitudes and not between geocentric phis
    // integral cos(phi) dPhi
    std::vector<Double> dy(phi.size(), 2.);
    if(phi.size() > 1)
    {
      dy.front() = std::fabs(std::sin(std::min(std::max(static_cast<Double>(phi.front()+0.5*(phi.at(1)-phi.at(0))), -PI), PI))-
                             std::sin(std::min(std::max(static_cast<Double>(phi.front()-0.5*(phi.at(1)-phi.at(0))), -PI), PI)));
      for(UInt i=1; i<phi.size()-1; i++)
        dy.at(i) = std::fabs(std::sin(0.5*(phi.at(i+1)+phi.at(i)))-std::sin(0.5*(phi.at(i)+phi.at(i-1))));
      dy.back() = std::fabs(std::sin(std::min(std::max(static_cast<Double>(phi.back()+0.5*(phi.at(phi.size()-1)-phi.at(phi.size()-2))), -PI), PI))-
                            std::sin(std::min(std::max(static_cast<Double>(phi.back()-0.5*(phi.at(phi.size()-1)-phi.at(phi.size()-2))), -PI), PI)));
    }

    areas.resize(points.size());
    for(UInt i=0; i<phi.size(); i++)
      for(UInt k=0; k<lambda.size(); k++)
        areas.at(i*lambda.size()+k) = dx.at(k)*dy.at(i);

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
    if(areas.size() && (areas.size() != points.size()))
      return FALSE;
    for(UInt idx=0; idx<values.size(); idx++)
      if(values.at(idx).size() != points.size())
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

    ellipsoid  = grid.ellipsoid;
    longitudes = lambda;
    latitudes.resize(phi.size());
    heights.resize(phi.size());
    Angle longitude;
    for(UInt i=0; i<phi.size(); i++)
      ellipsoid(polar(Angle(0), phi.at(i), radius.at(i)), longitude, latitudes.at(i), heights.at(i));

    values.resize( grid.values.size() );
    for(UInt idx=0; idx<values.size(); idx++)
    {
      values.at(idx) = Matrix(phi.size(), lambda.size());
      for(UInt i=0; i<phi.size(); i++)
        for(UInt k=0; k<lambda.size(); k++)
          values.at(idx)(i,k) = grid.values.at(idx).at(i*lambda.size()+k);
    }

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GriddedDataRectangular::geocentric(std::vector<Angle> &lambda, std::vector<Angle> &phi, std::vector<Double> &radius) const
{
  try
  {
    lambda = longitudes;
    phi.resize(latitudes.size());
    radius.resize(heights.size());
    for(UInt i=0; i<phi.size(); i++)
    {
      const Vector3d points = ellipsoid(Angle(0), latitudes.at(i), heights.at(i));
      phi.at(i)    = points.phi();
      radius.at(i) = points.r();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GriddedDataRectangular::cellBorders(std::vector<Double> &lambda, std::vector<Double> &phi) const
{
  try
  {
    lambda.clear();
    if(longitudes.size() > 1)
    {
      lambda.resize(longitudes.size()+1);
      lambda.front() = std::remainder(longitudes.front() - 0.5*std::remainder(longitudes.at(1)-longitudes.at(0), 2*PI), 2*PI);
      for(UInt k=0; k<longitudes.size()-1; k++)
        lambda.at(k+1) = std::remainder(longitudes.at(k) + 0.5*std::remainder(longitudes.at(k+1)-longitudes.at(k), 2*PI), 2*PI);
      lambda.back() = std::remainder(longitudes.back() + 0.5*std::remainder(longitudes.at(longitudes.size()-1)-longitudes.at(longitudes.size()-2), 2*PI), 2*PI);
    }
    else if(longitudes.size() == 1)
      lambda = {std::remainder(longitudes.front()-PI, 2*PI), std::remainder(longitudes.front()+PI, 2*PI)};

    phi.clear();
    if(latitudes.size() > 1)
    {
      phi.resize(latitudes.size()+1);
      phi.front() = std::min(std::max(static_cast<Double>(latitudes.front()-0.5*(latitudes.at(1)-latitudes.at(0))), -PI), PI);
      for(UInt i=0; i<latitudes.size()-1; i++)
        phi.at(i+1) = 0.5 * (latitudes.at(i) + latitudes.at(i+1));
      phi.back() = std::min(std::max(static_cast<Double>(latitudes.back()+0.5*(latitudes.at(latitudes.size()-1)-latitudes.at(latitudes.size()-2))), -PI), PI);
    }
    else if(latitudes.size() == 1)
      phi = {PI/2, -PI/2};
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GriddedDataRectangular::areaElements(std::vector<Double> &dLambda, std::vector<Double> &dPhi) const
{
  try
  {
    std::vector<Double> lambda, phi;
    cellBorders(lambda, phi);
    for(UInt i=0; i<phi.size(); i++)
      phi.at(i) = ellipsoid(Angle(0.), Angle(phi.at(i)), 0.).phi();  // geocentric

    dLambda.resize(longitudes.size());
    for(UInt k=0; k<dLambda.size(); k++)
      dLambda.at(k) = std::fabs(std::remainder(lambda.at(k+1)-lambda.at(k), 2*PI));
    if(longitudes.size() == 1)
      dLambda.front() = 2*PI;

    dPhi.resize(latitudes.size());
    for(UInt i=0; i<dPhi.size(); i++)
      dPhi.at(i) = std::fabs(std::sin(phi.at(i+1))-std::sin(phi.at(i))); // area = integral cos(phi) dPhi

    // total area
    return std::accumulate(dLambda.begin(), dLambda.end(), 0.) * std::accumulate(dPhi.begin(), dPhi.end(), 0.);
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
    if(heights.size() != latitudes.size())
      return FALSE;
    for(UInt idx=0; idx<values.size(); idx++)
      if((values.at(idx).rows() != latitudes.size()) || (values.at(idx).columns() != longitudes.size()))
        return FALSE;
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
