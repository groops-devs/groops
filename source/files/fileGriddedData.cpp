/***********************************************/
/**
* @file fileGriddedData.cpp
*
* @brief Read/write gridded values.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_GriddedData

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileGriddedData.h"

GROOPS_REGISTER_FILEFORMAT(GriddedData, FILE_GRIDDEDDATA_TYPE)

/***********************************************/

template<> void save(OutArchive &ar, const GriddedData &x)
{
  // can save efficiently?
  std::vector<Angle>  lambda, phi;
  std::vector<Double> radius;
  if((ar.archiveType() == OutArchive::BINARY) && x.isRectangle(lambda, phi, radius))
  {
    // points
    ar<<nameValue("isRectangle", TRUE);
    ar<<nameValue("ellipsoid",   x.ellipsoid);
    ar<<nameValue("lambda",      lambda);
    ar<<nameValue("phi",         phi);
    ar<<nameValue("radius",      radius);

    // areas elements
    // -------------
    const UInt cols = lambda.size();
    std::vector<Double> dLambda(cols);
    dLambda.at(0) = std::fabs(lambda.at(1)-lambda.at(0));
    for(UInt s=1; s<cols-1; s++)
      dLambda.at(s) = std::fabs(0.5*(lambda.at(s+1)-lambda.at(s-1)));
    dLambda.at(cols-1) = std::fabs(lambda.at(cols-1)-lambda.at(cols-2));

    // \int_{B0-dB/2}^{B0+dB/2} cosB dB = cosB0 * 2*sin(dB/2)
    const UInt rows = phi.size();
    std::vector<Double> dPhi(rows);
    dPhi.at(0) = std::fabs(2*std::sin((phi.at(0)-phi.at(1))/2));
    for(UInt i=1; i<rows-1; i++)
      dPhi.at(i) = std::fabs(2*std::sin((phi.at(i-1)-phi.at(i+1))/4));
    dPhi.at(rows-1) = std::fabs(2*std::sin((phi.at(rows-2)-phi.at(rows-1))/2));

    // compare with given areas
    if(x.areas.size())
    {
      Bool differ = FALSE;
      for(UInt i=0; i<dPhi.size(); i++)
        for(UInt k=0; k<dLambda.size(); k++)
          if(std::fabs(x.areas.at(i*dLambda.size()+k) - dPhi.at(i) * dLambda.at(k)) > 1e-8)
            differ = TRUE;
      if(differ)
      {
        dLambda = x.areas;
        dPhi = {1.};
      }
    }

    ar<<nameValue("dLambda", dLambda);
    ar<<nameValue("dPhi",    dPhi);

    // values
    ar<<nameValue("valueCount", x.values.size());
    for(UInt id=0; id<x.values.size(); id++)
      for(UInt i=0; i<x.values.at(id).size(); i++)
        ar<<nameValue("value", x.values.at(id).at(i));
    return;
  } // if(isRectangle)

  const Bool hasArea = (x.areas.size() != 0);

  if(ar.archiveType() == OutArchive::BINARY)
    ar<<nameValue("isRectangle",  FALSE);
  ar<<nameValue("hasArea",    hasArea);
  ar<<nameValue("valueCount", x.values.size());
  ar<<nameValue("ellipsoid",  x.ellipsoid);
  ar<<nameValue("pointCount", x.points.size());

  // comment
  std::string str = "longitude [deg]           latitude [deg]            height [m]              ";
  if(hasArea)
    str += "  unit areas [-]           ";
   for(UInt i=0; i<x.values.size(); i++)
   {
     std::string str2 = "  data"+i%"%i"s;
     str += str2 + std::string(26-str2.size(), ' ');
   }
   ar.comment(str);
   ar.comment(std::string(str.size(), '='));

  for(UInt i=0; i<x.points.size(); i++)
  {
    Angle  L,B;
    Double h;
    x.ellipsoid(x.points.at(i), L,B,h);

    ar<<beginGroup("points");
    ar<<nameValue("longitude", L);
    ar<<nameValue("latitude",  B);
    ar<<nameValue("height",    h);
    if(hasArea)  ar<<nameValue("areas",  x.areas.at(i));
    for(UInt k=0; k<x.values.size(); k++)
      ar<<nameValue("value", x.values.at(k).at(i));
    ar<<endGroup("points");
  }
}

/***********************************************/

template<> void save(OutArchive &ar, const GriddedDataRectangular &x)
{
  // can save efficiently?
  if(ar.archiveType() == OutArchive::BINARY)
  {
    std::vector<Angle>  lambda, phi;
    std::vector<Double> radius;
    std::vector<Double> dLambda, dPhi;
    x.geocentric(lambda, phi, radius, dLambda, dPhi);
    ar<<nameValue("isRectangle", TRUE);
    ar<<nameValue("ellipsoid",   x.ellipsoid);
    ar<<nameValue("lambda",      lambda);
    ar<<nameValue("phi",         phi);
    ar<<nameValue("radius",      radius);
    ar<<nameValue("dLambda",     dLambda);
    ar<<nameValue("dPhi",        dPhi);
    ar<<nameValue("valueCount",  x.values.size());
    for(UInt id=0; id<x.values.size(); id++)
      for(UInt i=0; i<x.values.at(id).rows(); i++)
        for(UInt k=0; k<x.values.at(id).columns(); k++)
          ar<<nameValue("value", x.values.at(id)(i, k));
    return;
  }

  GriddedData griddedData;
  x.convert(griddedData);
  save(ar, griddedData);
}

/***********************************************/

template<> void load(InArchive &ar, GriddedData &x)
{
  // saved efficiently?
  if((ar.archiveType() == InArchive::BINARY) && (ar.version() >= 20200123))
  {
    Bool isRectangle;
    ar>>nameValue("isRectangle", isRectangle);
    if(isRectangle)
    {
      // points
      std::vector<Angle>  lambda, phi;
      std::vector<Double> radius;
      std::vector<Double> dLambda, dPhi;
      ar>>nameValue("ellipsoid", x.ellipsoid);
      ar>>nameValue("lambda",    lambda);
      ar>>nameValue("phi",       phi);
      ar>>nameValue("radius",    radius);
      ar>>nameValue("dLambda",   dLambda);
      ar>>nameValue("dPhi",      dPhi);

      std::vector<Double> cosL(lambda.size()), sinL(lambda.size());
      for(UInt s=0; s<lambda.size(); s++)
      {
        cosL[s] = std::cos(lambda[s]);
        sinL[s] = std::sin(lambda[s]);
      }

      x.points.resize(phi.size()*lambda.size());
      for(UInt z=0; z<phi.size(); z++)
      {
        const Double cosB = std::cos(phi[z]);
        const Double sinB = std::sin(phi[z]);
        for(UInt s=0; s<lambda.size(); s++)
          x.points[z*lambda.size()+s] = Vector3d(radius[z]*cosB*cosL[s], radius[z]*cosB*sinL[s], radius[z]*sinB);
      }

      // areas
      x.areas.resize(dLambda.size()*dPhi.size());
      for(UInt i=0; i<dPhi.size(); i++)
        for(UInt k=0; k<dLambda.size(); k++)
          x.areas.at(i*dLambda.size()+k) = dPhi.at(i) * dLambda.at(k);

      // values
      UInt valueCount;
      ar>>nameValue("valueCount", valueCount);
      x.values.resize(valueCount);
      for(UInt id=0; id<x.values.size(); id++)
      {
        x.values.at(id).resize(x.points.size());
        for(UInt i=0; i<x.values.at(id).size(); i++)
          ar>>nameValue("value", x.values.at(id).at(i));
      }
      return;
    }
  }

  Bool hasArea;
  UInt valueCount;
  UInt pointCount;

  ar>>nameValue("hasArea",    hasArea);
  ar>>nameValue("valueCount", valueCount);
  ar>>nameValue("ellipsoid",  x.ellipsoid);
  ar>>nameValue("pointCount", pointCount);

  x.points.resize(pointCount);
  x.areas.resize ((hasArea)  ? pointCount : 0);
  x.values.resize(valueCount);
  for(UInt k=0; k<x.values.size(); k++)
    x.values.at(k).resize(pointCount);

  Angle  lat, lon;
  Double h;
  for(UInt i=0; i<pointCount; i++)
  {
    ar>>beginGroup("points");
    ar>>nameValue("longitude", lat);
    ar>>nameValue("latitude",  lon);
    ar>>nameValue("height",    h);
    x.points.at(i) = x.ellipsoid(lat, lon, h);
    if(hasArea)
      ar>>nameValue("areas",  x.areas.at(i));
    for(UInt k=0; k<x.values.size(); k++)
      ar>>nameValue("value", x.values.at(k).at(i));
    ar>>endGroup("points");
  }
}

/***********************************************/

template<> void load(InArchive &ar, GriddedDataRectangular &x)
{
  // saved efficiently?
  if((ar.archiveType() == InArchive::BINARY) && (ar.version() >= 20200123))
  {
    Bool isRectangle;
    ar>>nameValue("isRectangle", isRectangle);
    if(isRectangle)
    {
      std::vector<Double> dLambda, dPhi;
      ar>>nameValue("ellipsoid", x.ellipsoid);
      ar>>nameValue("lambda",    x.longitudes);
      ar>>nameValue("phi",       x.latitudes);
      ar>>nameValue("radius",    x.heights);
      ar>>nameValue("dLambda",   dLambda);
      ar>>nameValue("dPhi",      dPhi);
      // points
      Angle lon;
      for(UInt i=0; i<x.latitudes.size(); i++)
        x.ellipsoid(polar(Angle(0), x.latitudes.at(i), x.heights.at(i)), lon, x.latitudes.at(i), x.heights.at(i));
      // values
      UInt valueCount;
      ar>>nameValue("valueCount", valueCount);
      x.values.resize(valueCount);
      for(UInt id=0; id<x.values.size(); id++)
      {
        x.values.at(id) = Matrix(x.latitudes.size(), x.longitudes.size());
        for(UInt i=0; i<x.values.at(id).rows(); i++)
          for(UInt k=0; k<x.values.at(id).columns(); k++)
            ar>>nameValue("value", x.values.at(id)(i, k));
      }
      return;
    }
  }

  GriddedData griddedData;
  load(ar, griddedData);
  if(!x.init(griddedData))
    throw(Exception("GriddedData must be a rectangle grid"));
}

/***********************************************/
/***********************************************/

void writeFileGriddedData(const FileName &fileName, const GriddedData &x)
{
  try
  {
    if(!x.isValid())
      throw(Exception("GriddedData is not valid"));
    OutFileArchive file(fileName, FILE_GRIDDEDDATA_TYPE);
    file<<nameValue("grid", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/

void writeFileGriddedData(const FileName &fileName, const GriddedDataRectangular &x)
{
  try
  {
    if(!x.isValid())
      throw(Exception("GriddedData is not valid"));
    OutFileArchive file(fileName, FILE_GRIDDEDDATA_TYPE);
    file<<nameValue("grid", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/

// backward compatibility
const char *const FILE_GRIDRECTANGULAR_TYPE = "gridRectangular";

// backward compatibility
static void readFileGriddedDataRectangular(InFileArchive &file, GriddedDataRectangular &grid)
{
  file>>beginGroup("gridRectangular");
  file>>nameValue("ellipsoid",  grid.ellipsoid);
  file>>nameValue("longitude",  grid.longitudes);
  file>>nameValue("latitude",   grid.latitudes);
  file>>nameValue("height",     grid.heights);
  file>>nameValue("value",      grid.values);
  file>>endGroup("gridRectangular");
}

/***********************************************/

void readFileGriddedData(const FileName &fileName, GriddedData &x)
{
  try
  {
    InFileArchive file(fileName, "");

    // Contain file an old rectangular grid?
    if(file.type() == FILE_GRIDRECTANGULAR_TYPE)
    {
      GriddedDataRectangular rectangular;
      readFileGriddedDataRectangular(file, rectangular);
      rectangular.convert(x);
      return;
    }

    // check type
    if(!file.type().empty() && (file.type() != FILE_GRIDDEDDATA_TYPE))
      throw(Exception("file type is '"+file.type()+"' but must be '"+FILE_GRIDDEDDATA_TYPE+"'"));

    file>>nameValue("grid", x);
    if(!x.isValid())
      throw(Exception("GriddedData is not valid"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/

void readFileGriddedData(const FileName &fileName, GriddedDataRectangular &x)
{
  try
  {
    InFileArchive file(fileName, "");

    // Contain file an old rectangular grid?
    if(file.type() == FILE_GRIDRECTANGULAR_TYPE)
    {
      readFileGriddedDataRectangular(file, x);
      return;
    }

    // check type
    if(!file.type().empty() && (file.type() != FILE_GRIDDEDDATA_TYPE))
      throw(Exception("file type is '"+file.type()+"' but must be '"+FILE_GRIDDEDDATA_TYPE+"'"));

    file>>nameValue("grid", x);
    if(!x.isValid())
      throw(Exception("GriddedData is not valid"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/
