/***********************************************/
/**
* @file filePolygon.cpp
*
* @brief Polygons on Earth's surface.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_Polygon

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/filePolygon.h"

GROOPS_REGISTER_FILEFORMAT(Polygon, FILE_POLYGON_TYPE)

/***********************************************/

template<> void save(OutArchive &ar, const std::vector<Polygon> &polygon)
{
  UInt polygonCount = polygon.size();
  ar<<nameValue("polygonCount", polygonCount);
  for(UInt i=0; i<polygonCount; i++)
  {
    ar<<beginGroup("polygon");
    UInt pointCount = polygon.at(i).L.rows();
    ar<<nameValue("pointCount", pointCount);
    ar.comment("longitude [deg]           latitude [deg]          ");
    ar.comment("==================================================");
    for(UInt k=0; k<pointCount; k++)
    {
      Angle L(polygon.at(i).L(k));
      Angle B(polygon.at(i).B(k));
      ar<<beginGroup("point");
      ar<<nameValue("longitude", L);
      ar<<nameValue("latitude",  B);
      ar<<endGroup("point");
    }
    ar<<endGroup("polygon");
  }
}

/***********************************************/

template<> void load(InArchive &ar, std::vector<Polygon> &polygon)
{
  UInt polygonCount;
  ar>>nameValue("polygonCount", polygonCount);
  polygon.resize(polygonCount);
  for(UInt i=0; i<polygonCount; i++)
  {
    ar>>beginGroup("polygon");
    UInt pointCount;
    ar>>nameValue("pointCount", pointCount);
    polygon.at(i) = Polygon(pointCount);
    for(UInt k=0; k<pointCount; k++)
    {
      Angle L,B;
      ar>>beginGroup("point");
      ar>>nameValue("longitude", L);
      ar>>nameValue("latitude",  B);
      ar>>endGroup("point");
      polygon.at(i).L(k) = L;
      polygon.at(i).B(k) = B;
    }
    ar>>endGroup("polygon");
  }
}

/***********************************************/

void writeFilePolygon(const FileName &fileName, const std::vector<Polygon> &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_POLYGON_TYPE, FILE_POLYGON_VERSION);
    file<<nameValue("polygonList", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFilePolygon(const FileName &fileName, std::vector<Polygon> &x)
{
  try
  {
    InFileArchive file(fileName, FILE_POLYGON_TYPE, FILE_POLYGON_VERSION);
    file>>nameValue("polygonList", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
