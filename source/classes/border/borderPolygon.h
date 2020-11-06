/***********************************************/
/**
* @file borderPolygon.h
*
* @brief Borders described by spherical polygons.
* @see Border
*
* @author Andreas Kvas
* @date 2020-09-28
*
*/
/***********************************************/

#ifndef __GROOPS_BORDERPOLYGON__
#define __GROOPS_BORDERPOLYGON__

// Latex documentation
#ifdef DOCSTRING_Border
static const char *docstringBorderPolygon = R"(
\subsection{Polygon}\label{borderType:polygon}
The region is defined by \configFile{inputfilePolygon}{polygon}
containing one or more polygons given in longitude and latitude.
An additional \config{buffer} around the polygon can be defined.
Use a negative value to shrink the polygon area.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/filePolygon.h"
#include "classes/border/border.h"

/***** CLASS ***********************************/

/** @brief polygon.
* @ingroup borderGroup
* @see Border */
class BorderPolygon : public BorderBase
{
  Ellipsoid             ellipsoid;
  std::vector<std::vector<Vector3d>> vertices;
  std::vector<Vector3d> centroid;
  std::vector<Double>   capThreshold;
  Bool                  exclude;
  Double                buffer;

  Bool inPolygon(const Vector3d &testPoint, UInt polyNo) const;
  Bool inBuffer (const Vector3d &testPoint, UInt polyNo) const;

public:
  BorderPolygon(Config &config);

  Bool isInnerPoint(Angle lambda, Angle phi) const;
  Bool isExclude() const {return exclude;}
};

/***********************************************/

inline BorderPolygon::BorderPolygon(Config &config)
{
  FileName polygonName;

  readConfig(config, "inputfilePolygon", polygonName, Config::MUSTSET, "{groopsDataDir}/border/", "");
  readConfig(config, "buffer",           buffer,      Config::DEFAULT, "0", "buffer around polygon [km], <0: inside");
  readConfig(config, "exclude",          exclude,     Config::DEFAULT, "0", "dismiss points inside");
  if(isCreateSchema(config)) return;

  std::vector<Polygon> polygon;
  readFilePolygon(polygonName, polygon);

  ellipsoid = Ellipsoid(DEFAULT_GRS80_a, DEFAULT_GRS80_f);

  // find vertices and centroid for every polygon:
  centroid.resize(polygon.size());
  vertices.resize(polygon.size());
  capThreshold.resize(polygon.size(), 1.0);
  for(UInt i=0; i<polygon.size(); i++)
  {
    const UInt vertexCount = polygon.at(i).L.rows();
    vertices.at(i).resize(vertexCount);
    for(UInt k=0; k<vertexCount; k++)
    {
      vertices.at(i).at(k) = normalize(ellipsoid(Angle(polygon.at(i).L(k)), Angle(polygon.at(i).B(k)), 0.0));
      centroid.at(i)      += vertices.at(i).at(k);
    }
    centroid.at(i).normalize();

    for(auto &v : vertices.at(i))
      capThreshold.at(i) = std::min(capThreshold.at(i), std::cos(std::acos(inner(centroid.at(i), v)) + std::fabs(buffer)/DEFAULT_R*1e3));
  }
}

/***********************************************/

inline Bool BorderPolygon::inPolygon(const Vector3d &testPoint, UInt polyNo) const
{
  if(inner(centroid.at(polyNo), testPoint) < capThreshold.at(polyNo))
    return FALSE;

  const Vector3d p   = crossProduct(testPoint, -centroid.at(polyNo));
  const Vector3d axp = crossProduct(testPoint, p);
  const Vector3d cxp = crossProduct(-centroid.at(polyNo), p);
  const UInt     vertexCount = vertices.at(polyNo).size();
  UInt crossingCount = 0;
  for(UInt k=0; k<vertexCount; k++)
  {
    const Vector3d q = crossProduct(vertices.at(polyNo).at(k), vertices.at(polyNo).at((k+1)%vertexCount));
    const Vector3d t = crossProduct(p, q);
    if(t.norm() == 0.0)
      continue;

    Bool sign = std::signbit(-inner(t, axp));
    if((sign == std::signbit( inner(t, cxp))) &&
       (sign == std::signbit(-inner(t, crossProduct(vertices.at(polyNo).at(k), q)))) &&
       (sign == std::signbit( inner(t, crossProduct(vertices.at(polyNo).at((k+1)%vertexCount), q)))))
      crossingCount++;
  }

  return Bool(crossingCount%2);
}

/***********************************************/

inline Bool BorderPolygon::inBuffer(const Vector3d &testPoint, UInt polyNo) const
{
  if(inner(centroid.at(polyNo), testPoint) < capThreshold.at(polyNo))
    return FALSE;

  const Double cosBuffer   = std::cos(buffer*1e3/DEFAULT_R);
  const UInt   vertexCount = vertices.at(polyNo).size();
  for(UInt k=0; k<vertexCount; k++)
    if(cosBuffer <= inner(vertices.at(polyNo).at(k), testPoint))
      return TRUE;

  const Double sinBuffer = std::sin(std::fabs(buffer)*1e3/DEFAULT_R);
  for(UInt k=0; k<vertexCount; k++)
  {
    const Vector3d n = crossProduct(vertices.at(polyNo).at(k), vertices.at(polyNo).at((k+1)%vertexCount));
    if((inner(n, crossProduct(testPoint, vertices.at(polyNo).at(k))) <= 0) &&
       (inner(n, crossProduct(testPoint, vertices.at(polyNo).at((k+1)%vertexCount))) >= 0) &&
       (std::fabs(inner(testPoint, normalize(n))) <= sinBuffer))
      return TRUE;
  }

  return FALSE;
}

/***********************************************/

inline Bool BorderPolygon::isInnerPoint(Angle lambda, Angle phi) const
{
  const Vector3d testPoint = normalize(ellipsoid(lambda, phi, 0));

  if(buffer)
    for(UInt polyNo=0; polyNo<vertices.size(); polyNo++)
      if(inBuffer(testPoint, polyNo))
        return (buffer > 0);

  UInt count = 0;
  for(UInt polyNo=0; polyNo<vertices.size(); polyNo++)
    if(inPolygon(testPoint, polyNo))
      count++;
  return ((count%2)==1);
}

/***********************************************/

#endif /* __GROOPS_BORDER__ */
