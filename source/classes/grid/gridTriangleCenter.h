/***********************************************/
/**
* @file gridTriangleCenter.h
*
* @brief Triangle grid (center points).
* @see Grid
*
* @author Torsten Mayer-Guerr
* @date 2002-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDTRIANGLECENTER__
#define __GROOPS_GRIDTRIANGLECENTER__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridTriangleCenter = R"(
\subsection{TriangleCenter}
The points of the zeroth level are located at the centers of the icosahedron triangles.
To achieve a finer grid, each of the triangles is divided into four smaller triangles by
connecting the midpoints of the triangle edges. The refined grid points are again located
at the center of the triangles. Subsequently, the triangles can be further densified up to
the desired level of densification $n$, which is defined by \config{level}.

The number of global grid points for a certain level can be determined by
\begin{equation}\label{eq:numberCenter}
I=20\cdot 4^n.
\end{equation}
Thus the quantity of grid points depends exponentially on the level $n$, as with
every additional level the number of grid points quadruplicates.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/border/border.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Triangle grid (center points).
* @ingroup gridGroup
* @see Grid */
class GridTriangleCenter : public GridBase
{
public:
  GridTriangleCenter(Config &config);
};

/***********************************************/

class Triangle
{
  Vector3d p1, p2, p3; // Eckpunkte
  public:
    Triangle() {}
    Triangle(Vector3d _p1, Vector3d _p2, Vector3d _p3) : p1(_p1), p2(_p2), p3(_p3) {}
    // Durch Halbieren der Seiten 4 neue Dreiecke erzeugen
    void divide(Triangle &d1, Triangle &d2, Triangle &d3, Triangle &d4) const;
    // Schwerpunkt des Dreiecks
    Vector3d center() const {return normalize(p1+p2+p3);}
    // Flaeche des Dreiecks
    Double   area() const;
};

/***********************************************/

inline void Triangle::divide(Triangle &d1, Triangle &d2, Triangle &d3, Triangle &d4) const
{
  Vector3d p12 = normalize(p1+p2);
  Vector3d p23 = normalize(p2+p3);
  Vector3d p31 = normalize(p3+p1);

  d1 = Triangle(p1, p12, p31);
  d2 = Triangle(p2, p23, p12);
  d3 = Triangle(p3, p31, p23);
  d4 = Triangle(p12,p23, p31);
}

/***********************************************/

inline Double Triangle::area() const
{
  Vector3d e1 = normalize(crossProduct(p1,p2));
  Vector3d e2 = normalize(crossProduct(p2,p3));
  Vector3d e3 = normalize(crossProduct(p3,p1));

  Double w1 = PI-acos(inner(e1,e2));
  Double w2 = PI-acos(inner(e2,e3));
  Double w3 = PI-acos(inner(e3,e1));

  return w1+w2+w3-PI;
}

/***********************************************/

// Einen Vector von Dreiecken unterteilen
// und 4mal soviel kleinere Dreiecke zurueckliefern
inline std::vector<Triangle> divide(const std::vector<Triangle> &triangle)
{
  std::vector<Triangle> small(4*triangle.size());

  for(UInt i=0; i<triangle.size(); i++)
    triangle.at(i).divide(small.at(4*i), small.at(4*i+1), small.at(4*i+2), small.at(4*i+3));

  return small;
}

/***********************************************/
/***********************************************/

inline GridTriangleCenter::GridTriangleCenter(Config &config)
{
  try
  {
    UInt      level;
    Double    a, f;
    BorderPtr border;

    readConfig(config, "level",             level,      Config::MUSTSET,  "",                     "division of icosahedron, point count = 5*4**(n+1)");
    readConfig(config, "R",                 a,          Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "major axsis of the ellipsoid/sphere");
    readConfig(config, "inverseFlattening", f,          Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "flattening of the ellipsoid, 0: sphere");
    readConfig(config, "border",            border,     Config::DEFAULT,  "",                     "");
    if(isCreateSchema(config)) return;

    // Create global icosaedron grid
    // ----------------------------
    // 12 vertexes of icosaedron
    std::vector<Vector3d> p0(12);
    Double phi  = PI/2.0-acos((cos(72*DEG2RAD)+cos(72*DEG2RAD)*cos(72*DEG2RAD))/(sin(72*DEG2RAD)*sin(72*DEG2RAD)));
    p0.at(0)  = polar(Angle(  0.0*DEG2RAD), Angle( 90*DEG2RAD), 1.0);
    p0.at(1)  = polar(Angle(  0.0*DEG2RAD), Angle( phi), 1.0);
    p0.at(2)  = polar(Angle( 72.0*DEG2RAD), Angle( phi), 1.0);
    p0.at(3)  = polar(Angle(144.0*DEG2RAD), Angle( phi), 1.0);
    p0.at(4)  = polar(Angle(216.0*DEG2RAD), Angle( phi), 1.0);
    p0.at(5)  = polar(Angle(288.0*DEG2RAD), Angle( phi), 1.0);
    p0.at(6)  = polar(Angle( 36.0*DEG2RAD), Angle(-phi), 1.0);
    p0.at(7)  = polar(Angle(108.0*DEG2RAD), Angle(-phi), 1.0);
    p0.at(8)  = polar(Angle(180.0*DEG2RAD), Angle(-phi), 1.0);
    p0.at(9)  = polar(Angle(252.0*DEG2RAD), Angle(-phi), 1.0);
    p0.at(10) = polar(Angle(324.0*DEG2RAD), Angle(-phi), 1.0);
    p0.at(11) = polar(Angle(  0.0*DEG2RAD), Angle(-90*DEG2RAD), 1.0);

    // 3 vertexes of 20 triangles
    const UInt tri[20][3] = {{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},
                             {2,1,6},{3,2,7},{4,3,8},{5,4,9},{1,5,10},
                             {6,7,2},{7,8,3},{8,9,4},{9,10,5},{10,6,1},
                             {11,7,6},{11,8,7},{11,9,8},{11,10,9},{11,6,10}};

    // create triangles of the icosahedron
    std::vector<Triangle> triangle(20);
    for(UInt i=0; i<20; i++)
      triangle.at(i) = Triangle(p0.at(tri[i][0]), p0.at(tri[i][1]), p0.at(tri[i][2]));

    // divide triangles
    for(UInt i=0; i<level; i++)
      triangle = divide(triangle);

    // create points
    Ellipsoid ellipsoid(a,f);
    for(UInt i=0; i<triangle.size(); i++)
    {
      Vector3d point = triangle.at(i).center();
      // from sphere to ellipsoidal surface
      Double r = ellipsoid.b()/sqrt(1-ellipsoid.e()*ellipsoid.e()*cos(point.phi())*cos(point.phi()));
      point *= r;

      Angle  L,B;
      Double h;
      ellipsoid(point, L,B,h);

      if(border->isInnerPoint(L,B))
      {
        points.push_back(point);
        areas.push_back(triangle.at(i).area());
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
