/***********************************************/
/**
* @file gridTriangleVertex.h
*
* @brief Triangle grid (Vertcies)
* @see Grid
*
* @author Torsten Mayer-Guerr
* @date 2002-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDTRIANGLEVERTEX__
#define __GROOPS_GRIDTRIANGLEVERTEX__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridTriangleVertex = R"(
\subsection{TriangleVertex}
The zeroth level of densification
coincides with the 12 icosahedron vertices, as displayed in the upper left part
of Fig.~\ref{fig:triangle_grid}. Then, depending on the envisaged densification,
each triangle edge is divided into $n$ parts, illustrated in the upper right
part of Fig.~\ref{fig:triangle_grid}. The new nodes on the edges are then connected
by arcs of great circles parallel to the triangle edges. The intersections of
each three corresponding parallel lines become nodes of the densified grid as well.
As in case of a spherical triangle those three connecting lines do not exactly
intersect in one point, the center of the resulting triangle is used as location
for the new node (lower left part of Fig.~\ref{fig:triangle_grid}). The lower right
side of Fig.~\ref{fig:triangle_grid} finally shows the densified triangle vertex
grid for a level of $n=3$. The number of grid points in dependence of the chosen
level of densification can be calculated by
\begin{equation}\label{eq:numberVertex}
I=10\cdot(n+1)^2+2.
\end{equation}

\fig{!hb}{0.6}{icogrid}{fig:triangle_grid}{TriangleVertex grid.}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/border/border.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Triangle grid (Vertcies).
* @ingroup gridGroup
* @see Grid */
class GridTriangleVertex : public GridBase
{
  static void divideEdge(Vector3d p1, Vector3d p2, UInt level, std::vector<Vector3d> &points);
  static void divideTriangle(Vector3d p1, Vector3d p2, Vector3d p3, UInt level, std::vector<Vector3d> &points);

public:
  GridTriangleVertex(Config &config);
};

/***********************************************/

inline GridTriangleVertex::GridTriangleVertex(Config &config)
{
  try
  {
    UInt      level;
    Double    a, f;
    BorderPtr border;

    readConfig(config, "level",             level,      Config::MUSTSET,  "",                     "division of icosahedron, point count = 10*(n+1)**2+2");
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

    // start and end vertex of 30 edges
    const UInt edge[30][2] = {{0,1},{0,2},{0,3},{0,4},{0,5},
                              {1,2},{2,3},{3,4},{4,5},{5,1},
                              {1,6},{2,7},{3,8},{4,9},{5,10},
                              {6,2},{7,3},{8,4},{9,5},{10,1},
                              {6,7},{7,8},{8,9},{9,10},{10,6},
                              {11,6},{11,7},{11,8},{11,9},{11,10}};

    // divides edges
    for(UInt i=0; i<30; i++)
      divideEdge(p0.at(edge[i][0]), p0.at(edge[i][1]), level, p0);

    // create inner points
    for(UInt i=0; i<20; i++)
      divideTriangle( p0.at(tri[i][0]), p0.at(tri[i][1]), p0.at(tri[i][2]), level, p0);

    // p0 contain now the global point distribution on unit sphere

    Ellipsoid ellipsoid(a,f);
    for(UInt i=0; i<p0.size(); i++)
    {
      // from sphere to ellipsoidal surface
      Double r = ellipsoid.b()/sqrt(1-ellipsoid.e()*ellipsoid.e()*cos(p0.at(i).phi())*cos(p0.at(i).phi()));
      p0.at(i) *= r;

      Angle  L,B;
      Double h;
      ellipsoid(p0.at(i), L,B,h);

      if(border->isInnerPoint(L,B))
      {
        points.push_back(p0.at(i));
        areas.push_back(4*PI / p0.size());
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Erzeugt level gleichabstandige Punkte zwischen Anfpkt und Endpkt
inline void GridTriangleVertex::divideEdge(Vector3d p1, Vector3d p2, UInt level, std::vector<Vector3d> &points)
{
  Double psi = acos(inner(p1, p2))/(level+1);

  Vector3d pPerpendicular = normalize(crossProduct(crossProduct(p1, p2), p1));

  for(UInt i=1; i<=level; i++)
    points.push_back( cos(i*psi)*p1 + sin(i*psi)*pPerpendicular );
}

/***********************************************/

// Innere Punkte im Dreieck erzeugen
inline void GridTriangleVertex::divideTriangle(Vector3d p1, Vector3d p2, Vector3d p3, UInt level, std::vector<Vector3d> &points)
{
  std::vector<Vector3d> edge1, edge2, edge3;

  divideEdge(p1, p2, level, edge1);
  divideEdge(p2, p3, level, edge2);
  divideEdge(p3, p1, level, edge3);

  // innere punkte erzeugen
  for(UInt i=1; i<level; i++)
    for(UInt k=0; k<i; k++)
    {
      // Verbindungsgraden erzeugen
      Vector3d line13 = crossProduct(edge1.at(i),     edge3.at(level-1-i));
      Vector3d line12 = crossProduct(edge1.at(i-1-k), edge2.at(level-i+k));
      Vector3d line23 = crossProduct(edge2.at(k),     edge3.at(level-1-k));

      // Geraden schneiden (Schnittpunkt)
      Vector3d p1 = normalize(crossProduct(line13, line12));
      Vector3d p2 = normalize(crossProduct(line23, line13));
      Vector3d p3 = normalize(crossProduct(line23, line12));

      // Punkte mitteln
      Vector3d pMean = normalize(p1 + p2 + p3);
      points.push_back( -pMean );
    }
}

/***********************************************/

#endif /* __GROOPS__ */
