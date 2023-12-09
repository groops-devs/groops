/***********************************************/
/**
* @file gridReuter.h
*
* @brief Reuter grid.
* @see Grid
*
* @author Annette Eicker
* @date 2002-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDREUTER__
#define __GROOPS_GRIDREUTER__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridReuter = R"(
\subsection{Reuter}
The Reuter grid features equi-distant spacing along the meridians determined
by the control parameter~$\gamma$ according to
\begin{equation}
\Delta\vartheta=\frac{\pi}{\gamma}\qquad\Rightarrow\vartheta_j=j\Delta\vartheta,\qquad\mbox{with}\qquad 1\leq j\leq \gamma-1.
\end{equation}
Thus $\gamma+1$ denotes the number of points per meridian, as the two poles
are included in the point distribution as well. Along the circles of latitude,
the number of grid points decreases with increasing latitude in order to achieve
an evenly distributed point pattern. This number is chosen, so that the points
along each circle of latitude have the same spherical distance as two adjacent
latitudes. The resulting relationship is given by
\begin{equation}\label{eq:sphericalDistance}
\Delta\vartheta=\arccos\left( \cos^2\vartheta_j+\sin^2\vartheta_j\cos\Delta\lambda_j\right).
\end{equation}
The left hand side of this equation is the spherical distance between adjacent
latitudes, the right hand side stands for the spherical distance between two points
with the same polar distance $\vartheta_j$ and a longitudinal difference of
$\Delta\lambda_i$. This longitudinal distance can be adjusted depending on
$\vartheta_j$ to fulfill Eq.~\eqref{eq:sphericalDistance}. The resulting
formula for $\Delta\lambda_i$ is
\begin{equation}\label{eq:deltaLambdai}
\Delta\lambda_j=\arccos\left( \frac{\sin\Delta\vartheta -\cos^2\vartheta_j}{\sin^2\vartheta_j}\right).
\end{equation}
The number of points~$\gamma_j$ for each circle of latitude can then be determined by
\begin{equation}\label{eq:gammai}
\gamma_j=\left[ \frac{2\pi}{\Delta\lambda_j}\right] .
\end{equation}
Here the Gauss bracket $[x]$ specifies the largest integer equal to or less than $x$.
The longitudes are subsequently determined by
\begin{equation}
\lambda_{ij}=\frac{\Delta\lambda_j}{2}+i\cdot(2\pi/\gamma_j),\qquad\mbox{with}\qquad 0\leq i< \gamma_j.
\end{equation}
The number of grid points can be estimated by
\begin{equation}\label{eq:numberReuter}
I=\leq 2+\frac{4}{\pi}\gamma^2,
\end{equation}
The $\leq$ results from the fact that the $\gamma_j$ are restricted to integer values.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/border/border.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Reuter grid.
* @ingroup gridGroup
* @see Grid */
class GridReuter : public GridBase
{
public:
  GridReuter(Config &config);
};

/***********************************************/

inline GridReuter::GridReuter(Config &config)
{
  try
  {
    UInt      gamma;
    Double    a, f, height;
    BorderPtr border;

    readConfig(config, "gamma",             gamma,  Config::MUSTSET,  "",                     "number of parallels");
    readConfig(config, "height",            height, Config::DEFAULT,  "0.0",                  "ellipsoidal height");
    readConfig(config, "R",                 a,      Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "major axsis of the ellipsoid/sphere");
    readConfig(config, "inverseFlattening", f,      Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "flattening of the ellipsoid, 0: sphere");
    readConfig(config, "border",            border, Config::DEFAULT,  "", "");
    if(isCreateSchema(config)) return;

    Angle     deltaB(PI/gamma);
    Ellipsoid ellipsoid(a,f);

    // north pole
    Angle L(0.0);
    Angle B(PI/2);
    if(border->isInnerPoint(L,B))
    {
      points.push_back(ellipsoid(L, B, height));
      areas.push_back(2*PI*(1-std::cos(deltaB/2.0)));
    }

    // other points
    for(UInt i=1; i<gamma; i++)
    {
      Double theta   = i* PI/gamma;
      UInt   gamma_i = static_cast<UInt>(std::floor(2*PI/std::acos((std::cos(PI/gamma)-std::cos(theta)*std::cos(theta))/(std::sin(theta)*std::sin(theta)))));
      Angle  B(PI/2-theta);
      Angle  deltaL(2*PI/gamma_i);

      for(UInt j=1; j<=gamma_i; j++)
      {
        Angle L(fmod((j+0.5)*(2*PI/gamma_i)+PI, 2*PI)-PI);
        if(border->isInnerPoint(L,B))
        {
          points.push_back(ellipsoid(L, B, height));
          areas.push_back(deltaL * 2.0 * std::sin(deltaB/2.0) * std::cos(B));
        }
      }
    }

    // south pole
    L = Angle(0.0);
    B = Angle(-PI/2);
    if(border->isInnerPoint(L,B))
    {
      points.push_back(ellipsoid(L, B, height));
      areas.push_back(2*PI*(1-std::cos(deltaB/2.0)));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
