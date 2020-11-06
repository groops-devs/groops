/***********************************************/
/**
* @file gridCorput.h
*
* @brief Corput distribution.
* @see Grid
* Pseudo random distribution.
*
* @author Torsten Mayer-Guerr
* @date 2002-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDCORPUT__
#define __GROOPS_GRIDCORPUT__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridCorput = R"(
\subsection{Corput}
This kind of grid distributes an arbitrarily chosen number of $I$ points
(defined by \config{globalPointsCount}) following a recursive, quasi random sequence.
In longitudinal direction the pattern follows
\begin{equation}
\Delta\lambda=\frac{2\pi}{I}\qquad\Rightarrow\qquad\frac{\Delta\lambda}{2}+\lambda_i=i\cdot\Delta\lambda\qquad\mbox{with}\qquad 1\leq i\leq I.
\end{equation}
This implies that every grid point features a unique longitude, with equi-angular
longitudinal differences.

The polar distance in the form $t_i=\cos\vartheta_i$ for each point is determined
by the following recursive sequence:
\begin{itemize}
\item Starting from an interval $t\in[-1,1]$.
\item If $I=1$, then the midpoint of the interval is returned as result of
the sequence, and the sequence is terminated.
\item If the number of points is uneven, the  midpoint is included into the list of $t_i$.
\item Subsequently, the interval is bisected into an upper and lower half,
       and the sequence is called for both halves.
\item $t$ from upper and lower half are alternately sorted into the list of $t_i$.
\item The polar distances are calculated by
\begin{equation}
\vartheta_i=\arccos\, t_i.
\end{equation}
\end{itemize}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/border/border.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Corput distribution.
* Pseudo random distribution.
* @ingroup gridGroup
* @see Grid */
class GridCorput : public GridBase
{
  static std::vector<Double> sequence(Double min, Double max, UInt n);

public:
  GridCorput(Config &config);
};

/***********************************************/

inline GridCorput::GridCorput(Config &config)
{
  try
  {
    UInt      numberOfGlobalPoints;
    Double    a, f, height;
    BorderPtr border;

    readConfig(config, "globalPointsCount", numberOfGlobalPoints, Config::MUSTSET,  "", "");
    readConfig(config, "height",            height,     Config::DEFAULT,  "0.0",                  "ellipsoidal height");
    readConfig(config, "R",                 a,          Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "major axsis of the ellipsoid/sphere");
    readConfig(config, "inverseFlattening", f,          Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "flattening of the ellipsoid, 0: sphere");
    readConfig(config, "border",            border,     Config::DEFAULT,  "", "");
    if(isCreateSchema(config)) return;

    Ellipsoid ellipsoid(a,f);
    std::vector<Double> t = sequence(-1.0, 1.0, numberOfGlobalPoints);
    for(UInt i=0; i<numberOfGlobalPoints; i++)
    {
      Angle L(fmod(2*PI*(i+0.5)/numberOfGlobalPoints+PI,2*PI)-PI);
      Angle B(PI/2-acos(t.at(i)));
      if(border->isInnerPoint(L,B))
      {
        points.push_back(ellipsoid(L, B, height));
        areas.push_back(4*PI/numberOfGlobalPoints); // average area
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Double> GridCorput::sequence(Double min, Double max, UInt n)
{
  Double middle = (max+min)/2.0;

  UInt nr=0;
  std::vector<Double> t(n);

  // Ist ein Punkt in der Mitte?
  if(n%2==1) t.at(nr++) = middle;

  // Wenn nur ein Punkt, dann schon fertig
  if(n==1) return t;

  // Links und rechts erzeugen
  std::vector<Double> left  = sequence(min, middle, n/2);
  std::vector<Double> right = sequence(middle, max, n/2);

  // Abwechselnd sortieren
  for(UInt i=0; i<n/2; i++)
  {
    t.at(nr++) = left.at(i);
    t.at(nr++) = right.at(i);
  }

  return t;
}

/***********************************************/

#endif /* __GROOPS__ */
