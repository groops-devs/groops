/***********************************************/
/**
* @file gridDriscoll.h
*
* @brief Driscoll-Healy grid.
* @see Grid
*
* @author Annette Eicker
* @date 2002-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDDISCROLL__
#define __GROOPS_GRIDDISCROLL__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridDriscoll = R"(
\subsection{Driscoll}
The Driscoll-Healy grid, has equiangular spacing along the meridians as well
as along the circles of latitude. In longitudinal direction (along the parallels),
these angular differences for a given \config{dimension} $L$ coincide with those
described for the corresponding geographical grid and Gauss grid. Along the meridians,
the size of the latitudinal differences is half the size compared to the geographical
grid. This results in the following point pattern,
\begin{equation}
\begin{split}
\Delta\lambda=\frac{\pi}{L}\qquad&\Rightarrow\qquad\lambda_i=\frac{\Delta\lambda}{2}+i\cdot\Delta\lambda\qquad&\mbox{with}\qquad 0\leq i< 2L, \\
\Delta\vartheta=\frac{\pi}{2L}\qquad&\Rightarrow\qquad\vartheta_j=j\cdot\Delta\vartheta\qquad&\mbox{with}\qquad 1\leq j\leq 2L.
\end{split}
\end{equation}
Consequently, the number of grid points is
\begin{equation}
I=4\cdot L^2.
\end{equation}
The weights are given by
\begin{equation}
w_i=\Delta\lambda\frac{4}{2L}\sin(\vartheta_i)\sum_{l=0}^{L-1}\frac{\sin\left[ (2l+1)\;\vartheta_i\right] }{2l+1}.
\end{equation}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/border/border.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Driscoll-Healy grid.
* @ingroup gridGroup
* @see Grid */
class GridDriscoll : public GridBase
{
public:
  GridDriscoll(Config &config);
};

/***********************************************/

inline GridDriscoll::GridDriscoll(Config &config)
{
  try
  {
    UInt      dim;
    Double    a, f, height;
    BorderPtr border;

    readConfig(config, "dimension",         dim,        Config::MUSTSET,  "",                     "number of parallels = 2*dimension");
    readConfig(config, "height",            height,     Config::DEFAULT,  "0.0",                  "ellipsoidal height");
    readConfig(config, "R",                 a,          Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "major axsis of the ellipsoid/sphere");
    readConfig(config, "inverseFlattening", f,          Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "flattening of the ellipsoid, 0: sphere");
    readConfig(config, "border",            border,     Config::DEFAULT,  "", "");
    if(isCreateSchema(config)) return;

    Double deltaLambda = 2*PI/(2*dim); // # Meridiane gleich # Parallelkreise
    Double deltaPhi = PI/(2*dim); // # of Parallels muss gleich 2*n sein
    Double lambda;
    Double phi;
    Vector weights(2*dim);

    Ellipsoid ellipsoid(a,f);
    for(UInt i=0; i<2*dim; i++)   //2*dim Abschnitte entlang des Meridians
    {
      Double sum = 0;
      for(UInt m=0; m<dim; m++)
        sum += 1.0/(2*m+1)*sin((2*m+1)*i*PI/(2*dim));
      weights(i) = deltaLambda*4.0/(2*dim)*sin(i*deltaPhi)*sum;

      for(UInt j=0; j<2*dim; j++)  //2*dim Abschnitte entlang des Breitenkreises
      {
        lambda = j*deltaLambda+0.5*deltaLambda;
        phi    = PI/2-i*deltaPhi;

        if(border->isInnerPoint(Angle(lambda), Angle(phi)))
        {
          points.push_back(ellipsoid(Angle(lambda), Angle(phi), height));
          areas.push_back(weights(i));
        }
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
