/***********************************************/
/**
* @file gridGeograph.h
*
* @brief Geographical grid.
* @see Grid
*
* @author Torsten Mayer-Guerr
* @date 2002-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDGEOGRAPH__
#define __GROOPS_GRIDGEOGRAPH__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridGeograph = R"(
\subsection{Geograph}
The geographical grid is an equal-angular point distribution with points
located along meridians and along circles of latitude. \config{deltaLambda}
denotes the angular difference between adjacent points along meridians and
\config{deltaPhi} describes the angular difference between adjacent points
along circles of latitude. The point setting results as follows:
\begin{equation}
\lambda_i=\frac{\Delta\lambda}{2}+i\cdot\Delta\lambda\qquad\mbox{with}\qquad 0\leq i< \frac{360^\circ}{\Delta\lambda},
\end{equation}
\begin{equation}
\varphi_j=-90^\circ+\frac{\Delta\varphi}{2}+j\cdot\Delta\varphi\qquad\mbox{with}\qquad 0\leq j<\frac{180^\circ}{\Delta\varphi}.
\end{equation}
The number of grid points can be determined by
\begin{equation}
I=\frac{360^\circ}{\Delta\lambda}\cdot\frac{180^\circ}{\Delta\varphi}.
\end{equation}
The weights are calculated according to
\begin{equation}
w_i=\int\limits_{\lambda_i-\frac{\Delta\lambda}{2}}^{\lambda_i+\frac{\Delta\lambda}{2}}\int\limits_{\vartheta_i-\frac{\Delta\vartheta}{2}}^{\vartheta_i+\frac{\Delta\vartheta}{2}}=2\cdot\Delta\lambda\sin(\Delta\vartheta)\sin(\vartheta_i).
\end{equation}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/border/border.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Geographical grid.
* @ingroup gridGroup
* @see Grid */
class GridGeograph : public GridBase
{
public:
  GridGeograph(Config &config);
};

/***********************************************/

inline GridGeograph::GridGeograph(Config &config)
{
  try
  {
    Angle       deltaLambda, deltaPhi;
    Double      a, f, height;
    BorderPtr   border;

    readConfig(config, "deltaLambda",       deltaLambda, Config::MUSTSET,  "1", "");
    readConfig(config, "deltaPhi",          deltaPhi,    Config::MUSTSET,  "1", "");
    readConfig(config, "height",            height,      Config::DEFAULT,  "0.0",                  "ellipsoidal height expression (variables 'height', 'L', 'B')");
    readConfig(config, "R",                 a,           Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "major axsis of the ellipsoid/sphere");
    readConfig(config, "inverseFlattening", f,           Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "flattening of the ellipsoid, 0: sphere");
    readConfig(config, "border",            border,      Config::DEFAULT,  "", "");
    if(isCreateSchema(config)) return;

    Ellipsoid ellipsoid(a,f);
    for(Double phi=PI/2-deltaPhi/2.0; phi>=-PI/2; phi-=deltaPhi)
      for(Double lambda=-PI+deltaLambda/2.0; lambda<=PI; lambda+=deltaLambda)
      {
        if(border->isInnerPoint(Angle(lambda), Angle(phi)))
        {
          points.push_back(ellipsoid(Angle(lambda), Angle(phi), height));
          areas.push_back(deltaLambda * 2.0 * sin(deltaPhi/2.0) * cos(phi));
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
