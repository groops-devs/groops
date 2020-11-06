/***********************************************/
/**
* @file gridGauss.h
*
* @brief Gauss-Legendre grid.
* @see Grid
*
* @author Annette Eicker
* @date 2002-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDGAUSS__
#define __GROOPS_GRIDGAUSS__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridGauss = R"(
\subsection{Gauss}
 The grid features equiangular spacing along circles of latitude with
 \config{parallelsCount} defining the number $L$ of the parallels.
\begin{equation}
\Delta\lambda=\frac{\pi}{L}\qquad\Rightarrow\qquad\lambda_i=\frac{\Delta\lambda}{2}+i\cdot\Delta\lambda\qquad\mbox{with}\qquad 0\leq i< 2L.
\end{equation}
Along the meridians the points are located at $L$ parallels at
the $L$ zeros $\vartheta_j$ of the Legendre polynomial of degree~$L$,
\begin{equation}
P_L(\cos\vartheta_j)=0.
\end{equation}
Consequently, the number of grid points sums up to
\begin{equation}
I=2\cdot L^2.
\end{equation}
The weights can be calculated according to
\begin{equation}
w_i(L)=\Delta\lambda\frac{2}{(1-t_i^2)(P'_{L}(\cos(\vartheta _i)))^2},\label{weights}
\end{equation}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/legendrePolynomial.h"
#include "config/config.h"
#include "classes/border/border.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Gauss-Legendre grid.
* @ingroup gridGroup
* @see Grid */
class GridGauss : public GridBase
{
public:
  GridGauss(Config &config);
};

/***********************************************/

inline GridGauss::GridGauss(Config &config)
{
  try
  {
    UInt      numberOfParallels;
    Double    a, f;
    BorderPtr border;

    readConfig(config, "parallelsCount",    numberOfParallels, Config::MUSTSET, "", "");
    readConfig(config, "R",                 a,                 Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "major axsis of the ellipsoid/sphere");
    readConfig(config, "inverseFlattening", f,                 Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "flattening of the ellipsoid, 0: sphere");
    readConfig(config, "border",            border,            Config::DEFAULT,  "", "");
    if(isCreateSchema(config)) return;

    Ellipsoid ellipsoid(a,f);
    Double deltaLambda = 2*PI/(2*numberOfParallels); // # Meridiane: 2*anz_nullst
    // number of zeros is number of latitudes
    Vector zeros, weights;
    LegendrePolynomial::zeros(numberOfParallels, zeros, weights);
    for(UInt i=0; i<numberOfParallels; i++)
      for(Double lambda=-PI+deltaLambda/2.0; lambda<=PI; lambda+=deltaLambda)
      {
        // create point on unit sphere
        Vector3d point = polar(Angle(lambda), Angle(PI/2-acos(zeros(i))), 1.0);
        // from sphere to ellipsoidal surface
        Double r = ellipsoid.b()/sqrt(1-ellipsoid.e()*ellipsoid.e()*cos(point.phi())*cos(point.phi()));
        point *= r;

        Angle  L,B;
        Double h;
        ellipsoid(point, L,B,h);

        if(border->isInnerPoint(L,B))
        {
          points.push_back(point);
          areas.push_back(deltaLambda * weights(i));
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
