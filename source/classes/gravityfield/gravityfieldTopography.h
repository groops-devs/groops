/***********************************************/
/**
* @file gravityfieldTopography.h
*
* @brief Gravity field from topographic masses.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2013-02-10
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDTOPOGRAPHY__
#define __GROOPS_GRAVITYFIELDTOPOGRAPHY__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldTopography = R"(
\subsection{Topography}\label{gravityfieldType:topography}
The gravity is integrated from a topographic mass distribution.
For each grid point in \configFile{inputfileGridRectangular}{griddedData} a prisma with
\config{density} is assumed. The horizontal extension is computed from the grid spacing
and the vertical extension is given by \config{radialLowerBound}
and \config{radialUpperBound} above ellipsoid. All values are expressions and computed
for each point with given data in the grid file. The standard variables for grids
are available, see~\reference{dataVariables}{general.parser:dataVariables}.

Example: The grid file contains the orthometric height of the topography in the first
column, the geoid height in the second and the mean density of each prism in the third
column. In this case the following settings should be used:
\begin{itemize}
\item \config{radialUpperBound} = \verb|data0+data1|,
\item \config{radialLowerBound} = \verb|data1|,
\item \config{density} = \verb|data2|.
\end{itemize}

As the prim computation is time consuming a maximum distance around the evaluation point
can defined with \config{distancePrism}. Afterwards a simplified radial line
(the prism mass is concentrated to a line in the center) is used up to
a distance of \config{distanceLine}. At last the prim is approximated by a point mass
in the center up to a distance \config{distanceMax} (if set). Prisms nearby the evaluation
point can be excluded with \config{distanceMin}.
)";
#endif

/***********************************************/

#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Gravity field from topographic masses.
* @ingroup gravityfieldGroup
* @see Gravityfield */
class GravityfieldTopography : public GravityfieldBase
{
  Double              factor;
  UInt                rows, cols;
  std::vector<Angle>  lambda, phi;
  std::vector<Double> dLambda, dPhi;
  Vector              sinL, cosL, sinB, cosB;
  Matrix              rLower, rUpper, rho;
  Double              cosPsiMin, cosPsiPrism, cosPsiLine, cosPsiMax;
  Bool                testRectangle;
  Bool                isPhiAscending, isLambdaAscending;

  inline void findRectangle(const Vector3d &point, UInt &colsMin, UInt &rowsMin, UInt &colsMax, UInt &rowsMax) const;

public:
  GravityfieldTopography(Config &config);

  // Einfluss des Referenzfieldes im Aufpunkt:
  Double   potential      (const Time &time, const Vector3d &point) const;
  Double   radialGradient (const Time &time, const Vector3d &point) const;
  Vector3d gravity        (const Time &time, const Vector3d &point) const;
  Tensor3d gravityGradient(const Time &time, const Vector3d &point) const;
  Vector3d deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const;
  void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                           const Vector &hn, const Vector &ln, std::vector< std::vector<Vector3d> > &disp) const;
  void     variance       (const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const;
  SphericalHarmonics sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const;
};

/***********************************************/

#endif /* __GROOPS_GRAVITYFIELD__ */
