/***********************************************/
/**
* @file earthRotationIers2003.h
*
* @brief According to IERS2003 conventions.
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATIONIERS2003__
#define __GROOPS_EARTHROTATIONIERS2003__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotationIers2003 = R"(
\subsection{Iers2003}
This class realize the transformation according to IERS2003 conventions
given by the \emph{International Earth Rotation and Reference Systems Service} (IERS).
A file with the earth orientation parameter is needed (\configFile{inputfileEOP}{earthOrientationParameter}).

The following subroutines are used:
\begin{itemize}
\item BPN2000.f,
\item ERA2000.f,
\item pmsdnut.f,
\item POM2000.f,
\item SP2000.f,
\item T2C2000.f,
\item XYS2000A.f
\end{itemize}
from \url{ftp://maia.usno.navy.mil/conv2000/chapter5/} and
\begin{itemize}
\item orthoeop.f
\end{itemize}
from \url{ftp://maia.usno.navy.mil/conv2000/chapter8/}
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief According to IERS2003 conventions.
 * In addtion to those interpolated values from the input file, 
 * - diurnal and semi-diurnal variations from ocean tides are considered by the IERS routine \c ORTHO_EOP.F for
 *   x-pole, y-pole and dUT1
 * - librations are considered by the routine \c PMSDNUT2.F for x-pole and y-pole
 * - librations are considered by the routine \c UTLIBR.F for dUT1 and LOD 
 * 
 * sp is calculated with the ERFA function \c eraSp00. X,Y, and S are calculated by the 
 * precession & nutation model with the ERFA function \c eraXys00a.
* @ingroup earthRotationGroup
* @see EarthRotation */
class EarthRotationIers2003 : public EarthRotation
{
  Polynomial        polynomial;
  Matrix            EOP;

public:
  EarthRotationIers2003(Config &config);

  void earthOrientationParameter(const Time &timeGPS, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const;
};

/***********************************************/

#endif /* __GROOPS_EARTHROTATION__ */
