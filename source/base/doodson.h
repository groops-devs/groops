/***********************************************/
/**
* @file doodson.h
*
* @brief Doodson arguments and multipliers.
*
* @author Torsten Mayer-Guerr
* @author Daniel Rieser
* @date 2005-07-15
*/
/***********************************************/


#ifndef __GROOPS_DOODSON__
#define __GROOPS_DOODSON__

// Latex documentation
#ifdef DOCSTRING_Doodson
static const char *docstringDoodson = R"(
\section{Doodson}\label{doodson}
This is a string which describes a tidal frequency either coded as Doodson number
or using Darwin´s name, e.g. \verb|255.555| or \verb|M2|.

The following names are defined:
\begin{itemize}
\item \verb|055.565|: \verb|om1|  \item \verb|055.575|: \verb|om2|     \item \verb|056.554|: \verb|sa|
\item \verb|056.555|: \verb|sa|   \item \verb|057.555|: \verb|ssa|     \item \verb|058.554|: \verb|sta|
\item \verb|063.655|: \verb|msm|  \item \verb|065.455|: \verb|mm|      \item \verb|073.555|: \verb|msf|
\item \verb|075.555|: \verb|mf|   \item \verb|083.655|: \verb|mstm|    \item \verb|085.455|: \verb|mtm|
\item \verb|093.555|: \verb|msq|  \item \verb|093.555|: \verb|msqm|    \item \verb|125.755|: \verb|2q1|
\item \verb|127.555|: \verb|sig1| \item \verb|127.555|: \verb|sigma1|  \item \verb|135.655|: \verb|q1|
\item \verb|137.455|: \verb|ro1|  \item \verb|137.455|: \verb|rho1|    \item \verb|145.555|: \verb|o1|
\item \verb|147.555|: \verb|tau1| \item \verb|155.655|: \verb|m1|      \item \verb|157.455|: \verb|chi1|
\item \verb|162.556|: \verb|pi1|  \item \verb|163.555|: \verb|p1|      \item \verb|164.555|: \verb|s1|
\item \verb|165.555|: \verb|k1|   \item \verb|166.554|: \verb|psi1|    \item \verb|167.555|: \verb|fi1|
\item \verb|167.555|: \verb|phi1| \item \verb|173.655|: \verb|the1|    \item \verb|173.655|: \verb|theta1|
\item \verb|175.455|: \verb|j1|   \item \verb|183.555|: \verb|so1|     \item \verb|185.555|: \verb|oo1|
\item \verb|195.455|: \verb|v1|   \item \verb|225.855|: \verb|3n2|     \item \verb|227.655|: \verb|eps2|
\item \verb|235.755|: \verb|2n2|  \item \verb|237.555|: \verb|mu2|     \item \verb|237.555|: \verb|mi2|
\item \verb|245.655|: \verb|n2|   \item \verb|247.455|: \verb|nu2|     \item \verb|247.455|: \verb|ni2|
\item \verb|253.755|: \verb|gam2| \item \verb|254.556|: \verb|alf2|    \item \verb|255.555|: \verb|m2|
\item \verb|256.554|: \verb|bet2| \item \verb|257.555|: \verb|dlt2|    \item \verb|263.655|: \verb|la2|
\item \verb|263.655|: \verb|lmb2| \item \verb|263.655|: \verb|lambda2| \item \verb|265.455|: \verb|l2|
\item \verb|271.557|: \verb|2t2|  \item \verb|272.556|: \verb|t2|      \item \verb|273.555|: \verb|s2|
\item \verb|274.554|: \verb|r2|   \item \verb|275.555|: \verb|k2|      \item \verb|283.655|: \verb|ksi2|
\item \verb|285.455|: \verb|eta2| \item \verb|355.555|: \verb|m3|      \item \verb|381.555|: \verb|t3|
\item \verb|382.555|: \verb|s3|   \item \verb|383.555|: \verb|r3|      \item \verb|435.755|: \verb|n4|
\item \verb|445.655|: \verb|mn4|  \item \verb|455.555|: \verb|m4|      \item \verb|473.555|: \verb|ms4|
\item \verb|491.555|: \verb|s4|   \item \verb|655.555|: \verb|m6|      \item \verb|855.555|: \verb|m8|
\end{itemize}
)";
#endif

/***********************************************/

#include "base/importStd.h"
#include "base/constants.h"
#include "base/matrix.h"
#include "base/time.h"

/***** CLASS ***********************************/

/** @brief Doodson arguments.
* @ingroup base */
class Doodson
{
public:
  /** @brief Multipliers. */
  std::array<Int, 6> d;

  /** @brief Default constructor. */
  Doodson() {d[0]=d[1]=d[2]=d[3]=d[4]=d[5]=0;}

  /** @brief Constructor with std::vector of multipliers. */
  explicit Doodson(const std::vector<Int> &v);

  /** @brief Constructor from Doodsoncode or Name.
  * You can give e.g. "255.555" or "M2". */
  explicit Doodson(const std::string &str);

  /** @brief Doodsoncode.
  * Coded by n1(n2+5)(n3+5).(n4+5)(n5+5)(n6+5).
  * Example: 255.555 for M2 Tide. */
  std::string code() const;

  /** @brief Name of tide.
  * If tide has no name the code is given instead.
  * Example: M2 for 255.555. */
  std::string name() const;

  /** @brief Angle of this constituent at time @a timeGPS. */
  Double thetaf(const Time &timeGPS) const;

  /** @brief Frequency of this constituent.
  * In rad/day. calling without given time J2000 is assumed. */
  Double frequency(Time timeGPS = mjd2time(J2000)) const;

  /** @brief Doodson arguments at time @a timeGPS. */
  static Vector arguments(const Time &timeGPS);

  /** @brief matrix of doodson multipliers.
  * To get a vector of angles of the constituents at time @a timeGPS call
  * @code
  * Matrix A      = Doodson::matrix(doodsonList);
  * Vector thetaf = A * Doodson::arguments(timeGPS);
  * @endcode */
  static Matrix matrix(const std::vector<Doodson> &doodson);

  /** @brief matrix with nodal corrections for constituents in &doodson at time @a timeGPS.
  * Following the specifications of IHO, HSSC Tidal and Water Level Working Group.
  * Corrections are implemented for Mm,Mtm,Mf,MSq,O1,Q1,K1,K2,M2,N2,2N2,M4.
  * For all others factor = 1 and phase correction = 0.
  *
  * Return of fu with:
  * fu(:,0) = node factors.
  * fu(:,1) = phase corrections (in radians). */
  static Matrix nodeCorr(const std::vector<Doodson> &doodson, const Time &timeGPS, const UInt &nCorr);

  /** @brief Implementation of Comparable (comparison of frequency). */
  Bool operator== (const Doodson &dood) const {return d == dood.d;}
  /** @brief Implementation of Comparable (comparison of frequency). */
  Bool operator!= (const Doodson &dood) const {return d != dood.d;}
  /** @brief Implementation of Comparable (comparison of frequency). */
  Bool operator<  (const Doodson &dood) const {return frequency() <  dood.frequency();}
};

/***********************************************/

#endif /* __GROOPS_DOODSON__ */
