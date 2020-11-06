/***********************************************/
/**
* @file tidesDoodsonHarmonic.h
*
* @brief Ocean tides.
* In terms of spherical harmonics.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @author Daniel Rieser
* @date 2002-12-13.
*/
/***********************************************/

#ifndef __GROOPS_TIDEDOODSONHARMONIC__
#define __GROOPS_TIDEDOODSONHARMONIC__

// Latex documentation
#ifdef DOCSTRING_Tides
static const char *docstringTidesDoodsonHarmonic = R"(
\subsection{DoodsonHarmonicTide}\label{tidesType:doodsonHarmonicTide}
The time variable potential of ocean tides is given by a fourier expansion
\begin{equation}
V(\M x,t) = \sum_{f} V_f^c(\M x)\cos(\Theta_f(t)) + V_f^s(\M x)\sin(\Theta_f(t)),
\end{equation}
where $V_f^c(\M x)$ and $V_f^s(\M x)$ are spherical harmonics expansions and are
read from the file \configFile{inputfileDoodsonHarmonic}{doodsonHarmonic}.
If set the expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusivly.
$\Theta_f(t)$ are the arguments of the tide constituents~$f$:
\begin{equation}
\Theta_f(t) = \sum_{i=1}^6 n_f^i\beta_i(t),
\end{equation}
where $\beta_i(t)$ are the Doodson's fundamental arguments ($\tau,s,h,p,N',p_s$)
and $n_f^i$ are the Doodson multipliers for the term at frequency~$f$.

The major constituents given by \configFile{inputfileDoodsonHarmonic}{doodsonHarmonic} can be used to
interpolate minor tidal constituents using the file \configFile{inputfileAdmittance}{admittance}.
This file can be created with \program{DoodsonHarmonicsCalculateAdmittance}.

After the interpolation step a selection of the computed constituents can be
choosen by \configClass{selectDoodson}{doodson}. Only these constiuents are considered for the results.
If no \configClass{selectDoodson}{doodson} is set all constituents will be used. The constituents can
be coded as Doodson number (e.g. 255.555) or as names intoduced by Darwin (e.g. M2).

The computed result is multiplied with \config{factor}.
)";
#endif

/***********************************************/

#include "base/doodson.h"
#include "classes/tides/tides.h"

/***** CLASS ***********************************/

/** @brief Ocean tides.
* @ingroup tidesGroup
* In terms of spherical harmonics.
* @see Tides */
class TidesDoodsonHarmonic : public TidesBase
{
  Double               GM, R;
  std::vector<Matrix>  cnmCos, cnmSin;
  std::vector<Matrix>  snmCos, snmSin;
  std::vector<Doodson> doodson;
  Matrix               doodsonMatrix;
  Matrix               admittance;
  UInt                 nCorr;

  Matrix interpolationFactors(const Time &time) const;

public:
  TidesDoodsonHarmonic(Config &config);

  SphericalHarmonics sphericalHarmonics(const Time &time, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                        UInt maxDegree, UInt minDegree, Double GM, Double R) const override;
  void deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Rotary3d> &rotEarth,
                   EarthRotationPtr rotation, EphemeridesPtr ephemerides, const std::vector<Double> &gravity, const Vector &hn, const Vector &ln,
                   std::vector<std::vector<Vector3d>> &disp) const override;
};

/***********************************************/

#endif
