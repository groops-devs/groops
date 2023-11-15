/***********************************************/
/**
* @file parametrizationTemporalDoodsonHarmonic.h
*
* @brief Tidal variations.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONTEMPORALDOODSONHARMONIC__
#define __GROOPS_PARAMETRIZATIONTEMPORALDOODSONHARMONIC__

// Latex documentation
#ifdef DOCSTRING_ParametrizationTemporal
static const char *docstringParametrizationTemporalDoodsonHarmonic = R"(
\subsection{DoodsonHarmonic}
The time variable function is given by a fourier expansion
\begin{equation}
  f(x,t) = \sum_{i} f_i^c(x)\cos(\Theta_i(t)) + f_i^s(x)\sin(\Theta_i(t)),
\end{equation}
where $\Theta_i(t)$ are the arguments of the tide constituents $i$
\begin{equation}
  \Theta_i(t) = \sum_{k=1}^6 n_i^k\beta_k(t),
\end{equation}
where $\beta_k(t)$ are the Doodson's fundamental arguments ($\tau,s,h,p,N',p_s$) and $n_i^k$
are the Doodson multipliers for the term at frequency~$i$.
The multipliers must be given by \configClass{doodson}{doodson} coded as Doodson number
(e.g. 255.555) or as names intoduced by Darwin (e.g. M2).

The major constituents given by \configClass{doodson}{doodson} can be used to interpolate minor tidal constituents
using the file \configFile{inputfileAdmittance}{admittance}. This file can be created with
\program{DoodsonHarmonicsCalculateAdmittance}.

The total parameter count is $2N$ with $N$ the number of doodson frequencies.
The parameters are sorted in following order: $f_1^c, f_1^s, f_2^c, \ldots$ with
the \file{parameter names}{parameterName} \verb|*:*:doodson.cos(<doodsonName>):*| and \verb|*:*:doodson.sin(<doodsonName>):*|.
)";
#endif

/***********************************************/

#include "base/doodson.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Tidal variations.
* @ingroup parametrizationTemporalGroup
* @see ParametrizationTemporal */
class ParametrizationTemporalDoodsonHarmonic : public ParametrizationTemporalBase
{
  UInt   _parameterCount;
  Matrix  doodsonMatrix;
  Matrix  admittance;
  std::vector<Doodson> majorDoodson;

public:
  ParametrizationTemporalDoodsonHarmonic(Config &config);

  Bool setInterval(const Time &/*timeStart*/, const Time &/*timeEnd*/, Bool /*estimatePerArc*/) {return FALSE;};
  UInt parameterCount() const {return _parameterCount;}
  void factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const;
  void parameterName(std::vector<ParameterName> &name) const;
};

/***********************************************/

#endif
