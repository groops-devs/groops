/***********************************************/
/**
* @file fileDoodsonHarmonic.h
*
* @brief Read/write harmonics coded with doodson frequencies.
*
* @author Torsten Mayer-Guerr
* @date 2010-06-10
*
*/
/***********************************************/

#ifndef __GROOPS_FILEDOODSONHARMONIC__
#define __GROOPS_FILEDOODSONHARMONIC__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_DoodsonHarmonic
static const char *docstringDoodsonHarmonic = R"(
Ocean tides are represented as time variable gravitational potential
and is given by a fourier expansion
\begin{equation}
V(\M x,t) = \sum_{f} V_f^c(\M x)\cos(\Theta_f(t)) + V_f^s(\M x)\sin(\Theta_f(t)),
\end{equation}
where $V_f^c(\M x)$ and $V_f^s(\M x)$ are spherical harmonics.
The $\Theta_f(t)$ are the arguments of the tide constituents $f$:
\begin{equation}
\Theta_f(t) = \sum_{i=1}^6 n_f^i\beta_i(t),
\end{equation}
where $\beta_i(t)$ are the Doodson's fundamental arguments ($\tau,s,h,p,N',p_s$) and $n_f^i$
are the Doodson multipliers for the term at frequency~$f$.

To extract the potential coefficients of $V_f^c$ and $V_f^s$
for each frequency $f$ use \program{DoodsonHarmonics2PotentialCoefficients}.

See also \program{PotentialCoefficients2DoodsonHarmonics}.
)";
#endif

/***********************************************/

#include "base/matrix.h"
#include "base/doodson.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_DOODSONHARMONIC_TYPE    = "doodsonHarmonic";
constexpr UInt    FILE_DOODSONHARMONIC_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** CLASS ***********************************/

/** @brief Tides: Spherical harmonics coded with doodson frequencies. */
class DoodsonHarmonic
{
public:
  Double GM, R;
  std::vector<Doodson> doodson;
  std::vector<Matrix>  cnmCos, snmCos;
  std::vector<Matrix>  cnmSin, snmSin;

/// Default Constructor.
DoodsonHarmonic() {}

/// Constructor.
explicit DoodsonHarmonic(Double _GM, Double _R, const std::vector<Doodson> &_doodson,
                         const std::vector<Matrix> &_cnmCos, const std::vector<Matrix> &_snmCos,
                         const std::vector<Matrix> &_cnmSin, const std::vector<Matrix> &_snmSin)
        : GM(_GM), R(_R), doodson(_doodson), cnmCos(_cnmCos), snmCos(_snmCos), cnmSin(_cnmSin), snmSin(_snmSin) {}
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const DoodsonHarmonic &x);
template<> void load(InArchive  &ar, DoodsonHarmonic &x);

/** @brief Write into a DoodsonHarmonic file. */
void writeFileDoodsonHarmonic(const FileName &fileName, const DoodsonHarmonic &x);

/** @brief Read from a DoodsonHarmonic file. */
void readFileDoodsonHarmonic(const FileName &fileName, DoodsonHarmonic &x);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
