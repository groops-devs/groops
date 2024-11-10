/***********************************************/
/**
* @file fileEarthOrientationParameter.h
*
* @brief Read/write EarthOrientationParameter.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-09
*
*/
/***********************************************/

#ifndef __GROOPS_FILEEARTHORIENTATIONPARAMETER__
#define __GROOPS_FILEEARTHORIENTATIONPARAMETER__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_EarthOrientationParameter
static const char *docstringEarthOrientationParameter = R"(
Earth Orientation Parameter (EOP) as provided by the International Earth Rotation and Reference Systems Service (IERS) (e.g \verb|EOP 14 C04 (IAU2000A)|).

See \program{IersC04IAU2000EarthOrientationParameter}, \program{IersRapidIAU2000EarthOrientationParameter}.

\begin{verbatim}
groops earthOrientationParameter version=20200123
       9641 # number of epochs
# UTC [MJD]                 xp [arcsec]               yp [arcsec]               deltUT [sec]              LOD [sec]                 dX [arcsec]               dY [arcsec]
# ====================================================================================================================================================================================
  5.894700000000000000e+04  5.690599999999999825e-02  4.099130000000000273e-01 -2.316246000000000138e-01  1.636400000000000094e-03 -2.900000000000000017e-05  5.800000000000000034e-05
  5.894800000000000000e+04  5.771400000000000141e-02  4.110159999999999925e-01 -2.332083000000000073e-01  1.520099999999999923e-03 -6.000000000000000152e-05  2.199999999999999943e-05
  5.894900000000000000e+04  5.813000000000000111e-02  4.120099999999999874e-01 -2.346157000000000104e-01  1.293099999999999935e-03 -7.200000000000000182e-05  3.199999999999999855e-05
  5.895000000000000000e+04  5.854100000000000276e-02  4.129849999999999910e-01 -2.357567999999999886e-01  9.872999999999999832e-04 -7.600000000000000418e-05  5.899999999999999754e-05
  5.895100000000000000e+04  5.908599999999999963e-02  4.139869999999999939e-01 -2.366149999999999920e-01  7.075000000000000126e-04 -8.000000000000000654e-05  8.600000000000000331e-05
  5.895200000000000000e+04  5.976900000000000268e-02  4.154180000000000095e-01 -2.372105999999999937e-01  4.798000000000000073e-04 -8.399999999999999535e-05  1.129999999999999955e-04
  5.895300000000000000e+04  6.095400000000000124e-02  4.167310000000000181e-01 -2.375994999999999913e-01  3.118999999999999919e-04 -8.700000000000000051e-05  1.399999999999999877e-04
  5.895400000000000000e+04  6.210199999999999748e-02  4.180929999999999924e-01 -2.378588000000000091e-01  1.710000000000000094e-04 -9.100000000000000287e-05  1.669999999999999935e-04
  5.895500000000000000e+04  6.290999999999999370e-02  4.196619999999999795e-01 -2.380454999999999932e-01  1.719000000000000042e-04 -6.000000000000000152e-06  8.100000000000000375e-05
  5.895600000000000000e+04  6.385599999999999610e-02  4.214060000000000028e-01 -2.382557999999999898e-01  2.683000000000000163e-04  1.029999999999999964e-04 -3.600000000000000091e-05
  5.895700000000000000e+04  6.455500000000000127e-02  4.229890000000000039e-01 -2.385857000000000117e-01  4.040000000000000080e-04  1.019999999999999992e-04 -2.000000000000000164e-05
  5.895800000000000000e+04  6.440300000000000191e-02  4.242549999999999932e-01 -2.390210000000000112e-01  4.910999999999999567e-04  5.899999999999999754e-05  4.399999999999999886e-05
\end{verbatim}
)";
#endif

/***********************************************/

#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_EARTHORIENTATIONPARAMETER_TYPE    = "earthOrientationParameter";
constexpr UInt    FILE_EARTHORIENTATIONPARAMETER_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** FUNCTIONS *******************************/

/** @brief Write into a EarthOrientationParameter file.
* each row: mjd(UTC), xp[arcsec], yp[arcsec], dUT1[s], lod[s], dX[arcsec], dY[arcsec]. */
void writeFileEarthOrientationParameter(const FileName &fileName, const const_MatrixSlice &EOP);

/** @brief Read from a EarthOrientationParameter file.
* each row: mjd(UTC), xp[arcsec], yp[arcsec], dUT1[s], lod[s], dX[arcsec], dY[arcsec]. */
void readFileEarthOrientationParameter(const FileName &fileName, Matrix &EOP);

/// @}

/***********************************************/

#endif /* __GROOPS_FILEEARTHORIENTATIONPARAMETER__ */
