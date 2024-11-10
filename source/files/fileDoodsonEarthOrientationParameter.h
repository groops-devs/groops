/***********************************************/
/**
* @file fileDoodsonEarthOrientationParameter.h
*
* @brief Read/write high frequent EOPs.
*
* @author Torsten Mayer-Guerr
* @date 2019-05-15
*
*/
/***********************************************/

#ifndef __GROOPS_FILEDOODSONEARTHORIENTATIONPARAMETER__
#define __GROOPS_FILEDOODSONEARTHORIENTATIONPARAMETER__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_DoodsonEarthOrientationParameter
static const char *docstringDoodsonEarthOrientationParameter = R"(
Corrections for Earth orientation parameters (EOP) ($x_p, y_p, UT1, LOD$)
as cos/sin oscillations for a list of doodson tidal frequencies.

\begin{verbatim}
groops doodsonEarthOrientationParameter version=20200123
        11 # number of constituents
# dood.   xpCos [arcsec]            xpSin [arcsec]            ypCos [arcsec]            ypSin [arcsec]            ut1Cos [sec]              ut1Sin [sec]              lodCos [sec]              lodSin [sec]             name
# ===========================================================================================================================================================================================================================
 155.645  2.399999999999999786e-07  1.799999999999999971e-07  1.799999999999999971e-07 -2.399999999999999786e-07 -2.000000000000000042e-08  2.000000000000000042e-08 -1.499999999999999932e-07 -1.400000000000000095e-07 ""
 155.655 -8.220000000000000920e-06 -6.280000000000000012e-06 -6.280000000000000012e-06  8.220000000000000920e-06  7.899999999999999537e-07 -8.599999999999999187e-07  5.219999999999999150e-06  4.829999999999999496e-06 m1
 155.665 -1.649999999999999872e-06 -1.260000000000000007e-06 -1.260000000000000007e-06  1.649999999999999872e-06  1.600000000000000033e-07 -1.700000000000000135e-07  1.049999999999999900e-06  9.700000000000000302e-07 ""
 157.455 -1.539999999999999867e-06 -1.199999999999999946e-06 -1.199999999999999946e-06  1.539999999999999867e-06  1.499999999999999932e-07 -1.600000000000000033e-07  9.799999999999999345e-07  8.899999999999999491e-07 chi1
 157.465 -3.400000000000000270e-07 -2.599999999999999988e-07 -2.599999999999999988e-07  3.400000000000000270e-07  2.999999999999999732e-08 -4.000000000000000084e-08  2.099999999999999746e-07  1.999999999999999909e-07 ""
 161.557  9.999999999999999547e-08  8.000000000000000167e-08  8.000000000000000167e-08 -9.999999999999999547e-08 -1.000000000000000021e-08  1.000000000000000021e-08 -7.000000000000000477e-08 -4.999999999999999774e-08 ""
 162.556  2.549999999999999726e-06  2.019999999999999718e-06  2.019999999999999718e-06 -2.549999999999999726e-06 -2.099999999999999746e-07  2.899999999999999763e-07 -1.799999999999999919e-06 -1.289999999999999931e-06 pi1
 163.545 -4.899999999999999672e-07 -3.799999999999999616e-07 -3.799999999999999616e-07  4.899999999999999672e-07  4.000000000000000084e-08 -5.999999999999999464e-08  3.499999999999999842e-07  2.399999999999999786e-07 ""
 163.555  4.272999999999999238e-05  3.010999999999999801e-05  3.010999999999999801e-05 -4.272999999999999238e-05 -3.079999999999999734e-06  5.219999999999999150e-06 -3.270999999999999684e-05 -1.929999999999999817e-05 p1
 164.554 -3.599999999999999943e-07 -2.800000000000000191e-07 -2.800000000000000191e-07  3.599999999999999943e-07  2.999999999999999732e-08 -4.000000000000000084e-08  2.700000000000000090e-07  1.700000000000000135e-07 ""
 164.556 -1.029999999999999879e-06 -7.999999999999999638e-07 -7.999999999999999638e-07  1.029999999999999879e-06  8.000000000000000167e-08 -1.199999999999999893e-07  7.599999999999999233e-07  4.899999999999999672e-07 ""
\end{verbatim}
)";
#endif

/***********************************************/

#include "base/exception.h"
#include "base/doodson.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_DOODSONEARTHORIENTATIONPARAMETER_TYPE    = "doodsonEarthOrientationParameter";
constexpr UInt    FILE_DOODSONEARTHORIENTATIONPARAMETER_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** CLASS ***********************************/

/** @brief Point list with (multiple) data. */
class DoodsonEop
{
public:
  std::vector<Doodson> doodson;
  Matrix coeff; // for each consituent (row), columns: xpCos,  xpSin, ypCos, ypSin [arcsec], ut1Cos, ut1Sin, lodCos, lodSin [sec];
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const DoodsonEop &x);
template<> void load(InArchive  &ar, DoodsonEop &x);

/** @brief Write into a DoodsonEarthOrientationParameter file. */
void writeFileDoodsonEarthOrientationParameter(const FileName &fileName, const DoodsonEop &x);

/** @brief Read from a DoodsonEarthOrientationParameter file. */
void readFileDoodsonEarthOrientationParameter(const FileName &fileName, DoodsonEop &x);

/// @}

/***********************************************/

#endif /* __GROOPS_FILEDOODSONEARTHORIENTATIONPARAMETER__ */
