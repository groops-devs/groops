/***********************************************/
/**
* @file fileGriddedData.h
*
* @brief Read/write gridded values.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#ifndef __GROOPS_FILEGRIDDEDDATA__
#define __GROOPS_FILEGRIDDEDDATA__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_GriddedData
static const char *docstringGriddedData = R"(
List of arbitrarily distributed points defined by geographic coordinates and ellipsoidal
height. Each point can also have an associated area
(projected on the unit sphere with a total area of $4\pi$).
This file format supports multiple values per point (called \verb|data0|, \verb|data1| and so on).

For regular gridded data and binary format (\verb|*.dat|) a more efficient storage scheme is used.

See also: \program{GriddedDataCreate}.

\begin{verbatim}
groops griddedData version=20200123
 1  2  6.378137000000000000e+06  6.356752314140356146e+06 72 # hasArea, data columns, ellipoid a, ellipoid b, data rows
# longitude [deg]           latitude [deg]            height [m]                unit areas [-]             data0                     data1
# ===========================================================================================================================================================
 -1.650000000000000000e+02  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
 -1.350000000000000000e+02  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
 -1.050000000000000142e+02  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
 -7.500000000000001421e+01  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
 -4.500000000000002132e+01  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
 -1.500000000000002132e+01  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
  1.499999999999997691e+01  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
  4.499999999999997868e+01  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
  7.499999999999997158e+01  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
  1.049999999999999574e+02  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
  1.349999999999999432e+02  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
  1.649999999999999432e+02  7.500000000000000000e+01  0.000000000000000000e+00  7.014893453974438420e-02  1.000000000000000000e+00  2.000000000000000000e+00
 -1.650000000000000000e+02  4.500000000000000711e+01  0.000000000000000000e+00  1.916504532594049681e-01  1.000000000000000000e+00  2.000000000000000000e+00
 -1.350000000000000000e+02  4.500000000000000711e+01  0.000000000000000000e+00  1.916504532594049681e-01  1.000000000000000000e+00  2.000000000000000000e+00
\end{verbatim}

)";
#endif

/***********************************************/

#include "base/exception.h"
#include "base/griddedData.h"
#include "inputOutput/fileName.h"
#include "inputOutput/archive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_GRIDDEDDATA_TYPE = "griddedData";

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const GriddedData &x);
template<> void save(OutArchive &ar, const GriddedDataRectangular &x);
template<> void load(InArchive  &ar, GriddedData &x);
template<> void load(InArchive  &ar, GriddedDataRectangular &x);

/** @brief Write into a GriddedData file. */
void writeFileGriddedData(const FileName &fileName, const GriddedData &x);
void writeFileGriddedData(const FileName &fileName, const GriddedDataRectangular &x);

/** @brief Read from a GriddedData file. */
void readFileGriddedData(const FileName &fileName, GriddedData &x);
void readFileGriddedData(const FileName &fileName, GriddedDataRectangular &x);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
