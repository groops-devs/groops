/***********************************************/
/**
* @file filePolygon.h
*
* @brief Polygons on Earth's surface.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#ifndef __GROOPS_FILEPOLYGON__
#define __GROOPS_FILEPOLYGON__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_Polygon
static const char *docstringPolygon = R"(
List of longitude and latitudes to describe borders, e.g. river basins or continents.
It is used in \configClass{border:polygon}{borderType:polygon}.

\begin{verbatim}
groops polygon version=20200123
          2  # number of polygons
          6  # number of points (1. polygon)
# longitude [deg]           latitude [deg]
# ==================================================
 -1.598200000000000216e+02  2.203000000000000114e+01
 -1.596200000000000045e+02  2.189999999999999858e+01
 -1.593799999999999955e+02  2.189999999999999858e+01
 -1.593000000000000114e+02  2.221999999999999886e+01
 -1.595800000000000125e+02  2.221999999999999886e+01
 -1.598200000000000216e+02  2.203000000000000114e+01
          5  # number of points (2. polygon)
# longitude [deg]           latitude [deg]
# ==================================================
 -7.900000000000000000e+01  2.669999999999999929e+01
 -7.870000000000000284e+01  2.650000000000000000e+01
 -7.823000000000000398e+01  2.667000000000000171e+01
 -7.793000000000000682e+01  2.667000000000000171e+01
 -7.779999999999999716e+01  2.646999999999999886e+01
\end{verbatim}

)";
#endif

/***********************************************/

#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_POLYGON_TYPE = "polygon";

/***** CLASS ***********************************/

/** @brief Polygons on Earth's surface.
* The Polygon is described as clockwise list of points (latitude, longitude). */
class Polygon
{
public:
  Vector L,B;

  Polygon() {}
  Polygon(UInt size) : L(size), B(size) {}
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const std::vector<Polygon> &x);
template<> void load(InArchive  &ar, std::vector<Polygon> &x);

/** @brief Write into a Polygon file. */
void writeFilePolygon(const FileName &fileName, const std::vector<Polygon> &x);

/** @brief Read from a Polygon file. */
void readFilePolygon(const FileName &fileName, std::vector<Polygon> &x);

/// @}

/***********************************************/

#endif /* __GROOPS_FILEPOLYGON__ */
