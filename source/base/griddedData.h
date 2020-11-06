/***********************************************/
/**
* @file griddedData.h
*
* @brief Gridded values.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDDEDDATA__
#define __GROOPS_GRIDDEDDATA__

#include "base/importStd.h"
#include "base/vector3d.h"
#include "base/ellipsoid.h"

/** @addtogroup base */
/// @{

/***** CLASS ***********************************/

/** @brief Point list with (multiple) data. */
class GriddedData
{
public:
  Ellipsoid                        ellipsoid; //!< Ellipsoid for conversion to ellipsoidal coordinates
  std::vector<Vector3d>            points;    //!< List of points
  std::vector<Double>              areas;     //!< Area element (projected to unit sphere).
  std::vector<std::vector<Double>> values;    //!< data.at(dataIdx).at(pointIdx)

  /// Default constructor.
  GriddedData() = default;

  /// Constructor with points, area elements, and mutiple values for each point.
  GriddedData(const Ellipsoid &ellip, const std::vector<Vector3d> &_points, const std::vector<Double> &_areas, const std::vector<std::vector<Double>> &_values) : ellipsoid(ellip), points(_points), areas(_areas), values(_values) {}

  /** @brief Sort points geographically (North/West->South/East). */
  void sort();

  /** @brief Define points a rectangular grid?.
  * if function returns FALSE, @a lambda, @a phi, @a radius contain garbage.
  * if function returns TRUE, points are in same order as:
  @code
  UInt idx = 0;
  for(UInt i=0; i<phi.size(); i++)
    for(UInt k=0; k<lambda.size(); k++)
      points.at(idx++) == polar(lambda.at(k), phi.at(i), r.at(i));
  @endcode */
  Bool isRectangle(std::vector<Angle> &lambda, std::vector<Angle> &phi, std::vector<Double> &radius) const;

  /** @brief Automatically area computation of rectangular grids (overwrite areas). */
  Bool computeArea();

  /** @brief Is GriddedData valid?.
  * Test dimensions of vectors. */
  Bool isValid() const;
};

/***** CLASS ***********************************/

/** @brief Rectangular grid with (multiple) data. */
class GriddedDataRectangular
{
public:
  Ellipsoid           ellipsoid;  //!< Ellipsoid for conversion to ellipsoidal coordinates
  std::vector<Angle>  longitudes; //!< Longitude (columns)
  std::vector<Angle>  latitudes;  //!< Latitude (rows)
  std::vector<Double> heights;    //!< Elliposoidal height (rows)
  std::vector<Matrix> values;     //!< Multiple values at each point.

  void geocentric(std::vector<Angle> &lambda, std::vector<Angle> &phi, std::vector<Double> &radius,
                  std::vector<Double> &dLambda, std::vector<Double> &dPhi) const;

  /// Convert to GriddedData.
  void convert(GriddedData &grid) const;

  /// Create from GriddedData.
  Bool init(const GriddedData &grid);

  /** @brief Is GriddedDataRectangular valid?.
  * Test dimensions of vectors. */
  Bool isValid() const;
};

/// @}

/***********************************************/

#endif /* __GROOPS__ */
