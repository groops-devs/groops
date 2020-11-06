/***********************************************/
/**
* @file fileGriddedDataTimeSeries.h
*
* @brief Read/write time series of gridded values.
*
* @author Torsten Mayer-Guerr
* @date 2020-07-28
*
*/
/***********************************************/

#ifndef __GROOPS_FILEGRIDDEDDATATIMESERIES__
#define __GROOPS_FILEGRIDDEDDATATIMESERIES__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_GriddedDataTimeSeries
static const char *docstringGriddedDataTimeSeries = R"(
Time series of data for arbitrarily distributed points defined by geographic coordinates and ellipsoidal
height. The data can be temporal interpolated by \reference{basis splines}{fundamentals.basisSplines}.
The file format consists of a \file{griddedData}{griddedData}, a time series, and
for each spatial point and spline node pair multiple values called \verb|data0|, \verb|data1|, \ldots.
)";
#endif

/***********************************************/

#include "base/exception.h"
#include "base/griddedData.h"
#include "inputOutput/fileName.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_GRIDDEDDATATIMESERIES_TYPE = "griddedDataTimeSeries";

/***** CLASS ***********************************/

/** @brief Time series of gridded values. */
class InFileGriddedDataTimeSeries
{
  InFileArchive       file;
  UInt                splineDegree_, dataCount_;
  GriddedData         grid_;
  std::vector<Time>   times_;
  UInt                indexFile, indexData;
  std::streampos      seekStart;
  std::streamoff      seekSize;
  std::vector<Matrix> data_;

public:
  InFileGriddedDataTimeSeries() : splineDegree_(0), dataCount_(0) {}
  InFileGriddedDataTimeSeries(const FileName &name) {open(name);}
 ~InFileGriddedDataTimeSeries() {close();}

  void open(const FileName &name);
  void close();

  UInt splineDegree() const {return splineDegree_;}
  UInt nodeCount()    const {return times_.size()+splineDegree_;}  //!< number of spline nodal points (agree with times().size() for spline degree 1)
  UInt dataCount()    const {return dataCount_;}                   //!< number of data columns

  const GriddedData       &grid()  const {return grid_;}
  const std::vector<Time> &times() const {return times_;}

  /// data(points x data columns) at spline nodal point
  Matrix data(UInt idNode);

  /// data(points x data columns) interpolated at @p time.
  Matrix data(const Time &time);
};

/***** FUNCTIONS *******************************/

/** @brief Write into a GriddedData time series file. */
void writeFileGriddedDataTimeSeries(const FileName &fileName, UInt splineDegree, const std::vector<Time> &times,
                                    const GriddedData &grid, const std::vector<Matrix> &data);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
