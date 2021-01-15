/***********************************************/
/**
* @file griddedDataInterpolate.cpp
*
* @brief Interpolate values of rectangular grids to new points.
*
* @author Torsten Mayer-Guerr
* @date 2019-05-11
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Interpolate values of a regular rectangular \configFile{inputfileGriddedData}{griddedData}
to new points given by \configClass{grid}{gridType} and write as \configFile{outputfileGriddedData}{griddedData}.
Only longitude and latitude of points are considered; the height is ignored for interpolation.

(Only nearest neighbor method is implemented at the moment.)

\fig{!hb}{0.8}{griddedDataInterpolate}{fig:griddedDataInterpolate}{Interpolation of point data from rectangular gridded data.}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Interpolate values of rectangular grids to new points.
* @ingroup programsGroup */
class GriddedDataInterpolate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedDataInterpolate, SINGLEPROCESS, "Interpolate values of rectangular grids to new points", Grid)

/***********************************************/

void GriddedDataInterpolate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameOutGrid, fileNameInGrid;
    GridPtr     gridPtr;
    std::string choice;

    readConfig(config, "outputfileGriddedData", fileNameOutGrid, Config::MUSTSET, "", "");
    readConfig(config, "inputfileGriddedData",  fileNameInGrid,  Config::MUSTSET, "", "must be rectangular");
    readConfig(config, "grid",                  gridPtr,         Config::MUSTSET, "", "");
    if(readConfigChoice(config, "method", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "nearestNeighbor", choice, "")) {};
      endChoice(config);
    }
    if(isCreateSchema(config)) return;

    // read grid
    // ---------
    logStatus<<"read grid from file <"<<fileNameInGrid<<">"<<Log::endl;
    GriddedDataRectangular grid;
    readFileGriddedData(fileNameInGrid, grid);
    MiscGriddedData::printStatistics(grid);

    logStatus<<"create grid"<<Log::endl;
    GriddedData pointList(grid.ellipsoid, gridPtr->points(), gridPtr->areas(), std::vector<std::vector<Double>>(grid.values.size(), std::vector<Double>(gridPtr->points().size(), 0.)));

    // interpolate
    // -----------
    for(UInt i=0; i<pointList.points.size(); i++)
    {
      Angle  L, B;
      Double h;
      grid.ellipsoid(pointList.points.at(i), L, B, h);

      // find nearest neighbor
      const UInt row = std::distance(grid.latitudes.begin(), std::min_element(grid.latitudes.begin(), grid.latitudes.end(),
                                     [&B](Angle B1, Angle B2) {return std::fabs(B-B1) < std::fabs(B-B2);}));
      const UInt col = std::distance(grid.longitudes.begin(), std::min_element(grid.longitudes.begin(), grid.longitudes.end(),
                                     [&L](Angle L1, Angle L2) {return std::fabs(L-L1) < std::fabs(L-L2);}));

      for(UInt k=0; k<grid.values.size(); k++)
        pointList.values.at(k).at(i) = grid.values.at(k)(row, col);
    }

    // write
    // -----
    logStatus<<"write interpolated values to file <"<<fileNameOutGrid<<">"<<Log::endl;
    writeFileGriddedData(fileNameOutGrid, pointList);
    MiscGriddedData::printStatistics(pointList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
