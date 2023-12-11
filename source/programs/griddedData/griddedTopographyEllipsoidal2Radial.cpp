/***********************************************/
/**
* @file griddedTopographyEllipsoidal2Radial.cpp
*
* @brief Interpolate digital terrain models from ellipoidal heights to radial heights.
*
* @author Torsten Mayer-Guerr
* @date 2013-10-22
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Interpolate digital terrain models from ellipoidal heights to radial heights.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Interpolate digital terrain models from ellipoidal heights to radial heights.
* @ingroup programsGroup */
class GriddedTopographyEllipsoidal2Radial
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedTopographyEllipsoidal2Radial, SINGLEPROCESS, "Interpolate digital terrain models from ellipoidal heights to radial heights", Grid)

/***********************************************/

void GriddedTopographyEllipsoidal2Radial::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutGrid, fileNameInGrid;

    readConfig(config, "outputfileGriddedData", fileNameOutGrid, Config::OPTIONAL, "", "");
    readConfig(config, "inputfileGriddedData",  fileNameInGrid,  Config::MUSTSET,  "", "Digital Terrain Model");
    if(isCreateSchema(config)) return;

    // read grid
    // ---------
    logStatus<<"read grid from file <"<<fileNameInGrid<<">"<<Log::endl;
    GriddedDataRectangular grid;
    readFileGriddedData(fileNameInGrid, grid);
    MiscGriddedData::printStatistics(grid);

    std::vector<Angle>  lambda, phi;
    std::vector<Double> radius;
    grid.geocentric(lambda, phi, radius);
    const Double direction = ((phi.back() > phi.front()) ? +1 : -1); // north to south or vice versa

    logStatus<<"linear interpolation of elevation data"<<Log::endl;
    Double maxPhiShift    = 0.0;
    Double maxHeightShift = 0.0;
    Single::forEach(grid.longitudes.size(), [&](UInt k)
    {
      Vector heightOld = grid.values.at(0).column(k);
      std::vector<Angle> phiTop(grid.latitudes.size()); // spherical latitude at top of topography
      for(UInt i=0; i<phi.size(); i++)
        phiTop.at(i) = grid.ellipsoid(grid.longitudes.at(k), grid.latitudes.at(i), heightOld(i)).phi();

      // linear interpolation
      for(UInt i=0; i<phi.size(); i++)
      {
        // find interval
        const UInt idx = std::min(i + (direction*(phi.at(i)-phiTop.at(i)) > 0), phiTop.size()-1);
        const Double w = (phi.at(i)-phiTop.at(idx-1)) / (phiTop.at(idx)-phiTop.at(idx-1));
        grid.values.at(0)(i,k) = w * heightOld(idx) + (w-1) * heightOld(idx-1);

        maxPhiShift    = std::max(maxPhiShift,    std::fabs(phiTop.at(i)-phi.at(i)));
        maxHeightShift = std::max(maxHeightShift, std::fabs(grid.values.at(0)(i,k)-heightOld(i)));
      }
    });

    logInfo<<"  max. shifted latitude: "<<maxPhiShift*RAD2DEG<<"°"<<Log::endl;
    logInfo<<"  max. height change:    "<<maxHeightShift<<" m"<<Log::endl;

    // write new grid
    // --------------
    if(!fileNameOutGrid.empty())
    {
      logStatus<<"write grid to file <"<<fileNameOutGrid<<">"<<Log::endl;
      writeFileGriddedData(fileNameOutGrid, grid);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
