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
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GriddedTopographyEllipsoidal2Radial, PARALLEL, "Interpolate digital terrain models from ellipoidal heights to radial heights", Grid)

/***********************************************/

void GriddedTopographyEllipsoidal2Radial::run(Config &config)
{
  try
  {
    if(!Parallel::isMaster())
      return;

    FileName fileNameOutGrid, fileNameInGrid;

    readConfig(config, "outputfileGriddedData", fileNameOutGrid, Config::MUSTSET, "", "");
    readConfig(config, "inputfileGriddedData",  fileNameInGrid,  Config::MUSTSET, "", "Digital Terrain Model");
    if(isCreateSchema(config)) return;

    // read grid
    // ---------
    logStatus<<"read grid from file <"<<fileNameInGrid<<">"<<Log::endl;
    GriddedDataRectangular grid;
    readFileGriddedData(fileNameInGrid, grid);
    MiscGriddedData::printStatistics(grid);

    std::vector<Double> radius, dLambda, dPhi;
    std::vector<Angle>  lambda;
    std::vector<Angle>  phi;
    grid.geocentric(lambda, phi, radius, dLambda, dPhi);
    const UInt rows = phi.size();
    const UInt cols = lambda.size();

    // Compute origonal top position
    // -----------------------------
    logStatus<<"compute original top position"<<Log::endl;
    const Double direction = ((phi.at(1)>phi.at(0)) ? +1 : -1); // north to south or vice versa
    Matrix heightOld = grid.values.at(0);
    Matrix phiShift(rows,cols); // difference of spherical latitude at top of topography
    Double maxPhiShift = 0.0;
    Double maxHeightShift = 0.0;
    Double maxWeight = 0.0;

    for(UInt z=0; z<rows; z++)
      for(UInt s=0; s<cols; s++)
      {
        // x,y,z of top of topography
        Vector3d pointTop = grid.ellipsoid(grid.longitudes.at(s), grid.latitudes.at(z), heightOld(z,s));
        phiShift(z,s) = direction * (static_cast<Double>(pointTop.phi()) - phi.at(z)); // >0: point is near z+1

        maxPhiShift    = std::max(maxPhiShift,    fabs(phiShift(z,s)));
        maxHeightShift = std::max(maxHeightShift, radius.at(z)*fabs(phiShift(z,s)));
        maxWeight      = std::max(maxWeight,      fabs(phiShift(z,s)/dPhi.at(z)));
      }
    logInfo<<"  max. shifted point:    "<<maxHeightShift<<" m"<<Log::endl;
    logInfo<<"  max. shifted latitude: "<<maxPhiShift*RAD2DEG<<"°"<<Log::endl;
    logInfo<<"  max. weight:           "<<maxWeight*100<<"%"<<Log::endl;

    // interpolate data
    // ----------------
    logStatus<<"interpolate data"<<Log::endl;
    for(UInt z=0; z<rows; z++)
      for(UInt s=0; s<cols; s++)
      {
        Double w = phiShift(z,s)/dPhi.at(z);
        if((w<0) && (z+1<rows))
          grid.values.at(0)(z,s) = (1+w) * heightOld(z,s) - w * heightOld(z+1,s);
        else if((w>0) && (z>0))
          grid.values.at(0)(z,s) = (1-w) * heightOld(z,s) + w * heightOld(z-1,s);
      }

    MiscGriddedData::printStatistics(grid);

    // write new grid
    // --------------
    logStatus<<"write grid to file <"<<fileNameOutGrid<<">"<<Log::endl;
    writeFileGriddedData(fileNameOutGrid, grid);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
