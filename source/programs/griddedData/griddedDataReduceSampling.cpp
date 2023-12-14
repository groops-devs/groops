/***********************************************/
/**
* @file griddedDataReduceSampling.cpp
*
* @brief Generate coarse grid by computing mean values.
*
* @author Torsten Mayer-Guerr
* @date 2013-10-24
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Generate coarse grid by computing area weighted mean values.
The number of points is decimated by averaging integer multiplies of grid points
(\config{multiplierLongitude}, \config{multiplierLatitude}).

if \config{volumeConserving} is set, data are interpreted as heights above ellipsoid
and the tesseroid volume
\begin{equation}
  V=\int_r^{r+H}\int_{\varphi_1}^{\varphi_2}\int_{\lambda_1}^{\lambda_2} r^2\cos\varphi\,d\varphi\,d\lambda\,dr
\end{equation}
is conserved, where $r$ is the radius of the ellipsoid at grid center and
$(\varphi_1-\varphi_2)\times(\lambda_1-\lambda_2)$ are the grid cell boundaries.
This is meaninful for Digital Elevation Models (DEM).

The fine grid can be written, where the first coarse grid values (data0) are additionally appended.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Generate coarse grid by computing mean values.
* @ingroup programsGroup */
class GriddedDataReduceSampling
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedDataReduceSampling, SINGLEPROCESS, "Generate coarse grid by computing mean values", Grid)

/***********************************************/

void GriddedDataReduceSampling::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutCoarseGrid, fileNameOutFineGrid;
    FileName fileNameInFineGrid;
    UInt     numberRows, numberCols;
    Bool     volumeConserving;

    readConfig(config, "outputfileCoarseGridRectangular", fileNameOutCoarseGrid, Config::OPTIONAL, "", "coarse grid");
    readConfig(config, "outputfileFineGridRectangular",   fileNameOutFineGrid,   Config::OPTIONAL, "", "fine grid with additional coarse grid values");
    readConfig(config, "inputfileFineGridRectangular",    fileNameInFineGrid,    Config::MUSTSET,  "", "Digital Terrain Model");
    readConfig(config, "multiplierLongitude",             numberCols,            Config::MUSTSET,  "8", "Generalizing factor");
    readConfig(config, "multiplierLatitude",              numberRows,            Config::MUSTSET,  "8", "Generalizing factor");
    readConfig(config, "volumeConserving",                volumeConserving,      Config::DEFAULT,  "0", "data are interpreted as heights above ellipsoid");
    if(isCreateSchema(config)) return;

    // read grid
    // ---------
    logStatus<<"read grid from file <"<<fileNameInFineGrid<<">"<<Log::endl;
    GriddedDataRectangular gridFine;
    readFileGriddedData(fileNameInFineGrid, gridFine);

    std::vector<Angle>  tmp;
    std::vector<Double> radiusFine, longitude, latitude, dLambda, dPhi;
    gridFine.geocentric(tmp, tmp, radiusFine);
    gridFine.cellBorders(longitude, latitude);
    gridFine.areaElements(dLambda, dPhi);

    // Generate coarse grid
    // --------------------
    logStatus<<"generate coarse grid"<<Log::endl;
    GriddedDataRectangular gridCoarse;
    gridCoarse.ellipsoid = gridFine.ellipsoid;

    gridCoarse.latitudes.resize(gridFine.latitudes.size()/numberRows);
    for(UInt i=0; i<gridCoarse.latitudes.size(); i++)
      gridCoarse.latitudes.at(i) = 0.5*(latitude.at(i*numberRows) + latitude.at((i+1)*numberRows));

    gridCoarse.longitudes.resize(gridFine.longitudes.size()/numberCols);
    for(UInt k=0; k<gridCoarse.longitudes.size(); k++)
      gridCoarse.longitudes.at(k) = std::remainder(longitude.at(k*numberCols) + 0.5*std::remainder(longitude.at((k+1)*numberCols)-longitude.at(k*numberCols), 2*PI), 2*PI);

    gridCoarse.heights.resize(gridCoarse.latitudes.size()); // Compute mean height
    for(UInt i=0; i<gridCoarse.heights.size(); i++)
    {
      Double height = 0, weight = 0;
      for(UInt ii=0; ii<numberRows; ii++)
      {
        height += dPhi.at(i) * gridFine.heights.at(i*numberRows+ii);
        weight += dPhi.at(i);
      }
      gridCoarse.heights.at(i) = height/weight;
    }

    gridCoarse.values.resize(gridFine.values.size(), Matrix(gridCoarse.latitudes.size(), gridCoarse.longitudes.size()));

    std::vector<Double> radiusCoarse, dLambdaCoarse, dPhiCoarse;
    gridCoarse.geocentric(tmp, tmp, radiusCoarse);
    gridCoarse.areaElements(dLambdaCoarse, dPhiCoarse);

    // Compute mean values
    // -------------------
    logStatus<<"compute mean values"<<Log::endl;
    if(!fileNameOutFineGrid.empty())
      gridFine.values.push_back(Matrix(gridFine.latitudes.size(), gridFine.longitudes.size()));
    Single::forEach(gridCoarse.latitudes.size(), [&](UInt i)
    {
      for(UInt k=0; k<gridCoarse.longitudes.size(); k++)
      {
        // compute volume of tesseroid V = (r2^3-r1^3)/3 * area with r2=r1+h
        Vector volumes(gridCoarse.values.size());
        for(UInt ii=0; ii<numberRows; ii++)
          for(UInt kk=0; kk<numberCols; kk++)
          {
            const Double area = dPhi.at(i*numberRows+ii) * dLambda.at(k*numberCols+kk);
            const Double r1   = radiusFine.at(i*numberRows+ii);
            for(UInt idx=0; idx<gridCoarse.values.size(); idx++)
              if(volumeConserving)
                volumes(idx) += area * (std::pow(r1+gridFine.values.at(idx)(i*numberRows+ii, k*numberCols+kk), 3) - std::pow(r1, 3))/3.;
            else
              volumes(idx) += area * gridFine.values.at(idx)(i*numberRows+ii, k*numberCols+kk);
          }

        // height of tesseroid = (V*3/area - r1^3)^(1/3) - r1
        const Double area = dPhiCoarse.at(i) * dLambdaCoarse.at(k);
        const Double r1   = radiusCoarse.at(i);
        for(UInt idx=0; idx<gridCoarse.values.size(); idx++)
          if(volumeConserving)
            gridCoarse.values.at(idx)(i, k) = std::pow(volumes(idx)*3/area + std::pow(r1, 3), 1./3.) - r1;
          else
            gridCoarse.values.at(idx)(i, k) = volumes(idx)/area;

        if(!fileNameOutFineGrid.empty())
          for(UInt ii=0; ii<numberRows; ii++)
            for(UInt kk=0; kk<numberCols; kk++)
              gridFine.values.back()(i*numberRows+ii, k*numberCols+kk) = gridCoarse.values.at(0)(i, k);
      }
    });

    // write new grid
    // --------------
    if(!fileNameOutCoarseGrid.empty())
    {
      logStatus<<"write coarse grid to file <"<<fileNameOutCoarseGrid<<">"<<Log::endl;
      writeFileGriddedData(fileNameOutCoarseGrid, gridCoarse);
      MiscGriddedData::printStatistics(gridCoarse);
    }

    if(!fileNameOutFineGrid.empty())
    {
      logStatus<<"write fine grid to file <"<<fileNameOutFineGrid<<">"<<Log::endl;
      writeFileGriddedData(fileNameOutFineGrid, gridFine);
      MiscGriddedData::printStatistics(gridFine);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
