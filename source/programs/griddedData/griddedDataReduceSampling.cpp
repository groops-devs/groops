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
Generate coarse grid by computing mean values.
The number of points is decimated by averaging integer multiplies of grid points
(\config{multiplierLongitude}, \config{multiplierLatitude}).
The fine grid can be written, where the coarse grid values are additionally appended.
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

    readConfig(config, "outputfileCoarseGridRectangular", fileNameOutCoarseGrid, Config::OPTIONAL, "", "coarse grid");
    readConfig(config, "outputfileFineGridRectangular",   fileNameOutFineGrid,   Config::OPTIONAL, "", "fine grid with additional coarse grid values");
    readConfig(config, "inputfileFineGridRectangular",    fileNameInFineGrid,    Config::MUSTSET,  "", "Digital Terrain Model");
    readConfig(config, "multiplierLongitude",             numberCols,            Config::MUSTSET,  "8", "Generalizing factor");
    readConfig(config, "multiplierLatitude",              numberRows,            Config::MUSTSET,  "8", "Generalizing factor");
    if(isCreateSchema(config)) return;

    // read grid
    // ---------
    logStatus<<"read grid from file <"<<fileNameInFineGrid<<">"<<Log::endl;
    GriddedDataRectangular grid1;
    readFileGriddedData(fileNameInFineGrid, grid1);
    MiscGriddedData::printStatistics(grid1);

    std::vector<Double> radius1, dLambda1, dPhi1;
    std::vector<Angle>  lambda1;
    std::vector<Angle>  phi1;
    grid1.geocentric(lambda1, phi1, radius1, dLambda1, dPhi1);
    const UInt rows1 = phi1.size();
    const UInt cols1 = lambda1.size();

   // if((rows1%numberRows!=0)||(cols1%numberCols!=0))
   //   throw(Exception("rows or columns cannot divided by multiplier"));

    // Generate coarse grid
    // --------------------
    logStatus<<"generate coarse grid"<<Log::endl;
    const UInt rows2 = rows1/numberRows;
    const UInt cols2 = cols1/numberCols;
    GriddedDataRectangular grid2;
    grid2.ellipsoid = grid1.ellipsoid;
    grid2.longitudes.resize(cols2);
    grid2.latitudes.resize(rows2);
    grid2.heights.resize(rows2);
    for(UInt idx=0; idx<grid1.values.size(); idx++)
      grid2.values.push_back( Matrix(rows2, cols2) );

    for(UInt s=0; s<cols2; s++)
    {
      Double sum    = 0;
      Double weight = 0;
      for(UInt k=0; k<numberCols; k++)
      {
        sum    += dLambda1.at(s*numberCols+k) * grid1.longitudes.at(s*numberCols+k);
        weight += dLambda1.at(s*numberCols+k);
      }
      grid2.longitudes.at(s) = Angle(sum/weight);
    }

    for(UInt z=0; z<rows2; z++)
    {
      Double sum    = 0;
      Double weight = 0;
      for(UInt i=0; i<numberRows; i++)
      {
        sum    += dPhi1.at(z*numberRows+i) * grid1.latitudes.at(z*numberRows+i);
        weight += dPhi1.at(z*numberRows+i);
      }
      grid2.latitudes.at(z) = Angle(sum/weight);
    }

    std::vector<Double> radius2, dLambda2, dPhi2;
    std::vector<Angle>  lambda2;
    std::vector<Angle>  phi2;
    grid2.geocentric(lambda2, phi2, radius2, dLambda2, dPhi2);

    // Compute mean values
    // -------------------
    logStatus<<"compute mean values"<<Log::endl;
    grid1.values.push_back( Matrix(rows1, cols1) );
    for(UInt idx=0; idx<grid2.values.size(); idx++)
      for(UInt z=0; z<rows2; z++)
        for(UInt s=0; s<cols2; s++)
        {
          Double sum    = 0;
          Double weight = 0;
          for(UInt i=0; i<numberRows; i++)
          {
            const Double dPhicosPhi = cos(phi1.at(z*numberRows+i)) * dPhi1.at(z*numberRows+i);
            for(UInt k=0; k<numberCols; k++)
            {
              const Double w = dPhicosPhi * dLambda1.at(s*numberCols+k);
              sum    += w * grid1.values.at(idx)(z*numberRows+i, s*numberCols+k);
              weight += w;
            }
          }
          sum /= weight;
          grid2.values.at(idx)(z, s) = sum;

          if(idx==0)
            for(UInt i=0; i<numberRows; i++)
              for(UInt k=0; k<numberCols; k++)
                grid1.values.back()(z*numberRows+i, s*numberCols+k) = sum;
        }

    MiscGriddedData::printStatistics(grid1);
    MiscGriddedData::printStatistics(grid2);

    // write new grid
    // --------------
    if(!fileNameOutCoarseGrid.empty())
    {
      logStatus<<"write coarse grid to file <"<<fileNameOutCoarseGrid<<">"<<Log::endl;
      writeFileGriddedData(fileNameOutCoarseGrid, grid2);
    }

    if(!fileNameOutFineGrid.empty())
    {
      logStatus<<"write fine grid to file <"<<fileNameOutFineGrid<<">"<<Log::endl;
      writeFileGriddedData(fileNameOutFineGrid, grid1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
