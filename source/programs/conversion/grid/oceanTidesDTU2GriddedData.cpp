/***********************************************/
/**
* @file oceanTidesDTU2GriddedData.cpp
*
* @brief Convert DTU ocean tide grids to griddedData (amplitude, phase).
*
* @author Torsten Mayer-Guerr
* @date 2024-11-09
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert DTU ocean tide grids to griddedData (amplitude, phase).
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Convert DTU ocean tide grids to griddedData (amplitude, phase).
* @ingroup programsConversionGroup */
class OceanTidesDTU2GriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(OceanTidesDTU2GriddedData, SINGLEPROCESS, "Convert DTU ocean tide grids to griddedData (amplitude, phase)", Conversion)

/***********************************************/

void OceanTidesDTU2GriddedData::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameIn;
    Double   a, f;

    readConfig(config, "outputfileGriddedData", fileNameOut, Config::MUSTSET, "", "data0=amplitude, data1=phase");
    readConfig(config, "inputfile",             fileNameIn,  Config::MUSTSET, "", "");
    readConfig(config, "R",                     a,           Config::DEFAULT, "6378136.3", "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",     f,           Config::DEFAULT, "298.257",   "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    // =====================================================

    // Create rectangular grid
    // -----------------------
    GriddedDataRectangular grid;
    grid.ellipsoid = Ellipsoid(a,f);

    logStatus<<"read file <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);
    for(UInt idx=0; idx<2; idx++) // amplitude, phase
    {
      // M2      DTU23 ocean tide constituent
      //  Phase   (degrees)
      //       2881      5761
      //    -90.000    90.000
      //      0.000   360.000
      //    999.000   999.000
      // (11f7.2)
      UInt   rows, columns;
      Double minLat, maxLat, minLon, maxLon;
      Double noValue;
      std::string line;
      getline(file, line); // constituent
      logInfo<<"'"<<line<<"'"<<Log::endl;
      getline(file, line); // Amplitude or phase
      file>>rows>>columns;
      file>>minLat>>maxLat;
      file>>minLon>>maxLon;
      file>>noValue>>noValue;
      getline(file, line);
      getline(file, line); // format

      grid.longitudes.resize(columns-1);
      grid.latitudes.resize(rows);
      grid.heights.resize(rows);
      grid.values.resize(2, Matrix(rows, columns-1));
      for(UInt k=0; k<columns; k++)
        grid.longitudes.at(k%(columns-1)) = Angle(DEG2RAD*(minLon + k*(maxLon-minLon)/(columns-1)));
      for(UInt i=0; i<rows; i++)
        grid.latitudes.at(i) = Angle(DEG2RAD*(minLat + i*(maxLat-minLat)/(rows-1)));

      for(UInt i=0; i<rows; i++)
        for(UInt k=0; k<columns; k++)
        {
          file>>grid.values.at(idx)(i,k%(columns-1));
          if(grid.values.at(idx)(i,k%(columns-1)) == noValue)
            grid.values.at(idx)(i,k%(columns-1)) = 0.;
        }
      getline(file, line);
    }

    logStatus<<"write gridded data to file <"<<fileNameOut<<">"<<Log::endl;
    writeFileGriddedData(fileNameOut, grid);
    MiscGriddedData::printStatistics(grid);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
