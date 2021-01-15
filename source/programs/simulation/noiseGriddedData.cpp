/***********************************************/
/**
* @file noiseGriddedData.cpp
*
* @brief Add noise to gridded data.
*
* @author Torsten Mayer-Guerr
* @date 2020-03-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program adds noise to \file{gridded data data}{griddedData}.
See \configClass{noiseGenerator}{noiseGeneratorType} for details on noise generation.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Add noise to gridded data.
* @ingroup programsGroup */
class NoiseGriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NoiseGriddedData, SINGLEPROCESS, "add noise to gridded data", Simulation, Noise, Grid)

/***********************************************/

void NoiseGriddedData::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName          fileNameOut, fileNameIn;
    NoiseGeneratorPtr noiseGenerator;
    UInt              startData, countData = MAX_UINT;

    readConfig(config, "outputfileGriddedData", fileNameOut,    Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileGriddedData",  fileNameIn,     Config::MUSTSET,  "",  "");
    readConfig(config, "noise",                 noiseGenerator, Config::MUSTSET,  "",  "");
    readConfig(config, "startDataFields",       startData,      Config::DEFAULT,  "0", "start");
    readConfig(config, "countDataFields",       countData,      Config::OPTIONAL, "",  "number of data fields (default: all after start)");
    if(isCreateSchema(config)) return;

    logStatus<<"add noise to data <"<<fileNameIn<<">"<<Log::endl;
    GriddedData grid;
    readFileGriddedData(fileNameIn, grid);
    for(UInt i=startData; i<std::min(grid.values.size(), startData+countData); i++)
    {
      const Vector e = noiseGenerator->noise(grid.points.size(), 1);
      for(UInt k=0; k<grid.points.size(); k++)
        grid.values.at(i).at(k) += e(k);
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
