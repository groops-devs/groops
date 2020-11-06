/***********************************************/
/**
* @file noiseInstrument.cpp
*
* @brief Add noise to instrument data.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-02
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program adds noise to \file{instrument data}{instrument}.
See \configClass{noiseGenerator}{noiseGeneratorType} for details on noise generation.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/***** CLASS ***********************************/

/** @brief Add noise to instrument data.
* @ingroup programsGroup */
class NoiseInstrument
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(NoiseInstrument, PARALLEL, "add noise to instrument data", Simulation, Noise, Instrument)

/***********************************************/

void NoiseInstrument::run(Config &config)
{
  try
  {
    FileName          fileNameOut, fileNameIn;
    NoiseGeneratorPtr noiseGenerator;
    UInt              startData, countData = MAX_UINT;

    readConfig(config, "outputfileInstrument", fileNameOut,    Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileInstrument",  fileNameIn,     Config::MUSTSET,  "",  "");
    readConfig(config, "noise",                noiseGenerator, Config::MUSTSET,  "",  "");
    readConfig(config, "startDataFields",      startData,      Config::DEFAULT,  "0", "start");
    readConfig(config, "countDataFields",      countData,      Config::OPTIONAL, "",  "number of data fields (default: all after start)");
    if(isCreateSchema(config)) return;

    logStatus<<"add noise to data <"<<fileNameIn<<">"<<Log::endl;
    InstrumentFile instrumentFile(fileNameIn);

    std::vector<Arc> arcList(instrumentFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      Arc    arc  = instrumentFile.readArc(arcNo);
      Matrix data = arc.matrix();
      countData   = std::min(countData, data.columns()-1-startData);
      data.column(1+startData, countData) += noiseGenerator->noise(data.rows(), countData);
      return Arc(arc.times(), data, arc.getType());
    });

    if(Parallel::isMaster())
    {
      logStatus<<"write instrument data to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
