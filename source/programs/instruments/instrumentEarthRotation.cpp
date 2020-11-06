/***********************************************/
/**
* @file instrumentEarthRotation.cpp
*
* @brief Precompute Earth rotation matrix (CRF -> TRF) and save as StarCamera file.
*
* @author Torsten Mayer-Guerr
* @date 2009-05-12
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Precompute Earth rotation matrix from celestial to terrestrial frame
and save as \file{StarCamera file}{instrument}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/timeSeries/timeSeries.h"

/***********************************************/

/** @brief Precompute Earth rotation matrix (CRF -> TRF) and save as StarCamera file.
* @ingroup programsGroup */
class InstrumentEarthRotation
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentEarthRotation, SINGLEPROCESS, "Precompute Earth rotation matrix (CRF -> TRF) and save as StarCamera file.", Instrument)

/***********************************************/

void InstrumentEarthRotation::run(Config &config)
{
  try
  {
    FileName         fileNameStarCamera;
    EarthRotationPtr earthRotation;
    TimeSeriesPtr    timeSeries;

    readConfig(config, "outputfileStarCamera", fileNameStarCamera, Config::MUSTSET, "", "rotation from CRF to TRF");
    readConfig(config, "earthRotation",        earthRotation,      Config::MUSTSET, "", "");
    readConfig(config, "timeSeries",           timeSeries,         Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    std::vector<Time> times = timeSeries->times();

    logStatus<<"computing earth rotation"<<Log::endl;
    Arc arc;
    logTimerStart;
    for(UInt i=0; i<times.size(); i++)
    {
      logTimerLoop(i, times.size());
      StarCameraEpoch epoch;
      epoch.time = times.at(i);
      epoch.rotary = earthRotation->rotaryMatrix(times.at(i));
      arc.push_back(epoch);
    }
    logTimerLoopEnd(times.size());

    logStatus<<"write star camera data to file <"<<fileNameStarCamera<<">"<<Log::endl;
    InstrumentFile::write(fileNameStarCamera, arc);
    Arc::printStatistics(arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
