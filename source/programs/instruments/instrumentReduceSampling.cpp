/***********************************************/
/**
* @file instrumentReduceSampling.cpp
*
* @brief Reduce sampling of instrument data.
*
* @author Torsten Mayer-Guerr
* @date 2010-01-06
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reduce the sampling of a instrument file. Only epochs with a time stamp
with a division by \config{sampling} without remainder are kept (inside \config{margin}).
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Reduce sampling of instrument data.
* @ingroup programsGroup */
class InstrumentReduceSampling
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentReduceSampling, SINGLEPROCESS, "reduce sampling of instrument data.", Instrument)

/***********************************************/

void InstrumentReduceSampling::run(Config &config)
{
  try
  {
    FileName inName, outName;
    Double   seconds, margin;
    Bool     relative2FirstEpoch;

    readConfig(config, "outputfileInstrument", outName, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrument",  inName,  Config::MUSTSET,  "", "");
    readConfig(config, "sampling",             seconds, Config::MUSTSET,   "",     "new sampling in seconds");
    readConfig(config, "margin",               margin,  Config::DEFAULT,   "1e-5", "margin around the new sampling in seconds");
    readConfig(config, "relative2FirstEpoch",  relative2FirstEpoch,  Config::DEFAULT,   "0", "compute sampling relative to time of first epoch");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument data <"<<inName<<"> and reduce sampling"<<Log::endl;
    InstrumentFile  inFile(inName);
    UInt arcCount = inFile.arcCount();
    std::list<Arc> arcList;

    logTimerStart;
    Double refTime=0;
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
    {
      logTimerLoop(arcNo, arcCount);

      std::vector<Time> times;
      Arc arc = inFile.readArc(arcNo);
      if(arcNo==0 && relative2FirstEpoch)
        refTime = arc.at(0).time.mjdMod();
      for(UInt i=0; i<arc.size(); i++)
      {

        Double mod = std::fmod((arc.at(i).time.mjdMod()-refTime)*86400, seconds);

        if(mod/seconds > 0.5)
          mod -= seconds;
        if(fabs(mod) <= margin)
          if((times.size()==0) || (arc.at(i).time-times.back()).seconds()>(seconds-margin))
            times.push_back(arc.at(i).time);
      }

      arc.synchronize(times, margin);
      if(arc.size()!=0)
        arcList.push_back(arc);
    } // for(arc)
    logTimerLoopEnd(arcCount);

    logStatus<<"write instrument data to <"<<outName<<">"<<Log::endl;
    InstrumentFile::write(outName, arcList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
