/***********************************************/
/**
* @file instrumentStarCameraMultiply.cpp
*
* @brief Mutiply instrument data with a factor and add them together.
*
* @author Torsten Mayer-Guerr
* @date 2012-06-24
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program applies several rotations given by \configFile{inputfileStarCamera}{instrument}.
The resulting rotation is written as \configFile{outputfileStarCamera}{instrument}.
All instrument files must be synchronized (\program{InstrumentSynchronize}).
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Mutiply instrument data with a factor and add them together.
* @ingroup programsGroup */
class InstrumentStarCameraMultiply
{
public:
  class Data
  {
    public:
    FileName fileName;
    Bool     inverse;
  };

  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentStarCameraMultiply, SINGLEPROCESS, "Mutiply instrument data with a factor and add them together", Instrument)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, InstrumentStarCameraMultiply::Data &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileStarCamera", var.fileName, Config::MUSTSET,  "", "");
  readConfig(config, "inverse",             var.inverse,  Config::DEFAULT,  "0", "");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void InstrumentStarCameraMultiply::run(Config &config)
{
  try
  {
    FileName outName;
    std::vector<Data> data;

    readConfig(config, "outputfileStarCamera", outName, Config::MUSTSET,  "", "");
    readConfig(config, "instrument",           data,    Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    std::vector<InstrumentFilePtr> instrumentFile;
    instrumentFile.resize(data.size());
    for(UInt i=0; i<instrumentFile.size(); i++)
      instrumentFile.at(i) = InstrumentFile::newFile(data.at(i).fileName);
    for(UInt i=1; i<instrumentFile.size(); i++)
      InstrumentFile::checkArcCount({*instrumentFile.at(0), *instrumentFile.at(i)});

    // read data
    // ---------
    logStatus<<"read starCamera data"<<Log::endl;
    UInt arcCount = instrumentFile.at(0)->arcCount();
    std::vector<Arc> arcList;

    logTimerStart;
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
    {
      logTimerLoop(arcNo, arcCount);

      StarCameraArc starCamera = instrumentFile.at(0)->readArc(arcNo);
      UInt epochCount = starCamera.size();
      if(data.at(0).inverse)
        for(UInt k=0; k<epochCount; k++)
          starCamera.at(k).rotary = inverse(starCamera.at(k).rotary);

      for(UInt i=1; i<instrumentFile.size(); i++)
      {
        StarCameraArc starCamera2 = instrumentFile.at(i)->readArc(arcNo);
        Arc::checkSynchronized({starCamera, starCamera2});
        for(UInt k=0; k<epochCount; k++)
        {
          if(data.at(i).inverse)
            starCamera.at(k).rotary = starCamera.at(k).rotary * inverse(starCamera2.at(k).rotary);
          else
            starCamera.at(k).rotary = starCamera.at(k).rotary * starCamera2.at(k).rotary;
        }
      }
      arcList.push_back(starCamera);
    }
    logTimerLoopEnd(arcCount);

    // save file
    // ---------
    logStatus<<"write instrument data to file <"<<outName<<">"<<Log::endl;
    InstrumentFile::write(outName, arcList);

    Arc::printStatistics(arcList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
