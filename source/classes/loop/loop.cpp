/***********************************************/
/**
* @file loop.cpp
*
* @brief Loop.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#define DOCSTRING_Loop

#include "base/import.h"
#include "config/configRegister.h"
#include "loopTimeSeries.h"
#include "loopTimeIntervals.h"
#include "loopManualList.h"
#include "loopManualTable.h"
#include "loopFileAscii.h"
#include "loopFileAsciiTable.h"
#include "loopMatrix.h"
#include "loopUniformSampling.h"
#include "loopCommandOutput.h"
#include "loopLoop.h"
#include "loopPlatformEquipment.h"
#include "loopFileGnssStationInfo.h"
#include "loop.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Loop, "loopType",
                      LoopTimeSeries,
                      LoopTimeIntervals,
                      LoopManualList,
                      LoopManualTable,
                      LoopFileAscii,
                      LoopFileAsciiTable,
                      LoopMatrix,
                      LoopUniformSampling,
                      LoopCommandOutput,
                      LoopLoop,
                      LoopPlatformEquipment,
                      LoopFileGnssStationInfo)


GROOPS_READCONFIG_CLASS(Loop, "loopType")

/***********************************************/

LoopPtr Loop::create(Config &config, const std::string &name)
{
  try
  {
    LoopPtr     loop;
    std::string type;
    readConfigChoice(config, name, type, Config::MUSTSET, "", "");

    renameDeprecatedChoice(config, type, "temporal", "timeIntervals", date2time(2018, 5, 31));

    if(readConfigChoiceElement(config, "timeSeries",          type, "Loop over points in time"))
      loop = LoopPtr(new LoopTimeSeries(config));
    if(readConfigChoiceElement(config, "timeIntervals",       type, "Loop over time intervals"))
      loop = LoopPtr(new LoopTimeIntervals(config));
    if(readConfigChoiceElement(config, "manualList",          type, "Loop over list of strings"))
      loop = LoopPtr(new LoopManualList(config));
    if(readConfigChoiceElement(config, "manualTable",         type, "Loop over table of strings"))
      loop = LoopPtr(new LoopManualTable(config));
    if(readConfigChoiceElement(config, "fileAscii",           type, "Loop over list of strings from file"))
      loop = LoopPtr(new LoopFileAscii(config));
    if(readConfigChoiceElement(config, "fileAsciiTable",      type, "Loop over of a table containing strings"))
      loop = LoopPtr(new LoopFileAsciiTable(config));
    if(readConfigChoiceElement(config, "matrix",              type, "Loop over rows of a matrix"))
      loop = LoopPtr(new LoopMatrix(config));
    if(readConfigChoiceElement(config, "uniformSampling",     type, "Loop over sequence of numbers"))
      loop = LoopPtr(new LoopUniformSampling(config));
    if(readConfigChoiceElement(config, "commandOutput",       type, "Loop over lines of command output"))
      loop = LoopPtr(new LoopCommandOutput(config));
    if(readConfigChoiceElement(config, "loop",                type, "Loop over loops"))
      loop = LoopPtr(new LoopLoop(config));
    if(readConfigChoiceElement(config, "platformEquipment",   type, "Loop over equipment of a platform file"))
      loop = LoopPtr(new LoopPlatformEquipment(config));
    if(readConfigChoiceElement(config, "fileGnssStationInfo", type, "DEPRECATED. Use platformEquipment instead"))
      loop = LoopPtr(new LoopFileGnssStationInfo(config));
    endChoice(config);

    return loop;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
