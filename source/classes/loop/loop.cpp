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
#include "loopStringList.h"
#include "loopStringTable.h"
#include "loopFileStringList.h"
#include "loopFileStringTable.h"
#include "loopFileTextLines.h"
#include "loopMatrix.h"
#include "loopNumberSequence.h"
#include "loopDirectoryListing.h"
#include "loopCommandOutput.h"
#include "loopPlatformEquipment.h"
#include "loopLoop.h"
#include "loopSortAndRemoveDuplicates.h"
#include "loopManualTable.h"
#include "loopFileGnssStationInfo.h"
#include "loop.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Loop, "loopType",
                      LoopTimeSeries,
                      LoopTimeIntervals,
                      LoopStringList,
                      LoopStringTable,
                      LoopFileStringList,
                      LoopFileStringTable,
                      LoopFileTextLines,
                      LoopMatrix,
                      LoopNumberSequence,
                      LoopDirectoryListing,
                      LoopCommandOutput,
                      LoopPlatformEquipment,
                      LoopLoop,
                      LoopSortAndRemoveDuplicates,
                      LoopManualTable,
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

    renameDeprecatedChoice(config, type, "temporal",        "timeIntervals",   date2time(2018, 5, 31));
    renameDeprecatedChoice(config, type, "manualList",      "stringList",      date2time(2025, 9, 26));
    // renameDeprecatedChoice(config, type, "manualTable",     "stringTable",     date2time(2025, 9, 26));
    renameDeprecatedChoice(config, type, "fileAscii",       "fileStringList",  date2time(2025, 9, 26));
    renameDeprecatedChoice(config, type, "fileAsciiTable",  "fileStringTable", date2time(2025, 9, 26));
    renameDeprecatedChoice(config, type, "fileLines",       "fileTextLines",   date2time(2025, 9, 26));
    renameDeprecatedChoice(config, type, "uniformSampling", "numberSequence",  date2time(2025, 9, 26));

    if(readConfigChoiceElement(config, "timeSeries",              type, "Loop over points in time"))
      loop = LoopPtr(new LoopTimeSeries(config));
    if(readConfigChoiceElement(config, "timeIntervals",           type, "Loop over time intervals"))
      loop = LoopPtr(new LoopTimeIntervals(config));
    if(readConfigChoiceElement(config, "stringList",              type, "Loop over list of strings"))
      loop = LoopPtr(new LoopStringList(config));
    if(readConfigChoiceElement(config, "stringTable",             type, "Loop over table of strings"))
      loop = LoopPtr(new LoopStringTable(config));
    if(readConfigChoiceElement(config, "fileStringList",          type, "Loop over list of strings from file"))
      loop = LoopPtr(new LoopFileStringList(config));
    if(readConfigChoiceElement(config, "fileStringTable",         type, "Loop over of a table containing strings"))
      loop = LoopPtr(new LoopFileStringTable(config));
    if(readConfigChoiceElement(config, "fileTextLines",           type, "Loop over of lines of a text file"))
      loop = LoopPtr(new LoopFileTextLines(config));
    if(readConfigChoiceElement(config, "matrix",                  type, "Loop over rows of a matrix"))
      loop = LoopPtr(new LoopMatrix(config));
    if(readConfigChoiceElement(config, "numberSequence",          type, "Loop over sequence of numbers"))
      loop = LoopPtr(new LoopNumberSequence(config));
    if(readConfigChoiceElement(config, "directoryListing",        type, "Loop over files of a directory"))
      loop = LoopPtr(new LoopDirectoryListing(config));
    if(readConfigChoiceElement(config, "commandOutput",           type, "Loop over lines of command output"))
      loop = LoopPtr(new LoopCommandOutput(config));
    if(readConfigChoiceElement(config, "platformEquipment",       type, "Loop over equipment of a platform file"))
      loop = LoopPtr(new LoopPlatformEquipment(config));
    if(readConfigChoiceElement(config, "loop",                    type, "Loop over loops"))
      loop = LoopPtr(new LoopLoop(config));
    if(readConfigChoiceElement(config, "sortAndRemoveDuplicates", type, "sort loops"))
      loop = LoopPtr(new LoopSortAndRemoveDuplicates(config));
    if(readConfigChoiceElement(config, "manualTable",             type, "DEPRECATED since 2025-09-27. Use stringTable instead"))
      loop = LoopPtr(new LoopManualTable(config));
    if(readConfigChoiceElement(config, "fileGnssStationInfo",     type, "DEPRECATED since 2022-11-11. Use platformEquipment instead"))
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

void Loop::readConfigCondition(Config &config)
{
  ConditionPtr condition;
  isCondition = readConfigLater(config, "condition", condition, conditionConfig, Config::OPTIONAL, "", "check before each loop step");
  index_ = 0;
}

/***********************************************/

Bool Loop::checkCondition(VariableList &varList)
{
  try
  {
    index_++;
    if(!isCondition)
      return TRUE;
    ConditionPtr condition;
    conditionConfig.read(condition, varList);
    // is condition false? -> perfrom directly next iteration
    return (!condition || condition->condition(varList)) ? TRUE : iteration(varList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
