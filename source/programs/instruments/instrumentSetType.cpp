/**
* @file instrumentSetType.cpp
*
* @brief Reinterpret data columns of instrument file.
*
* @author Torsten Mayer-Guerr
* @date 2016-11-25
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert \file{instrument data}{instrument} into instrument data with new \configClass{type}{instrumentTypeType}.
The selected number of data columns must agree with the \configClass{type}{instrumentTypeType}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Reinterpret data columns of instrument file.
* @ingroup programsGroup */
class InstrumentSetType
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentSetType, SINGLEPROCESS, "Reinterpret data columns of instrument file", Instrument, Matrix)
GROOPS_RENAMED_PROGRAM(InstrumentMatrix2Instrument, InstrumentSetType, date2time(2020, 02, 15))

/***********************************************/

void InstrumentSetType::run(Config &config)
{
  try
  {
    FileName    fileNameOut, fileNameIn;
    Epoch::Type type = Epoch::MISCVALUESOLD;
    UInt        startData, countData = MAX_UINT;

    renameDeprecatedConfig(config, "inputfileTimeSeries", "inputfileInstrument", date2time(2020, 02, 15));

    readConfig(config, "outputfileInstrument", fileNameOut, Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileInstrument",  fileNameIn,  Config::MUSTSET,  "",  "");
    readConfig(config, "type",                 type,        Config::OPTIONAL, "",  "");
    readConfig(config, "startDataFields",      startData,   Config::DEFAULT,  "0", "start");
    readConfig(config, "countDataFields",      countData,   Config::OPTIONAL, "",  "number of data fields (default: all after start)");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument data <"<<fileNameIn<<">"<<Log::endl;
    InstrumentFile file(fileNameIn);
    countData = std::min(countData, file.dataCount(TRUE/*mustDefined*/)-startData);
    if(type == Epoch::EMPTY)         type = static_cast<Epoch::Type>(countData);
    if(type == Epoch::MISCVALUESOLD) type = Epoch::EMPTY;
    if((Epoch::dataCount(type) != NULLINDEX) && (countData != Epoch::dataCount(type)))
      throw(Exception("Selected type expects "+Epoch::dataCount(type)%"%i data columns, but "s+countData%"%i provided"s));

    std::vector<Arc> arcs;
    for(UInt arcNo=0; arcNo<file.arcCount(); arcNo++)
    {
      Arc arc = file.readArc(arcNo);
      arcs.push_back(Arc(arc.times(), arc.matrix().column(startData, 1+countData), type));
    }

    logStatus<<"write instrument data <"<<fileNameOut<<">"<<Log::endl;
    InstrumentFile::write(fileNameOut, arcs);
    Arc::printStatistics(arcs);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
