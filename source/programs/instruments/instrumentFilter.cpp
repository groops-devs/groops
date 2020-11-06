/***********************************************/
/**
* @file instrumentFilter.cpp
*
* @brief Filter any instrument data arc wise.
*
* @author Andreas Kvas
* @date 2016-06-20
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program filter selected data columns of \configFile{inputfileInstrument}{instrument}
with \configClass{digitalFilter}{digitalFilterType} arc wise.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileInstrument.h"
#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Filter any intrument file arc wise.
* @ingroup programsGroup */
class InstrumentFilter
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentFilter, PARALLEL, "filter any instrument file arc wise", Instrument)

/***********************************************/

void InstrumentFilter::run(Config &config)
{
  try
  {
    FileName         fileNameOut, fileNameIn;
    DigitalFilterPtr filter;
    UInt             startData, countData = MAX_UINT;

    readConfig(config, "outputfileInstrument", fileNameOut, Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileInstrument",  fileNameIn,  Config::MUSTSET,  "",  "");
    readConfig(config, "digitalFilter",        filter,      Config::MUSTSET,  "",  "");
    readConfig(config, "startDataFields",      startData,   Config::DEFAULT,  "0", "start");
    readConfig(config, "countDataFields",      countData,   Config::OPTIONAL, "",  "number of data fields (default: all after start)");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument data <"<<fileNameIn<<"> and filter"<<Log::endl;
    InstrumentFile instrumentFile(fileNameIn);

    std::vector<Arc> arcList(instrumentFile.arcCount(), instrumentFile.getType());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      Arc arc = instrumentFile.readArc(arcNo);
      if(arc.size() == 0)
        return arc;
      Matrix data = arc.matrix();
      countData = std::min(countData, data.columns()-1-startData);
      copy(filter->filter(data.column(1+startData, countData)), data.column(1+startData, countData));
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
