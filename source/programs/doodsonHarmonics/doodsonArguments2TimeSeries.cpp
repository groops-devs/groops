/***********************************************/
/**
* @file doodsonArguments2TimeSeries.cpp
*
* @brief time series of doodson/fundamental arguments.
*
* @author Torsten Mayer-Guerr
* @date 2014-10-10
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Time series of doodson/fundamental arguments.
The \configFile{outputfileTimeSeries}{instrument} contains the six Doodson arguments,
followed by the five fundamental arguments in radians.
)";

/***********************************************/

#include "programs/program.h"
#include "base/planets.h"
#include "base/doodson.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Time series of doodson/fundamental arguments.
* @ingroup programsGroup */
class DoodsonArguments2TimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(DoodsonArguments2TimeSeries, SINGLEPROCESS, "doodson/fundamental arguments.", DoodsonHarmonics, TimeSeries)

/***********************************************/

void DoodsonArguments2TimeSeries::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName      fileNameOut;
    TimeSeriesPtr timeSeries;

    readConfig(config, "outputfileTimeSeries", fileNameOut, Config::MUSTSET, "", "each epoch: 6 doodson args, 5 fundamental args [rad]");
    readConfig(config, "timeSeries",           timeSeries,  Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    // Computation
    // -----------
    logStatus<<"computing time series"<<Log::endl;
    std::vector<Time> times  = timeSeries->times();
    Matrix A(times.size(), 1+6+5);

    Single::forEach(times.size(), [&](UInt i)
    {
      A(i,0) = times.at(i).mjd();
      copy(Doodson::arguments(times.at(i)).trans(),    A.slice(i,1,1,6));
      copy(Planets::fundamentals(times.at(i)).trans(), A.slice(i,7,1,5));
    });

    // save results
    // ------------
    logStatus<<"write time series to file <"<<fileNameOut<<">"<<Log::endl;
    writeFileInstrument(fileNameOut, Arc(times, A));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
