/***********************************************/
/**
* @file timeSeriesCreate.cpp
*
* @brief Create a time series.
*
* @author Norbert Zehentner
* @date 2016-04-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program generates an \file{instrument file}{instrument},
containing a time series.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Create a time series.
* @ingroup programsGroup */
class TimeSeriesCreate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(TimeSeriesCreate, SINGLEPROCESS, "Create a time series", Misc, TimeSeries)
GROOPS_RENAMED_PROGRAM(InstrumentCreateTimeSeries, TimeSeriesCreate, date2time(2020, 01, 29))

/***********************************************/

void TimeSeriesCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName      fileNameTimeSeries;
    TimeSeriesPtr timeSeries;
    std::vector<ExpressionVariablePtr> dataExpr;

    renameDeprecatedConfig(config, "outputfileFile", "outputfileTimeSeries", date2time(2020, 01, 29));

    readConfig(config, "outputfileTimeSeries", fileNameTimeSeries, Config::MUSTSET,  "", "instrument file");
    readConfig(config, "timeSeries",           timeSeries,         Config::MUSTSET,  "", "time series to be created");
    readConfig(config, "data",                 dataExpr,           Config::OPTIONAL, "", "expression of output columns, extra 'epoch' variable");
    if(isCreateSchema(config)) return;

    logStatus<<"create time series"<<Log::endl;
    std::vector<Time> times = timeSeries->times();

    auto varList = config.getVarList();
    std::set<std::string> usedVariables;
    std::for_each(dataExpr.begin(), dataExpr.end(), [&](auto expr) {expr->usedVariables(varList, usedVariables);});
    addDataVariables("epoch", times, varList, usedVariables);
    std::for_each(dataExpr.begin(), dataExpr.end(), [&](auto expr) {expr->simplify(varList);});

    Matrix A(times.size(), 1+dataExpr.size());
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      varList["epoch"]->setValue(times.at(idEpoch).mjd());
      for(UInt k=0; k<dataExpr.size(); k++)
        A(idEpoch, 1+k) = dataExpr.at(k)->evaluate(varList);
    }

    logStatus<<"write instrument data to file <"<<fileNameTimeSeries<<">"<<Log::endl;
    Arc arc(times, A);
    InstrumentFile::write(fileNameTimeSeries, arc);
    Arc::printStatistics(arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
