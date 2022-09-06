/***********************************************/
/**
* @file instrumentResample.cpp
*
* @brief Resample data to given time series.
*
* @author Andreas Kvas
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2019-10-07
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program resamples \file{instrument data}{instrument} to a given
\configClass{timeSeries}{timeSeriesType} using a resampling
\configClass{method}{interpolatorTimeSeriesType}.

This program can also be used to reduce the sampling of an instrument file,
but a better way to reduce the sampling of noisy data with regular sampling
is to use a low pass filter first with \program{InstrumentFilter} and then thin
out the data with \program{InstrumentReduceSampling}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/interpolatorTimeSeries/interpolatorTimeSeries.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Resample data to given time series.
* @ingroup programsGroup */
class InstrumentResample
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentResample, SINGLEPROCESS, "Resample data to given time series", Instrument)

/***********************************************/

void InstrumentResample::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                  outName, inName;
    InterpolatorTimeSeriesPtr interpolator;
    TimeSeriesPtr             timeSeries;

    readConfig(config, "outputfileInstrument", outName,      Config::MUSTSET, "", "");
    readConfig(config, "inputfileInstrument",  inName,       Config::MUSTSET, "", "");
    readConfig(config, "method",               interpolator, Config::MUSTSET, "", "resampling method");
    readConfig(config, "timeSeries",           timeSeries,   Config::MUSTSET, "", "resampled points in time");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument data <"<<inName<<">"<<Log::endl;
    const Arc arc = InstrumentFile::read(inName);
    const std::vector<Time> timesNew = timeSeries->times();

    logStatus<<"resample data"<<Log::endl;
    interpolator->init(arc.times(), FALSE/*throwException*/);
    const Matrix A = interpolator->interpolate(timesNew, arc.matrix().column(1, Epoch::dataCount(arc.getType(), TRUE)));

    Arc arcNew;
    Epoch *epoch = Epoch::create(arc.getType());
    for(UInt i=0; i<timesNew.size(); i++)
      if(!std::isnan(A(i,0)))
      {
        epoch->time  = timesNew.at(i);
        epoch->setData(A.row(i).trans());
        arcNew.push_back(*epoch);
    }
    delete epoch;

    logStatus<<"write instrument data to file <"<<outName<<">"<<Log::endl;
    InstrumentFile::write(outName, arcNew);
    Arc::printStatistics(arcNew);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
