/***********************************************/
/**
* @file noiseTimeSeries.cpp
*
* @brief Generate one or more time series of noise.
*
* @author Matthias Ellmer
* @date 2015-11-03
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program generates \configFile{outputfileNoise}{instrument} with the requested characteristics.
See \configClass{noiseGenerator}{noiseGeneratorType} for details on noise options.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Generate one or more time series of noise.
* @ingroup programsGroup */
class NoiseTimeSeries
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(NoiseTimeSeries, SINGLEPROCESS, "generate random noise", Simulation, Noise, TimeSeries)

/***********************************************/

void NoiseTimeSeries::run(Config& config)
{
  try
  {
    FileName          outNameNoise, outNameCovariance;
    NoiseGeneratorPtr noiseGenerator;
    TimeSeriesPtr     timeSeries;
    UInt              series;

    readConfig(config, "outputfileNoise",              outNameNoise,      Config::MUSTSET,  "",  "");
    readConfig(config, "outputfileCovarianceFunction", outNameCovariance, Config::OPTIONAL, "",  "");
    readConfig(config, "noise",                        noiseGenerator,    Config::MUSTSET,  "",  "");
    readConfig(config, "timeSeries",                   timeSeries,        Config::MUSTSET,  "",  "");
    readConfig(config, "columns",                      series,            Config::MUSTSET,  "1", "number of noise series (columns)");
    if(isCreateSchema(config)) return;

    // Set up time series
    std::vector<Time> times = timeSeries->times();

    logStatus<<"generating noise"<<Log::endl;
    Matrix noise(times.size(), 1+series);
    copy(noiseGenerator->noise(noise.rows(), series), noise.column(1, series));

    logStatus<<"writing noise to file <"<< outNameNoise <<">"<<Log::endl;
    InstrumentFile::write(outNameNoise, Arc(times, noise));

    if(!outNameCovariance.empty())
    {
      logStatus<<"writing covariance function to file <"<< outNameCovariance <<">"<<Log::endl;
      writeFileMatrix(outNameCovariance, noiseGenerator->covarianceFunction(times.size(), medianSampling(times).seconds()));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
