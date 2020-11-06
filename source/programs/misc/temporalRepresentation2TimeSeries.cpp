/***********************************************/
/**
* @file temporalRepresentation2TimeSeries.cpp
*
* @brief Design matrix of temporal representation.
*
* @author Torsten Mayer-Guerr
* @date 2015-07-07
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the design matrix of temporal representation at a given time series.
The output matrix contains the time steps in MJD in the first column, the other columns contain the design matrix.
The intention of this program is to visualize the parametrization together with \program{PlotGraph}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Design matrix of temporal representation.
* @ingroup programsGroup */
class TemporalRepresentation2TimeSeries
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(TemporalRepresentation2TimeSeries, SINGLEPROCESS, "Design matrix of temporal representation.", Misc, TimeSeries)

/***********************************************/

void TemporalRepresentation2TimeSeries::run(Config &config)
{
  try
  {
    FileName                   fileNameMatrix;
    TimeSeriesPtr              timeSeries;
    ParametrizationTemporalPtr temporal;

    readConfig(config, "outputfileMatrix", fileNameMatrix, Config::MUSTSET, "", "Time (MJD) in first column, design matrix follows");
    readConfig(config, "timeSeries",       timeSeries,     Config::MUSTSET,  "", "");
    readConfig(config, "temporal",         temporal,       Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"Compute design matrix"<<Log::endl;
    std::vector<Time> times = timeSeries->times();
    Matrix A(times.size(), 1+temporal->parameterCount());
    for(UInt i=0; i<times.size(); i++)
    {
      A(i,0) = times.at(i).mjd();
      copy(temporal->factors(times.at(i)).trans(), A.slice(i,1,1,temporal->parameterCount()));
    }

    logStatus<<"write matrix to file <"<<fileNameMatrix<<">"<<Log::endl;
    writeFileMatrix(fileNameMatrix, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
