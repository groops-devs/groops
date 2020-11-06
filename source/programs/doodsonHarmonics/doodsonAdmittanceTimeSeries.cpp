/***********************************************/
/**
* @file doodsonAdmittanceTimeSeries.cpp
*
* @brief cos/sin multipliers of all major tides.
* Without admittance this would be a simple cos oscillation.
*
* @author Torsten Mayer-Guerr
* @date 2010-05-13
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
To visualize the interpolation of the minor tides it computes cosine multipliers of all major tides.
Without admittance this would be a simple cos oscillation.
The \config{outputfileTimeSeries} is an \file{instrument file}{instrument} containining the cos of all the major tides.

\fig{!hb}{0.8}{doodsonAdmittanceTimeSeries}{fig:doodsonAdmittanceTimeSeries}{Cosine of the Mf tidal frequency with modulation from the interpolated minor tides.}
)";

/***********************************************/

#include "programs/program.h"
#include "base/doodson.h"
#include "files/fileInstrument.h"
#include "files/fileAdmittance.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief cos/sin multipliers of the major tides.
* Without admittance this would be a simple cos oscillation.
* @ingroup programsGroup */
class DoodsonAdmittanceTimeSeries
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(DoodsonAdmittanceTimeSeries, SINGLEPROCESS, "cos/sin multipliers of the major tides.", DoodsonHarmonics, TimeSeries)

/***********************************************/

void DoodsonAdmittanceTimeSeries::run(Config &config)
{
  try
  {
    FileName      fileNameOut, fileNameAdmittance;
    TimeSeriesPtr timeSeries;

    readConfig(config, "outputfileTimeSeries", fileNameOut,        Config::MUSTSET, "", "columns: mjd, cos of major tides");
    readConfig(config, "inputfileAdmittance",  fileNameAdmittance, Config::MUSTSET, "", "cos/sin multipliers of the major tides");
    readConfig(config, "timeSeries",           timeSeries,         Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    // read admittace file
    // -------------------
    logStatus<<"read admittance file <"<<fileNameAdmittance<<">"<<Log::endl;
    Admittance admit;
    readFileAdmittance(fileNameAdmittance, admit);

    // Computation
    // -----------
    logStatus<<"computing time series"<<Log::endl;
    std::vector<Time> times  = timeSeries->times();
    Matrix doodsonMatrix = Doodson::matrix(admit.doodsonMinor); // Matrix with Doodson multiplicators
    Matrix A(times.size(), 1+admit.doodsonMajor.size());

    logTimerStart;
    for(UInt i=0; i<times.size(); i++)
    {
      logTimerLoop(i,times.size());

      Vector thetaf = doodsonMatrix * Doodson::arguments(times.at(i));
      Vector cosMinor(thetaf.rows());
      for(UInt k=0; k<thetaf.rows(); k++)
        cosMinor(k) = std::cos(thetaf(k));
      Vector cosMajor = admit.admittance * cosMinor;

      A(i,0) = times.at(i).mjd();
      for(UInt k=0; k<cosMajor.rows(); k++)
        A(i,1+k) = cosMajor(k);
    }
    logTimerLoopEnd(times.size());

    // save results
    // ------------
    logStatus<<"write time series to file <"<<fileNameOut<<">"<<Log::endl;
    for(UInt i=0; i<admit.doodsonMajor.size(); i++)
      logInfo<<"  "<<admit.doodsonMajor.at(i).code()<<" "<<admit.doodsonMajor.at(i).name()<<Log::endl;
    writeFileInstrument(fileNameOut, Arc(times, A));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
