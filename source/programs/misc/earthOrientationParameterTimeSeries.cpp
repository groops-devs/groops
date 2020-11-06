/***********************************************/
/**
* @file earthOrientationParameterTimeSeries.cpp
*
* @brief Time series of EOP.
*
* @author Torsten Mayer-Guerr
* @date 2017-05-26
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes a \configClass{timeSeries}{timeSeriesType} (GPS time) of Earth Orientation Parameter (EOP).
The \file{instrument file}{instrument} (MISCVALUES) contains the elements at each epoch in the following order:
\begin{itemize}
\item $x_p$ [rad]
\item $y_p$ [rad]
\item $s_p$ [rad]
\item $UT1-UTC$ [seconds]
\item length of day (LOD) [seconds]
\item $X$ [rad]
\item $Y$ [rad]
\item $S$ [rad]
\end{itemize}
The values are in situ values with all corrections and models applied. The time series can be used to
precompute Earth rotation with a low temporal resolution (e.g. 10 min) and reuse the file in
\configClass{earthRotation:file}{earthRotationType:file} to interpolate the data to the needed epochs
(e.g. to rotate orbit data). As some Earth rotation models are quite slow this can accelerate the computation.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Time series of EOP.
* @ingroup programsGroup */
class EarthOrientationParameterTimeSeries
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(EarthOrientationParameterTimeSeries, PARALLEL, "Time series of EOP", Misc, TimeSeries)

/***********************************************/

void EarthOrientationParameterTimeSeries::run(Config &config)
{
  try
  {
    FileName         fileNameEOP;
    EarthRotationPtr earthRotation;
    TimeSeriesPtr    timeSeries;

    readConfig(config, "outputfileEOP", fileNameEOP,   Config::MUSTSET,  "", "each row: mjd(GPS), xp, yp, sp, dUT1, LOD, X, Y, S");
    readConfig(config, "earthRotation", earthRotation, Config::MUSTSET,  "", "");
    readConfig(config, "timeSeries",    timeSeries,    Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"Computing Earth rotation"<<Log::endl;
    std::vector<Time> times = timeSeries->times();
    Matrix A(times.size(), 9);
    Parallel::forEach(times.size(), [&](UInt i) {earthRotation->earthOrientationParameter(times.at(i), A(i,1), A(i,2), A(i,3), A(i,4), A(i,5), A(i,6), A(i,7), A(i,8));});
    Parallel::reduceSum(A);

    if(Parallel::isMaster())
    {
      logStatus<<"writing EOP to file <"<<fileNameEOP<<">"<<Log::endl;
      InstrumentFile::write(fileNameEOP, Arc(times, A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
