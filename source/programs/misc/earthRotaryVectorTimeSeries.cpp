/***********************************************/
/**
* @file earthRotaryVectorTimeSeries.cpp
*
* @brief Time series of Earth's rotary axis.
*
* @author Torsten Mayer-Guerr
* @date 2017-06-08
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes a \configFile{outputfileTimeSeries}{instrument} of Earth's rotary axis
and its temporal derivative at \configClass{timeSeries}{timeSeriesType} (GPS time).
The \file{instrument file}{instrument} (MISCVALUES) contains the elements at each epoch in the following order:
\begin{itemize}
\item $\omega_x [rad/s]$
\item $\omega_y [rad/s]$
\item $\omega_z [rad/s]$
\item $\dot{\omega}_x [rad/s^2]$
\item $\dot{\omega}_y [rad/s^2]$
\item $\dot{\omega}_z [rad/s^2]$.
\end{itemize}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Time series of Earth's rotary axis.
* @ingroup programsGroup */
class EarthRotaryVectorTimeSeries
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(EarthRotaryVectorTimeSeries, PARALLEL, "Time series of Earth's rotary axis", Misc, TimeSeries)

/***********************************************/

void EarthRotaryVectorTimeSeries::run(Config &config)
{
  try
  {
    FileName         fileNameOut;
    EarthRotationPtr earthRotation;
    TimeSeriesPtr    timeSeries;
    Bool             inTRF = FALSE;

    readConfig(config, "outputfileTimeSeries", fileNameOut,   Config::MUSTSET,  "", "wx, wy, wz [rad], dwx, dwy, dwz [rad/s^2]");
    readConfig(config, "earthRotation",        earthRotation, Config::MUSTSET,  "", "");
    readConfig(config, "timeSeries",           timeSeries,    Config::MUSTSET,  "", "");
    readConfig(config, "inTRF",                inTRF,         Config::DEFAULT,  "0", "terrestrial reference frame, otherwise celestial");
    if(isCreateSchema(config)) return;


    logStatus<<"Computing Earth rotation"<<Log::endl;
    std::vector<Time> times = timeSeries->times();
    Matrix A(times.size(), 7);
    Parallel::forEach(times.size(), [&](UInt i)
    {
      Rotary3d rotEarth;
      if(inTRF)
        rotEarth = earthRotation->rotaryMatrix(times.at(i));
      const Vector3d axis    = rotEarth.rotate(earthRotation->rotaryAxis(times.at(i)));
      const Vector3d axisDot = rotEarth.rotate(earthRotation->rotaryAxisDerivate(times.at(i)));
      A(i,1) = axis.x();
      A(i,2) = axis.y();
      A(i,3) = axis.z();
      A(i,4) = axisDot.x();
      A(i,5) = axisDot.y();
      A(i,6) = axisDot.z();
    });
    Parallel::reduceSum(A);

    if(Parallel::isMaster())
    {
      logStatus<<"writing time series to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, Arc(times, A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
