/***********************************************/
/**
* @file instrument2RmsPlotGrid.cpp
*
* @brief Compute RMS plot grid from instrument file(s) containing 3D data (e.g. orbits, station positions).
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2016-11-20
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes an RMS plot grid from one or more \configFile{inputfileInstrument}{instrument}
containing 3D data (e.g. orbits or station positions), which can then be plotted as gridded data in \program{PlotGraph}.
The RMS is computed from the difference between \configFile{inputfileInstrument}{instrument} and
\configFile{inputfileInstrumentReference}{instrument}.
All instrument files must be synchronized (see \program{InstrumentSynchronize}).

Each separate \configFile{inputfileInstrument}{instrument} represents an entry (e.g. a satellite or station)
in the resulting grid. Therefore, providing, for example, 32 orbit files of GPS satellites
results in a grid with columns: mjd, id (0-31), rms.

The first three data columns of the instrument data are considered for computation of the RMS values.
The \config{factor} can be set to, for example, sqrt(3) to get 3D instead of 1D RMS values.

If \configClass{timeIntervals}{timeSeriesType} are provided, each \configFile{inputfileInstrument}{instrument}
and \configFile{inputfileInstrumentReference}{instrument} serves as a template with variable \verb|loopTime|.
This allows concatenation of instrument files, for example to create a month-long RMS plot grid from daily GPS
orbit files (see below).

Helmert parameters between the two frames can be estimated each epoch optionally if
\config{estimateShift}, \config{estimateScale}, or \config{estimateRotation} are set.
It uses a \reference{robust least squares adjustment}{fundamentals.robustLeastSquares}.

\fig{!hb}{0.8}{instrument2RmsPlotGrid}{fig:instrument2RmsPlotGrid}{Comparison of estimated GPS orbits with IGS final solution.}
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/varianceComponentEstimation.h"

/***** CLASS ***********************************/

/** @brief Compute RMS plot grid from instrument file(s) containing 3D data (e.g. orbits, station positions).
* @ingroup programsGroup */
class Instrument2RmsPlotGrid
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Instrument2RmsPlotGrid, SINGLEPROCESS, "Compute RMS plot grid from instrument file(s) containing 3D data (e.g. orbits, station positions).", Instrument, Plot, Gnss)

/***********************************************/

void Instrument2RmsPlotGrid::run(Config &config)
{
  try
  {
    FileName              fileNameOutRmsPlotGrid;
    FileName              fileNameHelmertTimeSeries;
    std::vector<FileName> fileNameInInstrument, fileNameInInstrumentRef;
    TimeSeriesPtr         timesIntervalPtr;
    Double                factor;
    Bool                  estimateShift, estimateScale, estimateRotation;
    Double                huber, huberPower;
    UInt                  iterCount;

    renameDeprecatedConfig(config, "intervals", "timeIntervals", date2time(2020, 7, 14));

    readConfig(config, "outputfileRmsPlotGrid",        fileNameOutRmsPlotGrid,    Config::OPTIONAL, "",     "columns: mjd, id, rms");
    readConfig(config, "outputfileHelmertTimeSeries",  fileNameHelmertTimeSeries, Config::OPTIONAL, "",     "columns: mjd, tx, ty, tz, scale, rx, ry, rz");
    readConfig(config, "inputfileInstrument",          fileNameInInstrument,      Config::MUSTSET,  "",     "one file per satellite/station");
    readConfig(config, "inputfileInstrumentReference", fileNameInInstrumentRef,   Config::MUSTSET,  "",     "one file per satellite/station, same order as above");
    readConfig(config, "timeIntervals",                timesIntervalPtr,          Config::DEFAULT,  "",     "for {loopTime} variable in inputfile");
    readConfig(config, "factor",                       factor,                    Config::DEFAULT,  "1",    "e.g. sqrt(3) for 3D RMS");
    readConfig(config, "estimateShift",                estimateShift,             Config::DEFAULT,  "1",    "coordinate center every epoch");
    readConfig(config, "estimateScale",                estimateScale,             Config::DEFAULT,  "1",    "scale factor of position every epoch");
    readConfig(config, "estimateRotation",             estimateRotation,          Config::DEFAULT,  "1",    "rotation every epoch");
    readConfig(config, "huber",                        huber,                     Config::DEFAULT,  "2.5", "for robust least squares");
    readConfig(config, "huberPower",                   huberPower,                Config::DEFAULT,  "1.5", "for robust least squares");
    readConfig(config, "huberMaxIteration",            iterCount,                 Config::DEFAULT,  "30",  "(maximum) number of iterations for robust estimation");
    if(isCreateSchema(config)) return;

    // ======================================================

    // init time intervals
    // -------------------
    std::vector<Time> timesInterval;
    if(timesIntervalPtr)
      timesInterval = timesIntervalPtr->times();

    VariableList fileNameVariableList;
    if(timesInterval.size())
      addTimeVariables(fileNameVariableList);

    const UInt countInterval   = timesInterval.size() ? timesInterval.size()-1 : 1;
    const UInt countInstrument = fileNameInInstrument.size();

    // ======================================================

    // read instrument data and compute RMS plot grid
    // -----------------------------------------
    logStatus<<"Reading instrument data and computing RMS plot grid"<<Log::endl;
    if(fileNameInInstrument.size() != fileNameInInstrumentRef.size())
      throw Exception("Number of instrument and instrument reference files does not match: " + fileNameInInstrument.size() % "%i"s + " != " + fileNameInInstrumentRef.size() % "%i"s);

    Matrix RmsPlotGrid(countInterval*countInstrument,3);
    std::vector<Vector> helmertTimeSeries;
    UInt countOutput = 0;

    logTimerStart;
    for(UInt idInterval = 0; idInterval < countInterval; idInterval++)
    {
      logTimerLoop(idInterval, countInterval);

      Time time;
      if(timesInterval.size())
      {
        evaluateTimeVariables(idInterval, timesInterval.at(idInterval), timesInterval.at(idInterval+1), fileNameVariableList);
        time = 0.5*(timesInterval.at(idInterval+1)+timesInterval.at(idInterval));
      }

      // read instrument data of current interval
      // ----------------------------------------
      UInt epochCount = MAX_UINT;
      std::vector<Arc> arcsInstrument(countInstrument), arcsInstrumentRef(countInstrument);
      for(UInt idInstrument=0; idInstrument<countInstrument; idInstrument++)
      {
        Arc arc;
        try
        {
          arc = InstrumentFile::read(fileNameInInstrument.at(idInstrument)(fileNameVariableList));
        }
        catch(std::exception &) {} // --> arc.size() == 0, see subsequent if()

        if(!arc.size())
        {
          //logWarning << "Instrument file <" << fileNameInInstrument.at(idInstrument)(fileNameVariableList) << "> not found or empty, skipping." << Log::endl;
          continue;
        }

        Arc arcRef;
        try
        {
          arcRef = InstrumentFile::read(fileNameInInstrumentRef.at(idInstrument)(fileNameVariableList));
        }
        catch(std::exception &) {} // --> arcRef.size() == 0, see subsequent if()

        if(!arcRef.size())
        {
          //logWarning << "Instrument reference file <" << fileNameInInstrumentRef.at(idInstrument)(fileNameVariableList) << "> not found or empty, skipping." << Log::endl;
          continue;
        }

        // check if data of reference file matches instrument file
        if(arc.size() != arcRef.size())
          throw Exception("Instrument and instrument reference epoch counts do not match: " + arc.size() % "%i"s + " != " + arcRef.size() % "%i"s + " for files <" +
                          fileNameInInstrument.at(idInstrument)(fileNameVariableList).str() + "> and <" + fileNameInInstrumentRef.at(idInstrument)(fileNameVariableList).str() + ">");

        arcsInstrument.at(idInstrument)    = arc;
        arcsInstrumentRef.at(idInstrument) = arcRef;
        epochCount = std::min(epochCount, arc.size());
      }

      if(epochCount == MAX_UINT)
      {
        logWarning<<"No instrument data found. continue with next epoch"<<Log::endl;
        continue;
      }

      // Estimate Helmert transformation
      // -------------------------------
      Matrix l(3*arcsInstrument.size(), epochCount);
      for(UInt idEpoch=0; idEpoch<epochCount; idEpoch++)
      {
        Time timeEpoch;
        Matrix A(3*arcsInstrument.size(), 3*estimateShift+estimateScale+3*estimateRotation);
        for(UInt idInstrument=0; idInstrument<countInstrument; idInstrument++)
        {
          if(arcsInstrument.at(idInstrument).size() == 0)
            continue;
          timeEpoch = arcsInstrument.at(idInstrument).at(idEpoch).time;

          const Vector3d pos(arcsInstrument.at(idInstrument).at(idEpoch).data().row(0,3));
          const Vector3d posRef(arcsInstrumentRef.at(idInstrument).at(idEpoch).data().row(0,3));
          const Vector3d diff = pos - posRef;
          l(3*idInstrument+0, idEpoch) = diff.x();
          l(3*idInstrument+1, idEpoch) = diff.y();
          l(3*idInstrument+2, idEpoch) = diff.z();

          UInt idx = 0;
          if(estimateShift)
          {
            A(3*idInstrument+0, idx++) = 1;
            A(3*idInstrument+1, idx++) = 1;
            A(3*idInstrument+2, idx++) = 1;
          }
          if(estimateScale)
          {
            A(3*idInstrument+0, idx)   = posRef.x()/DEFAULT_R;
            A(3*idInstrument+1, idx)   = posRef.y()/DEFAULT_R;
            A(3*idInstrument+2, idx++) = posRef.z()/DEFAULT_R;
          }
          if(estimateRotation)
          {
            A(3*idInstrument+1, idx)   =  posRef.z()/DEFAULT_R; //rx
            A(3*idInstrument+2, idx++) = -posRef.y()/DEFAULT_R; //rx
            A(3*idInstrument+0, idx)   = -posRef.z()/DEFAULT_R; //ry
            A(3*idInstrument+2, idx++) =  posRef.x()/DEFAULT_R; //ry
            A(3*idInstrument+0, idx)   =  posRef.y()/DEFAULT_R; //rz
            A(3*idInstrument+1, idx++) = -posRef.x()/DEFAULT_R; //rz
          }
        }

        Vector x;
        if(A.size())
        {
          Vector sigma;
          x = Vce::robustLeastSquares(A, l.column(idEpoch), 3, huber, huberPower, iterCount, sigma);
          l.column(idEpoch) -= A*x;
        }

        Vector h(8);
        h(0) = timeEpoch.mjd();
        UInt idx = 0;
        if(estimateShift)
        {
          h(1) = x(idx++);
          h(2) = x(idx++);
          h(3) = x(idx++);
        }
        if(estimateScale)
        {
          h(4) = x(idx++)/DEFAULT_R;
        }
        if(estimateRotation)
        {
          h(5) = x(idx++)/DEFAULT_R;
          h(6) = x(idx++)/DEFAULT_R;
          h(7) = x(idx++)/DEFAULT_R;
        }
        helmertTimeSeries.push_back(h);
      }

      // compute RMS values of current interval
      // -------------------------------------
      for(UInt idInstrument=0; idInstrument<countInstrument; idInstrument++)
        if(norm(l.row(3*idInstrument,3))>1e-10)
        {
          RmsPlotGrid(countOutput, 0) = time.mjd();
          RmsPlotGrid(countOutput, 1) = idInstrument;
          RmsPlotGrid(countOutput, 2) = factor * norm(l.row(3*idInstrument,3))/std::sqrt(3*l.columns());  // 1D RMS
          countOutput++;
        }
    }
    logTimerLoopEnd(countInterval);

    // ======================================================

    // save file
    // ---------
    if(!fileNameOutRmsPlotGrid.empty())
    {
      logStatus<<"Write RMS plot grid to file <"<<fileNameOutRmsPlotGrid<<">"<<Log::endl;
      writeFileMatrix(fileNameOutRmsPlotGrid, RmsPlotGrid.row(0, countOutput));
    }

    if(!fileNameHelmertTimeSeries.empty())
    {
      logStatus<<"Write Helmert time series to file <"<<fileNameHelmertTimeSeries<<">"<<Log::endl;
      Matrix H(helmertTimeSeries.size(), 8);
      for(UInt i=0; i<H.rows(); i++)
        copy(helmertTimeSeries.at(i).trans(), H.row(i));
      writeFileMatrix(fileNameHelmertTimeSeries, H);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
