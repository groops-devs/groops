/***********************************************/
/**
* @file instrumentDetrend.cpp
*
* @brief Reduces temporal parametrization (e.g. trend, polynomial) per arc from instrument file.
*
* @author Sebastian Strasser
* @date 2017-06-06
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Reduces \configClass{parametrizationTemporal}{parametrizationTemporalType} (e.g. const, trend, polynomial)
per arc from selected data columns of \configFile{inputfileInstrument}{instrument}
using a robust \reference{robust least squares adjustment}{fundamentals.robustLeastSquares}.

The \configFile{outputfileTimeSeriesArcParameters}{instrument} contains for every arc one (mid) epoch
with the estimated parameters. The order is: first all data (\verb|data0|, \verb|data1|, \ldots)
of first temporal parameter, followed by all data of the second temporal parameter and so on.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "misc/varianceComponentEstimation.h"

/***** CLASS ***********************************/

/** @brief Reduces temporal parametrization (e.g. trend, polynomial) from instrument file.
* @ingroup programsGroup */
class InstrumentDetrend
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentDetrend, PARALLEL, "Reduces temporal parametrization (e.g. trend, polynomial) per arc from instrument file.", Instrument)

/***********************************************/

void InstrumentDetrend::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName                   fileNameOut, fileNameParameters;
    FileName                   fileNameIn;
    ParametrizationTemporalPtr temporal;
    UInt                       startData, countData = MAX_UINT;
    Double                     huber, huberPower;
    UInt                       maxIter;

    renameDeprecatedConfig(config, "temporalRepresentation", "parametrizationTemporal", date2time(2020, 6, 3));

    readConfig(config, "outputfileInstrument",              fileNameOut,        Config::MUSTSET,  "",    "detrended instrument time series");
    readConfig(config, "outputfileTimeSeriesArcParameters", fileNameParameters, Config::OPTIONAL, "",    "time series of estimated parameters per arc");
    readConfig(config, "inputfileInstrument",               fileNameIn,         Config::MUSTSET,  "",    "");
    readConfig(config, "parametrizationTemporal",           temporal,           Config::MUSTSET,  "",    "per arc, data is reduced by temporal representation");
    readConfig(config, "startDataFields",                   startData,          Config::DEFAULT,  "0",   "start");
    readConfig(config, "countDataFields",                   countData,          Config::OPTIONAL, "",    "number of data fields (default: all after start)");
    readConfig(config, "huber",                             huber,              Config::DEFAULT,  "2.5", "for robust least squares");
    readConfig(config, "huberPower",                        huberPower,         Config::DEFAULT,  "1.5", "for robust least squares");
    readConfig(config, "huberMaxIteration",                 maxIter,            Config::DEFAULT,  "5",   "(maximum) number of iterations for robust estimation");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument files <"<<fileNameIn<<"> and detrending"<<Log::endl;
    InstrumentFile instrumentFile(fileNameIn);
    countData = std::min(countData, instrumentFile.dataCount(TRUE/*mustDefined*/)-startData);

    Matrix arcParameters;
    if(!fileNameParameters.empty())
      arcParameters = Matrix(instrumentFile.arcCount(), 1+countData*temporal->parameterCount());

    std::vector<Arc> arcList(instrumentFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      Arc arc = instrumentFile.readArc(arcNo);
      if(arc.size() == 0)
        return arc;

      // fill design matrix
      const std::vector<Time> times = arc.times();
      temporal->setInterval(times.front(), times.back()+medianSampling(times), TRUE);
      Matrix A(arc.size(), temporal->parameterCount());
      for(UInt idEpoch=0; idEpoch<arc.size(); idEpoch++)
        copy(temporal->factors(times.at(idEpoch)).trans(), A.row(idEpoch));

      Matrix l = arc.matrix();
      Vector sigma;
      Matrix x = Vce::robustLeastSquares(A, l.column(1+startData, countData), 1, huber, huberPower, maxIter, sigma);
      matMult(-1., A, x, l.column(1+startData, countData)); // residuals

      if(arcParameters.size())
      {
        arcParameters(arcNo, 0) = (0.5*(times.front()+times.back()+medianSampling(times))).mjd();
        copy(flatten(x.trans()).trans(), arcParameters.slice(arcNo, 1, 1, x.size()));
      }

      return Arc(times, l, arc.getType());
    }, comm); // forEach

    // ======================================================

    if(!fileNameOut.empty() && Parallel::isMaster(comm))
    {
      logStatus<<"write instrument data to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
      Arc::printStatistics(arcList);
    }

    if(!fileNameParameters.empty() && arcParameters.size())
    {
      Parallel::reduceSum(arcParameters, 0, comm);
      if(Parallel::isMaster(comm))
      {
        logStatus<<"write arc parameters to instrument file <"<<fileNameParameters<<">"<<Log::endl;
        InstrumentFile::write(fileNameParameters, Arc(arcParameters));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
