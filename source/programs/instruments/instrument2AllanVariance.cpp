/***********************************************/
/**
* @file instrument2AllanVariance.cpp
*
* @brief Compute Allan variance from instrument files.
*
* @author Andreas Kvas
* @date 2018-08-02
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the overlapping Allan variance from an
\configFile{inputfileInstrument}{instrument}.
The estimate is averaged over all arcs (arcs are assumed to contain no data gaps).

The overlapping Allan variance is defined as
\begin{equation}
  \sigma^2(m\tau_0) = \frac{1}{2(m\tau_0)^2(N-2m)} \sum_{n=1}^{N-2m}(x_{n+2m}-2x_{n+m}+x_n)^2,
\end{equation}
where $m\tau_0$ is the averaging interval defined by the median sampling $\tau_0$.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute Allan variance from instrument files.
* @ingroup programsGroup */
class Instrument2AllanVariance
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Instrument2AllanVariance, PARALLEL, "Compute Allan variance from instrument files.", Instrument, Covariance)
GROOPS_RENAMED_PROGRAM(InstrumentComputeAllanVariance, Instrument2AllanVariance, date2time(2020, 7, 7))

/***********************************************/

void Instrument2AllanVariance::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameInstrument, fileNameVariance;

    readConfig(config, "outputfileAllanVariance", fileNameVariance,   Config::MUSTSET, "", "column 0: averaging interval [seconds], column 1-(n-1): Allan variance for each data column");
    readConfig(config, "inputfileInstrument",     fileNameInstrument, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument data"<<Log::endl;
    InstrumentFile instrumentFile(fileNameInstrument);
    const UInt dataCount = instrumentFile.dataCount(TRUE/*mustDefined*/);

    UInt arcEpochCount, arcCount;
    Time sampling = seconds2time(1.0);
    if(Parallel::isMaster(comm))
    {
      arcCount = instrumentFile.arcCount();
      std::vector<Time> times;

      arcEpochCount = 0;
      for(UInt arcNo = 0; arcNo<arcCount; arcNo++)
      {
        Arc arc = instrumentFile.readArc(arcNo);
        if(arc.size() == 0)
          continue;
        auto arcTimes = arc.times();

        arcEpochCount = std::max(arc.size(), arcEpochCount);
        times.insert(times.end(), arcTimes.begin(), arcTimes.end());
      }
      sampling = medianSampling(times);

      logInfo<<"  maximum arc length: "<<arcEpochCount<<" epochs"<<Log::endl;
      logInfo<<"  median sampling:    "<<sampling.seconds()<<" seconds"<<Log::endl;
    }
    Parallel::broadCast(arcEpochCount, 0, comm);
    Parallel::broadCast(arcCount,      0, comm);

    logStatus<<"compute Allan variance"<<Log::endl;
    Matrix allanVariance(arcEpochCount/2, dataCount+1);
    Vector samplesPerInterval(allanVariance.rows());
    Parallel::forEach(arcCount, [&](UInt arcNo)
    {
      const Matrix data = instrumentFile.readArc(arcNo).matrix();
      for(UInt m=1; m<allanVariance.rows(); m++)
      {
        if(2*m >= data.rows())
          continue;

        const Double factor = 0.5/(std::pow(m*sampling.seconds(), 2) * (data.rows()-2*m));
        for(UInt i=1; i<data.columns(); i++)
          for(UInt n=0; n<data.rows()-2*m; n++)
          {
            allanVariance(m, i) += factor * std::pow(data(n+2*m, i) - 2*data(n+m, i) + data(n, i), 2);
            samplesPerInterval(m)++;
          }
      }
    }, comm);
    Parallel::reduceSum(allanVariance, 0, comm);
    Parallel::reduceSum(samplesPerInterval, 0, comm);

    if(Parallel::isMaster(comm))
    {
      for(UInt m=1; m<allanVariance.rows(); m++)
      {
        allanVariance.row(m) *= (1.0/samplesPerInterval(m));
        allanVariance(m, 0) = (m*sampling.seconds());
      }

      writeFileMatrix(fileNameVariance, allanVariance.row(1, allanVariance.rows()-1));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
