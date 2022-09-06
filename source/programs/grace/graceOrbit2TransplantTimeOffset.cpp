/***********************************************/
/**
* @file graceOrbit2TransplantTimeOffset.cpp
*
* @brief Compute the time shift between two co-orbiting satellites from their orbit data.
*
* @author Andreas Kvas
* @date 2022-04-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the time shift between two co-orbiting satellites based on dynamic orbit data.
When applied to data of the first satellite, the computed time shift virtually shifts data of first satellite into the location of the second satellite.
Note that \config{inputfileOrbit1} and \config{inputfileOrbit2} need velocity and acceleration data, which
can be computed with \program{OrbitAddVelocityAndAcceleration}.
The program tries to find a minimum of the objective function
\begin{equation}
  f(\Delta t) = \| r_1(t) - r_2(t + \Delta t) \|^2,
\end{equation}
by applying Newton's method to the first derivative, thus iteratively computing
\begin{equation}
  \Delta t_{k+1} = \Delta t_k + \frac{f'(\Delta t_k)}{f''(\Delta t_k)}.
\end{equation}
This iteration is stopped when the difference between to consecutive time shift values falls below \config{threshold} or
\config{maximumIterations} is reached. An \config{initialGuess} of the time shift can speed up convergence.

See also \program{OrbitAddVelocityAndAcceleration} and \program{InstrumentApplyTimeOffset}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute the time shift between two co-orbiting satellites from their orbit data.
* @ingroup programsGroup */
class GraceOrbit2TransplantTimeOffset
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceOrbit2TransplantTimeOffset, PARALLEL, "compute the time shift between two co-orbiting satellites from their orbit data", Grace, Instrument)

/***********************************************/

void GraceOrbit2TransplantTimeOffset::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameTimeOffset;
    FileName fileNameOrbit1,  fileNameOrbit2;
    UInt     interpolationDegree, maxIter;
    Double   initialGuess, threshold;

    readConfig(config, "outputfileTimeOffset", fileNameTimeOffset,  Config::MUSTSET,  "",     "estimated time offset in seconds (MISCVALUE)");
    readConfig(config, "inputfileOrbit1",      fileNameOrbit1,      Config::MUSTSET,  "",     "orbit data of satellite 1");
    readConfig(config, "inputfileOrbit2",      fileNameOrbit2,      Config::MUSTSET,  "",     "orbit data of satellite 2");
    readConfig(config, "interpolationDegree",  interpolationDegree, Config::DEFAULT,  "7",    "polynomial degree for the interpolation of position, velocity and acceleration");
    readConfig(config, "initialGuess",         initialGuess,        Config::DEFAULT,  "0.0",  "initial guess for the time shift [seconds]");
    readConfig(config, "maximumIterations",    maxIter,             Config::DEFAULT,  "50",   "maximum number of iterations");
    readConfig(config, "threshold",            threshold,           Config::DEFAULT,  "1e-5", "when the maximum difference between two iterations is below this value, stop [seconds]");
    if(isCreateSchema(config)) return;

    logStatus<<"read orbit data and estimate time shift"<<Log::endl;
    InstrumentFile  orbit1File(fileNameOrbit1);
    InstrumentFile  orbit2File(fileNameOrbit2);
    InstrumentFile::checkArcCount({orbit1File, orbit2File});

    std::vector<Arc> arcList(orbit1File.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit1 = orbit1File.readArc(arcNo);
      OrbitArc orbit2 = orbit2File.readArc(arcNo);
      Arc::checkSynchronized({orbit1, orbit2});

      const std::vector<Time> times2 = orbit2.times();
      Matrix dt(orbit2.size(), 2, initialGuess);

      Polynomial polynomial;
      polynomial.init(times2, interpolationDegree, FALSE/*throwException*/, FALSE/*isLeastSquares*/, -1.*(interpolationDegree+1), std::numeric_limits<Double>::infinity());

      Matrix A1 = orbit1.matrix();
      Matrix A2 = orbit2.matrix();

      if(maxabs(A1.column(4, 3)) == 0.0 || maxabs(A1.column(7, 3)) == 0.0 || maxabs(A2.column(4, 3)) == 0.0 || maxabs(A2.column(7, 3)) == 0.0)
        throw(Exception("dynamic orbits require both velocity and acceleration"));

      UInt iter;
      for(iter=0; iter<maxIter; iter++)
      {
        std::vector<Time> timesNew = times2;
        for(UInt k=0; k<timesNew.size(); k++)
          timesNew[k] += seconds2time(dt(k, 1));

        Vector dtk(dt.rows());
        Matrix A2_interp = polynomial.interpolate(timesNew, A2.column(1, A2.columns()-1));
        for(UInt k=0; k<timesNew.size(); k++)
        {
          Vector3d dr((A1.slice(k, 1, 1, 3) - A2_interp.slice(k, 0, 1, 3)).trans());
          Vector3d v2(A2_interp.slice(k, 3, 1, 3).trans());
          Vector3d a2(A2_interp.slice(k, 6, 1, 3).trans());

          dtk(k) = inner(dr, v2) / (inner(dr, a2) + inner(v2, v2));
          dt(k, 1) += dtk(k);
        }
        if(maxabs(dtk) < threshold)
          break;
      }
      if(iter == maxIter)
        logWarning<<"maximum iterations reached"<<Log::endl;

      return Arc(times2, dt, Epoch::Type::MISCVALUE);
    }, comm); // forEach

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write time offset to file <"<<fileNameTimeOffset<<">"<<Log::endl;
      InstrumentFile::write(fileNameTimeOffset, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
