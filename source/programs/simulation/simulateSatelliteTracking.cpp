/***********************************************/
/**
* @file simulateSatelliteTracking.cpp
*
* @brief Simulate tracking data (range, range-rate, range-accelerations) between 2 satellites.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-21
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates \file{tracking data}{instrument} (range, range-rate, range-accelerations)
between 2 satellites. The range is given by
\begin{equation}
\rho(t) = \left\lVert{\M r_B(t) - \M r_A(t)}\right\rVert = \M e_{AB}(t)\cdot\M r_{AB}(t),
\end{equation}
with $\M r_{AB} = \M r_B - \M r_A$ and the unit vector in line of sight (LOS) direction
\begin{equation}\label{sst.los}
\M e_{AB} = \frac{\M r_{AB}}{\left\lVert{\M r_{AB}}\right\rVert}=\frac{\M r_{AB}}{\rho}.
\end{equation}
Range-rates~$\dot{\rho}$ and range accelrations~$\ddot{\rho}$ are obtained by differentation
\begin{equation}\label{obsRangeRate}
\dot{\rho}  = \M e_{AB}\cdot\dot{\M r}_{AB} + \dot{\M e}_{AB}\cdot\M r_{AB}
            = \M e_{AB}\cdot\dot{\M r}_{AB},
\end{equation}
\begin{equation}\label{obsRangeAccl}
\begin{split}
\ddot{\rho} &= \M e_{AB}\cdot\ddot{\M r}_{AB} +\dot{\M e}_{AB}\cdot\dot{\M r}_{AB}
            = \M e_{AB}\cdot\ddot{\M r}_{AB} +
   \frac{1}{\rho}\left(\dot{\M r}_{AB}^2-\dot{\rho}^2\right). \\
\end{split}
\end{equation}
with the derivative of the unit vector
\begin{equation}
\dot{\M e}_{AB}=\frac{d}{dt}\left(\frac{\M r_{AB}}{\rho}\right)
=\frac{\dot{\M r}_{AB}}{\rho}-\frac{\dot{\rho}\cdot\M r_{AB}}{\rho^2}
=\frac{1}{\rho}\left({\dot{\M r}_{AB}-\dot{\rho}\cdot\M e_{AB}}\right).
\end{equation}
The \configFile{inputfileOrbit}{instrument}s must contain positions, velocities, and acceleration
(see \program{OrbitAddVelocityAndAcceleration}).
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Simulate tracking data (range, range-rate, range-accelerations) between 2 satellites.
* @ingroup programsGroup */
class SimulateSatelliteTracking
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateSatelliteTracking, PARALLEL, "simulate tracking data (range, range-rate, range-accelerations) between 2 satellites", Simulation, Instrument)

/***********************************************/

void SimulateSatelliteTracking::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName outName, orbit1Name, orbit2Name;

    readConfig(config, "outputfileSatelliteTracking", outName,    Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit1",             orbit1Name, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit2",             orbit2Name, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // Tracking-Daten erzeugen
    // -----------------------
    logStatus<<"read orbits and generate tracking data"<<Log::endl;
    InstrumentFile orbit1File(orbit1Name);
    InstrumentFile orbit2File(orbit2Name);
    InstrumentFile::checkArcCount({orbit1File, orbit2File});
    std::vector<Arc> arcList(orbit1File.arcCount());

    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit1 = orbit1File.readArc(arcNo);
      OrbitArc orbit2 = orbit2File.readArc(arcNo);
      Arc::checkSynchronized({orbit1, orbit2});

      SatelliteTrackingArc arc;
      for(UInt i=0; i<orbit1.size(); i++)
      {
        Vector3d vrange         = orbit2.at(i).position     - orbit1.at(i).position;
        Vector3d vvelocity      = orbit2.at(i).velocity     - orbit1.at(i).velocity;
        Vector3d vacceleration  = orbit2.at(i).acceleration - orbit1.at(i).acceleration;

        SatelliteTrackingEpoch epoch;
        epoch.time              = orbit1.at(i).time;
        epoch.range             = vrange.norm();
        epoch.rangeRate         = inner(vrange,vvelocity)/epoch.range;
        epoch.rangeAcceleration = (vvelocity.quadsum()-pow(epoch.rangeRate,2))/epoch.range
                                + inner(vrange,vacceleration)/epoch.range;
        arc.push_back(epoch);
      }
      return arc;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write tracking data to file <"<<outName<<">"<<Log::endl;
      InstrumentFile::write(outName, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
