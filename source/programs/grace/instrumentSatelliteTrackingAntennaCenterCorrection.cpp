/***********************************************/
/**
* @file instrumentSatelliteTrackingAntennaCenterCorrection.cpp
*
* @brief Compute antenna center correction from orbit configuration.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-21
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the correction due to offset of the antenna center relative the center of mass.
The offsets $\M c_A$ and $\M c_B$ in \configFile{inputfileAntennaCenters}{matrix} are given in the satellite
reference frame. These offsets are rotated into the the inertial frame with $\M D_A$ and $\M D_B$ from
\configFile{inputfileStarCamera}{instrument} and projected onto the line of sight (LOS)
\begin{equation}
  \rho_{AOC} = \M e_{AB}\cdot(\M D_A\,\M c_A - \M D_B\,\M c_B),
\end{equation}
with the unit vector in line of sight direction
\begin{equation}
  \M e_{AB} = \frac{\M r_B - \M r_A}{\left\lVert{\M r_B - \M r_A}\right\rVert}.
\end{equation}
The corrections for the range-rates and range-acceleration are computed by differentiating
an interpolation polynomial of degree \config{interpolationDegree}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute antenna center correction from orbit configuration.
* @ingroup programsGroup */
class InstrumentSatelliteTrackingAntennaCenterCorrection
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentSatelliteTrackingAntennaCenterCorrection, PARALLEL, "compute antenna center correction from orbit configuration", Grace, Instrument)

/***********************************************/

void InstrumentSatelliteTrackingAntennaCenterCorrection::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName    outSSTName;
    FileName    orbit1Name,  starCamera1Name;
    FileName    orbit2Name,  starCamera2Name;
    Vector3d    center1, center2;
    UInt        degree;
    std::string choice;

    readConfig(config, "outputfileSatelliteTracking", outSSTName,      Config::OPTIONAL, "", "corrections for range, range-rate, and range-accelerations");
    readConfig(config, "inputfileOrbit1",             orbit1Name,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit2",             orbit2Name,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera1",        starCamera1Name, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera2",        starCamera2Name, Config::MUSTSET,  "", "");
    if(readConfigChoice(config, "antennaCenters",     choice,          Config::MUSTSET, "", "KBR antenna phase center"))
    {
      if(readConfigChoiceElement(config, "value", choice, ""))
      {
        readConfig(config, "center1X", center1.x(), Config::DEFAULT,   "1.4451172588", "x-coordinate of antenna position in SRF [m] for GRACEA");
        readConfig(config, "center1Y", center1.y(), Config::DEFAULT,  "-0.0004233040", "y-coordinate of antenna position in SRF [m] for GRACEA");
        readConfig(config, "center1Z", center1.z(), Config::DEFAULT,   "0.0022786600", "z-coordinate of antenna position in SRF [m] for GRACEA");
        readConfig(config, "center2X", center2.x(), Config::DEFAULT,   "1.4443870350", "x-coordinate of antenna position in SRF [m] for GRACEB");
        readConfig(config, "center2Y", center2.y(), Config::DEFAULT,   "0.0005761203", "y-coordinate of antenna position in SRF [m] for GRACEB");
        readConfig(config, "center2Z", center2.z(), Config::DEFAULT,   "0.0033040887", "z-coordinate of antenna position in SRF [m] for GRACEB");
      }
      if(readConfigChoiceElement(config, "file",  choice, ""))
      {
        FileName fileName;
        readConfig(config, "inputAntennaCenters", fileName, Config::MUSTSET, "", "");
        if(!isCreateSchema(config))
        {
          Matrix x;
          readFileMatrix(fileName, x);
          center1 = Vector3d(x(0,0), x(1,0), x(2,0));
          center2 = Vector3d(x(3,0), x(4,0), x(5,0));
        }
      }
      endChoice(config);
    }
    readConfig(config, "interpolationDegree", degree, Config::DEFAULT,  "2", "differentiation by polynomial approximation of degree n");
    if(isCreateSchema(config)) return;

    logStatus<<"read orbit and star camera data and generate antenna offset corrections"<<Log::endl;
    InstrumentFile  orbit1File(orbit1Name);
    InstrumentFile  orbit2File(orbit2Name);
    InstrumentFile  starCamera1File(starCamera1Name);
    InstrumentFile  starCamera2File(starCamera2Name);
    InstrumentFile::checkArcCount({orbit1File, orbit2File, starCamera1File, starCamera2File});

    std::vector<Arc> arcList(orbit1File.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc      orbit1      = orbit1File.readArc(arcNo);
      OrbitArc      orbit2      = orbit2File.readArc(arcNo);
      StarCameraArc starCamera1 = starCamera1File.readArc(arcNo);
      StarCameraArc starCamera2 = starCamera2File.readArc(arcNo);
      Arc::checkSynchronized({orbit1, orbit2, starCamera1, starCamera2});

      Matrix A(orbit1.size(), 4);
      for(UInt i=0; i<A.rows(); i++)
      {
        const Vector3d u = (orbit2.at(i).position - orbit1.at(i).position);                                       // center of mass vector
        const Vector3d v = (starCamera2.at(i).rotary.rotate(center2) - starCamera1.at(i).rotary.rotate(center1)); // combined antenna offset vector
        A(i,1) = u.r() - (u+v).r();  // COM - ANT -> AOC
      }

      const std::vector<Time> times = orbit1.times();
      Polynomial p(times, degree);
      copy(p.derivative(   times, A.column(1)), A.column(2));
      copy(p.derivative2nd(times, A.column(1)), A.column(3));

      return Arc(times, A, Epoch::Type::SATELLITETRACKING);
    }, comm); // forEach

    if(Parallel::isMaster(comm) && !outSSTName.empty())
    {
      logStatus<<"write tracking data to file <"<<outSSTName<<">"<<Log::endl;
      InstrumentFile::write(outSSTName, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
