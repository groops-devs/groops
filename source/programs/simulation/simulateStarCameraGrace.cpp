/***********************************************/
/**
* @file simulateStarCameraGrace.cpp
*
* @brief Simulates the orientation of the two GRACE satellites.
*
* @author Torsten Mayer-Guerr
* @date 2009-08-02
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Simulates \file{star camera data}{instrument} of the two GRACE satellites.
\begin{itemize}
\item x: the antenna center pointing to the other satellite.
\item y: normal to line of sight and the radial direction.
\item z: forms a right handed system.
\end{itemize}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Simulates the orientation of the two GRACE satellites.
* @ingroup programsGroup */
class SimulateStarCameraGrace
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(SimulateStarCameraGrace, SINGLEPROCESS, "Simulates the orientation of the two GRACE satellites.", Simulation, Instrument)

/***********************************************/

void SimulateStarCameraGrace::run(Config &config)
{
  try
  {
    FileName    outStarCamera1Name, outStarCamera2Name;
    FileName    orbit1Name, orbit2Name;
    Vector3d    center1(1,0,0), center2(1,0,0);
    std::string choice;

    readConfig(config, "outputfileStarCamera1",   outStarCamera1Name, Config::MUSTSET,  "", "");
    readConfig(config, "outputfileStarCamera2",   outStarCamera2Name, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit1",         orbit1Name,         Config::MUSTSET,  "", "position define the orientation of the satellite at each epoch");
    readConfig(config, "inputfileOrbit2",         orbit2Name,         Config::MUSTSET,  "", "position define the orientation of the satellite at each epoch");
    if(readConfigChoice(config, "antennaCenters", choice,             Config::OPTIONAL, "", "KBR antenna phase center"))
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
    if(isCreateSchema(config)) return;

    const Rotary3d rotKFrame1 = inverse(Rotary3d(center1, Vector3d(0,1,0)));
    const Rotary3d rotKFrame2 = inverse(Rotary3d(center2, Vector3d(0,1,0)));

    logStatus<<"read orbit and generate star camera data"<<Log::endl;
    InstrumentFile orbit1File(orbit1Name);
    InstrumentFile orbit2File(orbit2Name);
    InstrumentFile::checkArcCount({orbit1File, orbit2File});

    std::list<Arc> arcList1, arcList2;
    for(UInt arc=0; arc<orbit1File.arcCount(); arc++)
    {
      OrbitArc orbit1 = orbit1File.readArc(arc);
      OrbitArc orbit2 = orbit2File.readArc(arc);
      Arc::checkSynchronized({orbit1, orbit2});

      StarCameraArc arc1, arc2;
      for(UInt i=0; i<orbit1.size(); i++)
      {
        // Line of sight
        const Vector3d e12 = normalize(orbit2.at(i).position - orbit1.at(i).position);

        StarCameraEpoch epoch1, epoch2;
        epoch1.time   = orbit1.at(i).time;
        epoch2.time   = orbit2.at(i).time;
        epoch1.rotary = Rotary3d( e12, crossProduct( e12, orbit1.at(i).position)) * rotKFrame1;
        epoch2.rotary = Rotary3d(-e12, crossProduct(-e12, orbit2.at(i).position)) * rotKFrame2;

        arc1.push_back(epoch1);
        arc2.push_back(epoch2);
      }
      arcList1.push_back(arc1);
      arcList2.push_back(arc2);
    }

    logStatus<<"write star camera data to file <"<<outStarCamera1Name<<">"<<Log::endl;
    InstrumentFile::write(outStarCamera1Name, arcList1);
    logStatus<<"write star camera data to file <"<<outStarCamera2Name<<">"<<Log::endl;
    InstrumentFile::write(outStarCamera2Name, arcList2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
