/***********************************************/
/**
* @file graceSstSpecialEvents.cpp
*
* @brief Detect GRACE events (eclipse, sun intrusion into star camera box).
*
* @author Saniya Behzadpour
* @date 2018-08-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Time-indexing deterministic signals in the GRACE K-Band measurements caused by Sun intrusions
into the star camera baffles of GRACE-A and eclipse transits of the satellites.
The events are determined by satellites' position (\configFile{inputfileOrbit1/2}{instrument})
and orientation (\configFile{inputfileStarCamera1/2}{instrument}). Each type of event is represented
by its mid-interval point per orbit revolution and is reported in \configFile{outputfileEvents}{instrument}.

The waveform of each event is nearly constant within one month and can be approximated by a polynomial.
For the purpose of gravity field recovery, each waveform is parameterized by a polynomial and the coefficients
of this polynomial are estimated as additional instrument calibration parameters in a common adjustment
with all other instrument, satellite, and gravity field parameters,
see \configClass{parametrizationSatelliteTracking:specialEffect}{parametrizationSatelliteTrackingType:specialEffect}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/eclipse/eclipse.h"

/***** CLASS ***********************************/

/** @brief Detect GRACE events (eclipse, sun intrusion into star camera box).
* @ingroup programsGroup */
class GraceSstSpecialEvents
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceSstSpecialEvents, SINGLEPROCESS, "Detect GRACE events (eclipse, sun intrusion into star camera box).", Grace)

/***********************************************/

void GraceSstSpecialEvents::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName       fileNameOut, fileNameOutInterval;
    FileName       fileNameOrbit1, fileNameOrbit2;
    FileName       fileNameStarCamera1, fileNameStarCamera2;
    EphemeridesPtr ephemerides;
    EclipsePtr     eclipse;
    Double         marginBefore, marginAfter;

    readConfig(config, "outputfileEvents",     fileNameOut,         Config::MUSTSET,  "", "");
    readConfig(config, "outputfileIntervals",  fileNameOutInterval, Config::OPTIONAL, "", "");
    readConfig(config, "inputfileOrbit1",      fileNameOrbit1,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit2",      fileNameOrbit2,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera1", fileNameStarCamera1, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera2", fileNameStarCamera2, Config::MUSTSET,  "", "");
    readConfig(config, "ephemerides",          ephemerides,         Config::MUSTSET,  "", "");
    readConfig(config, "eclipse",              eclipse,             Config::MUSTSET,  "", "");
    readConfig(config, "marginLeft",           marginBefore,        Config::DEFAULT,  "100", "margin size (on both sides) [seconds]");
    readConfig(config, "marginRight",          marginAfter,         Config::DEFAULT,  "100", "margin size (on both sides) [seconds]");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read star camera and orbit data"<<Log::endl;
    OrbitArc      orbit1      = InstrumentFile::read(fileNameOrbit1);
    OrbitArc      orbit2      = InstrumentFile::read(fileNameOrbit2);
    StarCameraArc starCamera1 = InstrumentFile::read(fileNameStarCamera1);
    Arc::checkSynchronized({orbit1, orbit2, starCamera1});

    MiscValueArc arc;
    for(UInt i=0; i<orbit1.size(); i++)
    {
      // eclipse transit
      // ---------------
      const Double delta = eclipse->factor(orbit2.at(i).time, orbit2.at(i).position, ephemerides)
                         - eclipse->factor(orbit1.at(i).time, orbit1.at(i).position, ephemerides);

      if(fabs(delta) > 1e-10)
      {
        MiscValueEpoch epoch;
        epoch.time  = orbit1.at(i).time;
        epoch.value = (delta>0) ? 1 : 2;
        arc.push_back(epoch);
        continue;
      }

      // star camera box
      // ----------------
      const Vector3d direction = normalize(orbit1.at(i).position - ephemerides->position(orbit1.at(i).time, Ephemerides::SUN)); // from sun to satellite
      const Vector3d normP     = rotaryZ(Angle(180*DEG2RAD)).rotate(rotaryX(Angle(+135*DEG2RAD)).rotate(starCamera1.at(i).rotary.inverseRotate(direction))); // Primary
      const Vector3d normS     = rotaryZ(Angle(180*DEG2RAD)).rotate(rotaryX(Angle(-135*DEG2RAD)).rotate(starCamera1.at(i).rotary.inverseRotate(direction))); // Secondary
      const Vector3d norm      = (normP.z() < normS.z()) ? normP : normS;
      const Double foreAft     = RAD2DEG * std::atan(norm.x()/norm.z());
      const Double upDown      = RAD2DEG * std::atan(norm.y()/norm.z());
      const Double width       = 0.6;
      if((std::fabs(foreAft) < 10+width) && (std::fabs(upDown) < 8.5+width) && ((std::fabs(std::fabs(foreAft)-10) < width) || (std::fabs(std::fabs(upDown)-8.5) < width)))
      {
        MiscValueEpoch epoch;
        epoch.time  = orbit1.at(i).time;
        if (std::fabs(std::fabs(upDown)-8.5) < width)
        {
          if (upDown>0)
            epoch.value = ((foreAft>0) ? 5 : 6);
          else
            epoch.value = ((foreAft>0) ? 7 : 8);
        }
        else
          epoch.value = ((foreAft>0) ? 3 : 4);
        arc.push_back(epoch);
      }
    }

    // =============================================

    // Replace by mid epoch
    // --------------------
    MiscValueArc arcNew;
    Double       type = -1;
    Time         timeStart;
    for(UInt i=0; i<arc.size(); i++)
    {
      if((arc.at(i).value != type) || (arc.at(i).time > arc.at(i-1).time+seconds2time(10*60)))
      {
        arcNew.push_back(arc.at(i));
        timeStart = arc.at(i).time;
        type      = arc.at(i).value;
      }
      if(arc.at(i).value == type)
        arcNew.at(arcNew.size()-1).time = timeStart + 0.5*(arc.at(i).time-timeStart);
    }

    logStatus<<"write events to file<"<<fileNameOut<<">"<<Log::endl;
    InstrumentFile::write(fileNameOut, arcNew);
    Arc::printStatistics(arcNew);

    // =============================================

    MiscValueArc arcInterval;
    std::vector<Time> time = orbit1.times();
    UInt idx = 0;
    for(UInt idEpoch=0; idEpoch<time.size(); idEpoch++)
    {
      // find event
      while((idx < arcNew.size()) && (arcNew.at(idx).time+seconds2time(marginAfter) < time.at(idEpoch)))
        idx++;
      if(idx >= arcNew.size())
        break;
      if(time.at(idEpoch) < arcNew.at(idx).time-seconds2time(marginBefore))
        continue;

      MiscValueEpoch epoch;
      epoch.time  = time.at(idEpoch);
      epoch.value = arcNew.at(idx).value;
      arcInterval.push_back(epoch);
    }//for

    logStatus<<"write intervals to file<"<<fileNameOutInterval<<">"<<Log::endl;
    InstrumentFile::write(fileNameOutInterval, arcInterval);
    Arc::printStatistics(arcInterval);

    // =============================================
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
