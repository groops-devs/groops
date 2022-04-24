/***********************************************/
/**
* @file simulateStarCameraSentinel1.cpp
*
* @brief simulate star camera data for Sentinel 1A.
*
* @author Norbert Zehentner
* @date 2016-03-23
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates \file{star camera}{instrument} measurements at each satellite's position for the Sentinel 1A satellite.
The \configFile{inputfileOrbit}{instrument} must contain positions and velocities (see \program{OrbitAddVelocityAndAcceleration}).
The resulting rotation matrices rotate from satellite frame to inertial frame.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#ifndef GROOPS_DISABLE_ERFA
#include <erfa.h>
#endif

/***** CLASS ***********************************/

/** @brief simulate star camera data for Sentinel 1A.
* @ingroup programsGroup */
class SimulateStarCameraSentinel1
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateStarCameraSentinel1, PARALLEL, "simulate star camera data for Sentinel 1", Simulation, Instrument)

/***********************************************/

void SimulateStarCameraSentinel1::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName orbitName, starCameraName;

    readConfig(config, "outputfileStarCamera", starCameraName, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit",       orbitName,      Config::MUSTSET, "", "position and velocity defines the orientation of the satellite at each epoch");
    if(isCreateSchema(config)) return;

    // StarCamera-Daten erzeugen
    // -----------------------
    logStatus<<"read orbit and generate star camera data"<<Log::endl;
    Ellipsoid      ellipsoid(DEFAULT_GRS80_a, DEFAULT_GRS80_f);
    InstrumentFile orbitFile(orbitName);

    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit = orbitFile.readArc(arcNo);
      StarCameraArc arc;
      for(UInt i=0; i<orbit.size(); i++)
      {
        // get rotation form inertial to TOD frame
        Double rc2i[3][3];
#ifdef GROOPS_DISABLE_ERFA
        throw(Exception("compiled without ERFA library"));
#else
        const Time timeTT = timeGPS2TT(orbit.at(i).time);
        Double X=0, Y=0, S=0;
        eraXys00b(2400000.5+timeTT.mjdInt(), timeTT.mjdMod(), &X, &Y, &S);
        eraC2ixys(X, Y, S, rc2i);
#endif
        Matrix R(3,3);
        R(0,0) = rc2i[0][0]; R(0,1) = rc2i[1][0]; R(0,2) = rc2i[2][0];
        R(1,0) = rc2i[0][1]; R(1,1) = rc2i[1][1]; R(1,2) = rc2i[2][1];
        R(2,0) = rc2i[0][2]; R(2,1) = rc2i[1][2]; R(2,2) = rc2i[2][2];
        Rotary3d rot(R);

        Vector3d pos = orbit.at(i).position;
        Vector3d vel = orbit.at(i).velocity;
        const Double omega = 0.729211585e-4;
        vel.x() +=  omega*pos.y(); // Correct velocity for earth rotation
        vel.y() += -omega*pos.x();
        // correct position
        // compute subsatellite point
        Angle  L, B;
        Double h;
        ellipsoid(pos, L, B, h); // FIXME rot.rotate(pos)???
        Vector3d posGround = ellipsoid(L, B, 0);
        Vector3d x         = normalize(vel);
        Rotary3d rotC2ZD   = Rotary3d(x, crossProduct(x, pos-posGround)) * rot;

        // part 2: rotation from zero-Doppler reference frame to sentinel 1 reference frame (body)
        ellipsoid(rot.rotate(pos), L, B, h); // compute height
        const Double   rThetaRef  = 0.513999;  // rad
        const Double   rAlphaRoll = 9.878564e-7; // rad/meter
        const Double   rHref      = 711700;
        const Rotary3d rot2       = rotaryX(Angle(rThetaRef - rAlphaRoll*(h-rHref)));

        StarCameraEpoch epoch;
        epoch.time   = orbit.at(i).time;
        epoch.rotary = rotC2ZD*rot2;
        arc.push_back(epoch);
      }
      return arc;
    }, comm); // forEach

    // write
    // -----
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write star camera data to file <"<<starCameraName<<">"<<Log::endl;
      InstrumentFile::write(starCameraName, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
