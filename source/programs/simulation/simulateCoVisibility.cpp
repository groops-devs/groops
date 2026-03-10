/***********************************************/
/**
* @file simulateSatelliteTracking.cpp
*
* @brief Simulate co-visibility from a ground station of 2 satellites.
*
* @author Schievano Giulia
* @date 2025-01-09
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
#include "classes/earthRotation/earthRotation.h"
#include "base/ellipsoid.h"

/***** CLASS ***********************************/

/** @brief Simulate tracking data (range, range-rate, range-accelerations) between 2 satellites.
* @ingroup programsGroup */
class SimulateCoVisibility
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateCoVisibility, PARALLEL, "simulate tracking data (range, range-rate, range-accelerations) between 2 satellites", Simulation, Instrument)

/***********************************************/

void SimulateCoVisibility::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName outName1, orbit1Name, orbit2Name;
    Vector3d OGSpos;
    Double minElev, maxElev, maxRange;
    EarthRotationPtr earthRotation;

    readConfig(config, "outputfileSatelliteTracking",  outName1,        Config::MUSTSET,  "file", "GROOPS co-visibility file");
    readConfig(config, "inputfileOrbit1",              orbit1Name,      Config::MUSTSET,  "file", "GROOPS orbit file");
    readConfig(config, "inputfileOrbit2",              orbit2Name,      Config::MUSTSET,  "file", "GROOPS orbit file");
    readConfig(config, "position0x",                   OGSpos.x(),      Config::MUSTSET, "", "[m] in ECEF");
    readConfig(config, "position0y",                   OGSpos.y(),      Config::MUSTSET, "", "[m] in ECEF");
    readConfig(config, "position0z",                   OGSpos.z(),      Config::MUSTSET, "", "[m] in ECEF");
    readConfig(config, "minElev",                      minElev,         Config::MUSTSET, "", "[deg] mask on horizon");
    readConfig(config, "maxElev",                      maxElev,         Config::MUSTSET, "", "[deg] mask on zenith");
    readConfig(config, "maxRange",                     maxRange,        Config::MUSTSET, "", "[m] mask on SAT-OGS range");
    readConfig(config, "earthRotation",                earthRotation,   Config::MUSTSET, "file", "rotate data into Earth-fixed frame");

    if(isCreateSchema(config)) return;

    // Tracking-Daten erzeugen
    // -----------------------
    logStatus<<"read orbits and generate tracking data"<<Log::endl;
    InstrumentFile orbit1File(orbit1Name);
    InstrumentFile orbit2File(orbit2Name);
    InstrumentFile::checkArcCount({orbit1File, orbit2File});
    std::vector<Arc> arcList(orbit1File.arcCount());

    Double theta_th_min = (PI/2-minElev*DEG2RAD)*RAD2DEG;
    Double theta_th_max = (PI/2-maxElev*DEG2RAD)*RAD2DEG;

    logInfo<<"  min angle           : "<<theta_th_min<<" degrees"<<Log::endl;
    logInfo<<"  max angle           : "<<theta_th_max<<" degrees"<<Log::endl;
 

    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit1 = orbit1File.readArc(arcNo);
      OrbitArc orbit2 = orbit2File.readArc(arcNo);
      Arc::checkSynchronized({orbit1, orbit2});


      UInt iter = 0;
      Double sampling = orbit1.at(1).time.seconds() - orbit1.at(0).time.seconds();
      

      SimulateCoVisibilityArc arc;
      for(UInt i=0; i<orbit1.size(); i++)
      
      { 

      Rotary3d rot;
      Vector3d omega;
      if(earthRotation)
      {
        rot   = earthRotation->rotaryMatrix(orbit1.at(i).time);

      }
       Vector3d Orbit1_pos          = rot.rotate(orbit1.at(i).position);
       Vector3d Orbit2_pos          = rot.rotate(orbit2.at(i).position);

       Vector3d Orbit1_vel          = rot.rotate(orbit1.at(i).velocity);
       Vector3d Orbit2_vel          = rot.rotate(orbit2.at(i).velocity);

       Vector3d Orbit1_acc          = rot.rotate(orbit1.at(i).acceleration);
       Vector3d Orbit2_acc          = rot.rotate(orbit2.at(i).acceleration);
     
     

       Vector3d vrange          = Orbit2_pos     - Orbit1_pos;
       Vector3d vvelocity       = Orbit2_vel     - Orbit1_vel;
       Vector3d vacceleration   = Orbit2_acc     - Orbit1_acc;

       Vector3d los_sat         = vrange/vrange.norm();
       Vector3d vrangeOGS_s1    = Orbit1_pos     - OGSpos;
       Vector3d vrangeOGS_s2    = Orbit2_pos     - OGSpos;
       Vector3d los_OGS_s1      = vrangeOGS_s1/vrangeOGS_s1.norm();
       Vector3d los_OGS_s2      = vrangeOGS_s2/vrangeOGS_s2.norm();

      Double L = 48.081333*DEG2RAD; // Latitude
      Double B = 11.283*DEG2RAD; // Longitude
      Double H = 400; // Altitude

      Matrix rMatrix(3,3);
      
      rMatrix(0,0) =  -std::sin(B); 
      rMatrix(0,1) =   std::cos(B);
      rMatrix(0,2) =   0;
      
      rMatrix(1,0) =  -std::cos(B)*std::sin(L);
      rMatrix(1,1) =  -std::sin(B)*std::sin(L);
      rMatrix(1,2) =   std::cos(L);

      rMatrix(2,0) =   std::cos(B)*std::cos(L);
      rMatrix(2,1) =   std::sin(B)*std::cos(L);
      rMatrix(2,2) =   std::sin(L);


      
      SimulateCoVisibilityEpoch epoch;
      epoch.time               = orbit1.at(i).time;
      epoch.range              = vrange.norm();
      epoch.rangeRateProj      = inner(vvelocity, los_sat);
      epoch.rangeRate          = vvelocity.norm();
      epoch.rangeAcceleration  = vacceleration.norm();
      


     // Double relAzimuth1 = std::atan2(inner(Vector3d(rMatrix(0,0), rMatrix(0,1), rMatrix(0,2)) ,los_OGS_s1), 
     //                                 inner(Vector3d(rMatrix(1,0), rMatrix(1,1), rMatrix(1,2)), los_OGS_s1));
      epoch.relAzimuth1 = std::atan2(inner(Vector3d(rMatrix(0,0), rMatrix(0,1), rMatrix(0,2)) ,los_OGS_s1), 
                                     inner(Vector3d(rMatrix(1,0), rMatrix(1,1), rMatrix(1,2)), los_OGS_s1))*RAD2DEG; // geocentric to ENU 

     // Double relAzimuth2 = std::atan2(inner(Vector3d(rMatrix(0,0), rMatrix(0,1), rMatrix(0,2)) ,los_OGS_s2), 
     //                                 inner(Vector3d(rMatrix(1,0), rMatrix(1,1), rMatrix(1,2)), los_OGS_s2));
      epoch.relAzimuth2 = std::atan2(inner(Vector3d(rMatrix(0,0), rMatrix(0,1), rMatrix(0,2)) ,los_OGS_s2), 
                                     inner(Vector3d(rMatrix(1,0), rMatrix(1,1), rMatrix(1,2)), los_OGS_s2))*RAD2DEG;

      
       // epoch.relElevation      = PI/2 -s
       //                            std::acos(inner(vrange, orbit1.at(i).position)/(vrange.norm(), 
       //                            orbit1.at(i).position.norm()));   

        epoch.relElevation1      = (std::acos(inner(OGSpos, vrangeOGS_s1)/(vrangeOGS_s1.norm()* 
                                    OGSpos.norm())))*RAD2DEG;

        
        
        if (epoch.relElevation1 < theta_th_min & epoch.relElevation1 > theta_th_max & vrangeOGS_s1.norm() < maxRange) {
          // block of code if condition is true
          epoch.coVis1  = 1;
        }
        else {
          // block of code if condition is false
          epoch.coVis1 = 0;
        }

        
        epoch.relElevation2      = (std::acos(inner(OGSpos, vrangeOGS_s2)/(vrangeOGS_s2.norm()* 
                                    OGSpos.norm())))*RAD2DEG;


        if (epoch.relElevation2 < theta_th_min & epoch.relElevation2 > theta_th_max & vrangeOGS_s2.norm() < maxRange) {
          // block of code if condition is true
          epoch.coVis2  = 1;
        }
        else {
          // block of code if condition is false
          epoch.coVis2 = 0;
        }


        if (epoch.coVis1  == 1 & epoch.coVis2  == 1) {
          // block of code if condition is true
          iter = iter + 1;
          epoch.coVis  = 1;
          epoch.coVisTot  = iter *sampling;
        }
        else {
          // block of code if condition is false
          epoch.coVis = 0;
          epoch.coVisTot = iter *sampling;
        }

       // 
        epoch.OGS_s1 = vrangeOGS_s1.norm();
        epoch.OGS_s2 = vrangeOGS_s2.norm();

        epoch.rangeRateProj_OGS_s1      = inner(Orbit1_vel, los_OGS_s1);
        

        arc.push_back(epoch);

        // Info
        // ----
      
        logInfo<<"Tracking"<<Log::endl;
        logInfo<<"  range           : "<<epoch.range<<"m"<<Log::endl;
        logInfo<<"  range rate      : "<<epoch.rangeRate<<"m/s"<<Log::endl;
        logInfo<<"  rel elev1       : "<<epoch.relElevation1<<" Degree"<<Log::endl;
        logInfo<<"  rel elev2       : "<<epoch.relElevation2<<" Degree"<<Log::endl;
        logInfo<<"  rel azim1       : "<<epoch.relAzimuth1<<" Degree"<<Log::endl;
        logInfo<<"  rel azim2       : "<<epoch.relAzimuth2<<" Degree"<<Log::endl;
        logInfo<<"  co-visibility   : "<<epoch.coVis<<"  "<<Log::endl;
        logInfo<<"  co-visibility   : "<<epoch.coVisTot<<" seconds"<<Log::endl;
        logInfo<<"  OGS Sat1 range  : "<<epoch.OGS_s1<<" km"<<Log::endl;
        logInfo<<"  OGS Sat2 range  : "<<epoch.OGS_s2<<" km"<<Log::endl;
        logInfo<<"  Projected vel   : "<<epoch.rangeRateProj<<" m/s"<<Log::endl;
        logInfo<<"  Projected v OGS : "<<epoch.rangeRateProj_OGS_s1<<" m/s"<<Log::endl;
        logInfo<<"  Posx sat1  : "<<Orbit1_pos.x()<<" km"<<Log::endl;
        logInfo<<"  Posy sat1  : "<<Orbit1_pos.y()<<" km"<<Log::endl;
        logInfo<<"  Posz sat1  : "<<Orbit1_pos.z()<<" km"<<Log::endl;
        logInfo<<"  Posx sat2  : "<<Orbit2_pos.x()<<" km"<<Log::endl;
        logInfo<<"  Posy sat2  : "<<Orbit2_pos.y()<<" km"<<Log::endl;
        logInfo<<"  Posz sat2  : "<<Orbit2_pos.z()<<" km"<<Log::endl;





      }
      return arc;
    }, comm);



    if(Parallel::isMaster(comm))
    {
      logStatus<<"write tracking data to file <"<<outName1<<">"<<Log::endl;
      InstrumentFile::write(outName1, arcList);
      Arc::printStatistics(arcList);
    }

    

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
