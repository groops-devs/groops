/***********************************************/
/**
* @file simulateStarCameraGnss.cpp
*
* @brief Simulate star camera data for GNSS satellites using nominal orientation or an attitude model.
*
* @author Sebastian Strasser
* @date 2016-07-14
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates \file{star camera}{instrument} measurements at each satellite position.
The resulting rotation matrices rotate from body frame to inertial frame.

The nominal orientation is simulated to be z-axis towards the center of Earth, y-axis is along solar panel axis,
and x-axis forms a right hand system (points roughly towards the sun).

GPS attitude model: Kouba (2009), 'A simplified yaw-attitude model for eclipsing GPS satellites' \url{https://doi.org/10.1007/s10291-008-0092-1}

Galileo attitude model: \url{https://www.gsc-europa.eu/support-to-developers/galileo-satellite-metadata}

See also \program{GnssAntex2AntennaDefinition} and \program{GnssYawBias2Instrument}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/kepler.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileGnssStationInfo.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/eclipse/eclipse.h"

/***** CLASS ***********************************/

/** @brief Simulate star camera data for GNSS satellites using nominal orientation or an attitude model.
* @ingroup programsGroup */
class SimulateStarCameraGnss
{
public:
  void run(Config &config);

private:
  Matrix yawRates;

  /// Helper class used for computation of true satellite attitude
  class Epoch
  {
  public:
    Time     time;
    Double   yawNominal;
    Double   yawTrue;
    Double   yawRateNominal;
    Double   yawRateTrue;
    Double   yawBias;
    Double   mu;
    Double   beta;
    Double   shadowFactor;
    Vector3d pos;
    Vector3d vel;
    Vector3d posSun;

    Epoch() : yawNominal(0), yawTrue(0), yawRateNominal(0), yawRateTrue(0), yawBias(0), mu(0), beta(0), shadowFactor(0) {}
  };

  /// GNSS block number
  enum GnssBlock
  {
    UNKNOWN    = 0,
    GPS_I      = 1,
    GPS_II     = 2,
    GPS_IIA    = 3,
    GPS_IIR    = 4,
    GPS_IIRA   = 5,
    GPS_IIRB   = 6,
    GPS_IIRM   = 7,
    GPS_IIF    = 8,
    GPS_IIIA   = 9,
    GLONASS    = 101,
    GLONASS_M  = 102,
    GLONASS_K1 = 103,
    GALILEO_1  = 201,
    GALILEO_2  = 202
//    COMPASS   = 301
  };

  /** @brief Convert GNSS block from string to enum. */
  static GnssBlock string2block(const std::string &string)
  {
    if(string == "BLOCK I")
      return GPS_I;
    else if(string == "BLOCK II")
      return GPS_II;
    else if(string == "BLOCK IIA")
      return GPS_IIA;
    else if(string == "BLOCK IIR-A")
      return GPS_IIRA;
    else if(string == "BLOCK IIR-B")
      return GPS_IIRB;
    else if(string == "BLOCK IIR-M")
      return GPS_IIRM;
    else if(string == "BLOCK IIF")
      return GPS_IIF;
    else if(string == "BLOCK IIIA")
      return GPS_IIIA;
    else if(string == "GLONASS")
      return GLONASS;
    else if(string == "GLONASS-M")
      return GLONASS_M;
    else if(string == "GLONASS-K1")
      return GLONASS_K1;
    else if(string == "GALILEO-1")
      return GALILEO_1;
    else if(string == "GALILEO-2")
      return GALILEO_2;
    else
      throw(Exception("cannot convert string to GnssBlock: " + string));
  }

  /** @brief Convert GNSS block from enum to string. */
  static std::string block2string(GnssBlock block)
  {
    if(block == GPS_I)
      return "BLOCK I";
    else if(block == GPS_II)
      return "BLOCK II";
    else if(block == GPS_IIA)
      return "BLOCK IIA";
    else if(block == GPS_IIRA)
      return "BLOCK IIR-A";
    else if(block == GPS_IIRB)
      return "BLOCK IIR-B";
    else if(block == GPS_IIRM)
      return "BLOCK IIR-M";
    else if(block == GPS_IIF)
      return "BLOCK IIF";
    else if(block == GPS_IIIA)
      return "BLOCK IIIA";
    else if(block == GLONASS)
      return "GLONASS";
    else if(block == GLONASS_M)
      return "GLONASS-M";
    else if(block == GLONASS_K1)
      return "GLONASS-K1";
    else if(block == GALILEO_1)
      return "GALILEO-1";
    else if(block == GALILEO_2)
      return "GALILEO-2";
    else
      throw(Exception("cannot convert GnssBlock to string: " + block%"%i"s));
  }

  /** @brief Computes nominal satellite orientation.
   *  Nominal attitude: z towards Earth, y along solar panel axis, x roughly towards sun (y cross z). */
  static Rotary3d orientationNominal(const Vector3d &posSat, const Vector3d &posSun);

  /** @brief Computes ECOM DYB satellite orientation.
   *  ECOM attitude: x (d) towards Sun, y along solar panel axis, z (b) roughly towards Earth (d cross y). */
  static Rotary3d orientationEcom(const Vector3d &posSat, const Vector3d &posSun);

  /** @brief Computes nominal satellite orientation for a full arc.
   *  Nominal attitude: z towards Earth, y along solar panel axis, x roughly towards sun (y cross z). */
  static StarCameraArc orientationArcNominal(const OrbitArc &orbitArc, EphemeridesPtr ephemerides);

  /**
   * @brief Models attitude for noon/midnight turns, shadow crossings, and post-shadow maneuvers based on attitude model by Kouba (2009).
   * 'True' attitude: z towards Earth, y along solar panel axis (adjusted for yaw maneuvers), x roughly towards sun (y cross z).
   * @param orbitArc             Orbit arc (must contain position and velocity data).
   * @param svn                  Space vehicle number (e.g. G023).
   * @param block                GNSS satellite block.
   * @param yawBiasArc           Yaw bias data.
   * @param eclipse              eclipse model.
   * @param shadowBuffer         [seconds] pre/post orbit time period to look for missing shadow entries/exits.
   * @param shadowThreshold      Shadow scaling factor for finding exact shadow entry/exit times, between 0 (shadow) and 1 (sun).
   * @return Star camera arc (rotation matrix per epoch).
   */
  StarCameraArc orientationArcKouba2009(const OrbitArc &orbitArc, const std::string &svn, GnssBlock block, const MiscValueArc &yawBiasArc, EphemeridesPtr ephemerides, EclipsePtr eclipse, Double shadowBuffer, Double shadowThreshold) const;

  /**
   * @brief Models attitude for noon/midnight turns based on GLONASS M attitude law derived by Dilssner et al. (2011).
   * 'True' attitude: z towards Earth, y along solar panel axis (adjusted for yaw maneuvers), x roughly towards sun (y cross z).
   * Source: F. Dilssner, T. Springer, G. Gienger, J. Dow. The GLONASS-M satellite yaw-attitude model. Adv. Space Res., 47 (1) (2011), pp. 160-171, 10.1016/j.asr.2010.09.007
   * @param orbitArc             Orbit arc (must contain position and velocity data).
   * @param block                GNSS satellite block.
   * @param eclipse              eclipse model.
   * @param shadowBuffer         [seconds] pre/post orbit time period to look for missing shadow entries/exits.
   * @param shadowThreshold      Shadow scaling factor for finding exact shadow entry/exit times, between 0 (shadow) and 1 (sun).
   * @return Star camera arc (rotation matrix per epoch).
   */
  static StarCameraArc orientationArcGlonassM(const OrbitArc &orbitArc, GnssBlock block, EphemeridesPtr ephemerides, EclipsePtr eclipse, Double shadowBuffer, Double shadowThreshold);

  /**
   * @brief Models attitude for noon/midnight turns based on official Galileo attitude law.
   * 'True' attitude: z towards Earth, y along solar panel axis (adjusted for yaw maneuvers), x roughly towards sun (y cross z).
   * Source: https://www.gsc-europa.eu/support-to-developers/galileo-satellite-metadata
   * @param orbitArc          Orbit arc (must contain position and velocity data).
   * @param integrationLimit  [seconds] pre orbit time period to look for last switch to modified attitude law
   * @return Star camera arc (rotation matrix per epoch).
   */
  static StarCameraArc orientationArcGalileo1(const OrbitArc &orbitArc, EphemeridesPtr ephemerides, Double integrationLimit);
  static StarCameraArc orientationArcGalileo2(const OrbitArc &orbitArc, EphemeridesPtr ephemerides, Double integrationLimit);

  /** @brief Computes nominal yaw angle, nominal yaw rate, sun angle (beta) and orbit angle (mu). */
  static void computeYaw(Epoch &epoch);

  /** @brief Computes yaw angle during noon/midnight turns and post-shadow recovery maneuvers where the satellite can't keep up with the nominal yaw rate. Returns TRUE if it catches up, FALSE otherwise. */
  /**
   * @brief catchUpYaw
   * @param epochStart
   * @param epochs
   * @param yawRateHardwareMax Maximum hardware yaw rate of the satellite.
   * @param allowYawReversal   TRUE to allow yaw reversal (e.g. IIA post-shadow maneuver), FALSE to keep current yaw direction (e.g.
   * @return TRUE if yaw angle has caught up with nominal value, FALSE otherwise.
   */
  static Bool catchUpYaw(const Epoch epochStart, std::vector<Epoch> &epochs, Double yawRateHardwareMax, Bool allowYawReversal);

  /** @brief Polynomial interpolation of velocities from positions. */
  static void interpolateVelocity(OrbitArc &arc, UInt interpolationDegree);

  /**
   * @brief Find exact epoch of shadow entry/exit by Kepler integration.
   * @param epochStart            Start epoch of integration (time, pos, vel).
   * @param integrationStep       Kepler integration step size in [seconds] (can be negative).
   * @param integrationLimit      Kepler integration limit in [seconds] (used as absolute value).
   * @param eclipse               eclipse model.
   * @param findShadowEntry       TRUE = find shadow entry, FALSE = find shadow exit
   * @param shadowFactorThreshold Shadow scaling factor for entry/exit determination (must be between 0 and 1).
   * @param[out] epoch            Epoch of shadow entry/exit.
   * @return TRUE if shadow entry/exit was found, FALSE otherwise.
   */
  static Bool findShadowBoundary(const Epoch &epochStart, Double integrationStep, Double integrationLimit,
                                 EphemeridesPtr ephemerides, EclipsePtr eclipse, Bool findShadowEntry, Double shadowFactorThreshold, const MiscValueArc &yawBiasArc, Epoch &epoch);

  /**
   * @brief Find exact epoch of noon/midnight turn start by Kepler integration.
   * @param epochStart            Start epoch of integration (time, pos, vel).
   * @param integrationStep       Kepler integration step size in [seconds] (can be negative).
   * @param integrationLimit      Kepler integration limit in [seconds] (used as absolute value).
   * @param yawRateHardwareMax    Maximum hardware yaw rate of the satellite.
   * @param yawBiasArc            Yaw bias data.
   * @param[out] epoch            Epoch of noon/midnight turn start.
   * @return TRUE if noon/midnight turn was found, FALSE otherwise.
   */
  static Bool findNoonMidnightTurnStart(EphemeridesPtr ephemerides, const Epoch &epochStart, Double integrationStep, Double integrationLimit, Double yawRateHardwareMax,
                                        const MiscValueArc &yawBiasArc, Epoch &epoch);

  /** @brief Returns yaw bias at @p time. */
  static Double yawBias(const Time &time, const MiscValueArc &yawBiasArc);

  /** @brief Returns @p angle wrapped to [-PI, PI]. */
  static Double wrapAngle(Double angle);

  static Bool blockIsIn(GnssBlock block, const std::vector<GnssBlock> &blockList) { return std::find(blockList.begin(), blockList.end(), block) != blockList.end(); }
};

GROOPS_REGISTER_PROGRAM(SimulateStarCameraGnss, SINGLEPROCESS, "Simulate star camera data for GNSS satellites using nominal orientation or an attitude model.", Simulation, Gnss, Instrument)

/***********************************************/

void SimulateStarCameraGnss::run(Config &config)
{
  try
  {
    FileName       fileNameOutStarCamera, fileNameInOrbit, fileNameTransmitterInfo, fileNameYawBias, fileNameYawRates;
    std::string    attitudeModelChoice;
    Double         shadowBuffer = 3600;
    Double         shadowThreshold = 0.5;
    UInt           interpolationDegree = 8;
    EphemeridesPtr ephemerides;
    EclipsePtr     eclipse;

    readConfig(config,          "outputfileStarCamera", fileNameOutStarCamera,     Config::MUSTSET, "", "rotation from body frame to inertial frame");
    readConfig(config,          "inputfileOrbit",       fileNameInOrbit,           Config::MUSTSET, "", "position defines the orientation of the satellite at each epoch");
    if(readConfigChoice(config, "attitudeModel",        attitudeModelChoice,       Config::MUSTSET, "", "GNSS satellite attitude model"))
    {
      readConfigChoiceElement(config,    "nominal",   attitudeModelChoice, "Nominal attitude: z towards Earth, y along solar panel axis, x roughly towards sun (y cross z)");
      readConfigChoiceElement(config,    "ecom",      attitudeModelChoice, "ECOM attitude: x (d) towards Sun, y along solar panel axis, z (b) roughly towards Earth (d cross y)");
      if(readConfigChoiceElement(config, "Kouba2009", attitudeModelChoice, "True attitude based on model by Kouba (2009): z towards Earth, y along solar panel axis (adjusted for yaw maneuvers), x roughly towards sun (y cross z)"))
      {
        readConfig(config, "inputfileTransmitterInfo", fileNameTransmitterInfo,  Config::MUSTSET,  "{groopsDataDir}/gnss/transmitterGps/transmitterInfo/igs/igs14/transmitterInfo_igs14.{prn}.xml", "GNSS satellite info file");
        readConfig(config, "inputfileYawBias",         fileNameYawBias,          Config::MUSTSET,  "{groopsDataDir}/gnss/transmitterGps/yawBias/yawBias.{prn}.txt", "yaw bias instrument file");
        readConfig(config, "inputfileYawRates",        fileNameYawRates,         Config::DEFAULT,  "{groopsDataDir}/gnss/transmitterGps/yawRates.gps.txt", "yaw rate per SVN in [deg/s]");
        readConfig(config, "shadowBuffer",             shadowBuffer,             Config::DEFAULT,  "3600", "[seconds] pre/post orbit time period to look for missing shadow entries/exits");
        readConfig(config, "shadowThreshold",          shadowThreshold,          Config::DEFAULT,  "0.5",  "shadow scaling factor for finding exact shadow entry/exit times, between 0 (shadow) and 1 (sun)");
        readConfig(config, "interpolationDegree",      interpolationDegree,      Config::DEFAULT,  "8",    "Polynomial degree for velocity interpolation, must be even!");
      }
      if(readConfigChoiceElement(config, "glonass", attitudeModelChoice, "True attitude based on official Galileo attitude law: z towards Earth, y along solar panel axis (adjusted for yaw maneuvers), x roughly towards sun (y cross z)"))
      {
        readConfig(config, "inputfileTransmitterInfo", fileNameTransmitterInfo,  Config::MUSTSET,  "{groopsDataDir}/gnss/transmitterGlonass/transmitterInfo/igs/igs14/transmitterInfo_igs14.{prn}.xml", "GNSS satellite info file");
        readConfig(config, "shadowBuffer",             shadowBuffer,             Config::DEFAULT,  "3600", "[seconds] pre orbit time period to look for last switch to modified attitude law");
        readConfig(config, "shadowThreshold",          shadowThreshold,          Config::DEFAULT,  "0.5",  "shadow scaling factor for finding exact shadow entry/exit times, between 0 (shadow) and 1 (sun)");
        readConfig(config, "interpolationDegree",      interpolationDegree,      Config::DEFAULT,  "8",    "Polynomial degree for velocity interpolation, must be even!");
      }
      if(readConfigChoiceElement(config, "galileo", attitudeModelChoice, "True attitude based on official Galileo attitude law: z towards Earth, y along solar panel axis (adjusted for yaw maneuvers), x roughly towards sun (y cross z)"))
      {
        readConfig(config, "inputfileTransmitterInfo", fileNameTransmitterInfo,  Config::MUSTSET,  "{groopsDataDir}/gnss/transmitterGalileo/transmitterInfo/igs/igs14/transmitterInfo_igs14.{prn}.xml", "GNSS satellite info file");
        readConfig(config, "shadowBuffer",             shadowBuffer,             Config::DEFAULT,  "3600", "[seconds] pre orbit time period to look for last switch to modified attitude law");
        readConfig(config, "interpolationDegree",      interpolationDegree,      Config::DEFAULT,  "8",    "Polynomial degree for velocity interpolation, must be even!");
      }
      endChoice(config);
    }
    readConfig(config, "ephemerides", ephemerides, Config::MUSTSET, "", "");
    readConfig(config, "eclipse",     eclipse,     Config::MUSTSET, "", "model to determine if satellite is in Earth's shadow");
    if(isCreateSchema(config)) return;

    std::list<Arc> arcListOut;

    // create star camera data for nominal or ECOM attitude
    // ----------------------------------------------------
    if(attitudeModelChoice == "nominal" || attitudeModelChoice == "ecom")
    {
      logStatus << "read orbit and generate star camera data" << Log::endl;
      InstrumentFile  orbitFile(fileNameInOrbit);
      UInt arcCount = orbitFile.arcCount();
      UInt epochCount = 0;

      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      {
        OrbitArc orbit = orbitFile.readArc(arcNo);
        UInt posCount  = orbit.size();

        StarCameraArc arc;
        for(UInt i=0; i<posCount; i++)
        {
          StarCameraEpoch epoch;
          epoch.time   = orbit.at(i).time;
          if(attitudeModelChoice == "ecom")
            epoch.rotary = orientationEcom(orbit.at(i).position, ephemerides->position(orbit.at(i).time, Ephemerides::SUN));
          else
            epoch.rotary = orientationNominal(orbit.at(i).position, ephemerides->position(orbit.at(i).time, Ephemerides::SUN));
          arc.push_back(epoch);
          epochCount++;
        }
        arcListOut.push_back(arc);
      }

      if(epochCount == 0)
        throw(Exception("empty input orbit"));
    }

    // create star camera data for attitude model by Kouba (2009)
    // ----------------------------------------------------------
    else if(attitudeModelChoice == "Kouba2009" || attitudeModelChoice == "glonass" || attitudeModelChoice == "galileo")
    {
      if(shadowThreshold < 0 || shadowThreshold > 1)
        throw(Exception("shadow threshold must be between 0 (shadow) and 1 (sun): " + shadowThreshold%"%f"s));

      if(interpolationDegree%2 != 0)
        throw(Exception("interpolation degree must be even"));

      // read orbit arc by arc to remember arc structure for output starCamera file
      logStatus << "read orbit file <" << fileNameInOrbit << ">" << Log::endl;
      InstrumentFile  orbitFile(fileNameInOrbit);
      OrbitArc orbitArc;
      std::vector<std::vector<Time>> timesOrbitArcOriginal;
      for(UInt arcNo = 0; arcNo < orbitFile.arcCount(); arcNo++)
      {
        OrbitArc arc = orbitFile.readArc(arcNo);
        orbitArc.append(arc);
        timesOrbitArcOriginal.push_back(arc.times());
      }
      if(!orbitArc.size())
        throw(Exception("empty input orbit"));

      // sort and remove duplicate epochs
      orbitArc.sort();
      for(UInt i = 1; i < orbitArc.size(); i++)
        if(orbitArc.at(i).time == orbitArc.at(i-1).time)
          orbitArc.remove(i--);

      // interpolate velocities if necessary
      if(orbitArc.at(0).velocity.r() == 0 || orbitArc.at(orbitArc.size()-1).velocity.r() == 0)
      {
        logStatus << "velocity data missing, interpolating from orbit" << Log::endl;
        interpolateVelocity(orbitArc, interpolationDegree);
      }

      logStatus << "read transmitter info file <" << fileNameTransmitterInfo << ">" << Log::endl;
      GnssStationInfo transmitterInfo;
      readFileGnssStationInfo(fileNameTransmitterInfo, transmitterInfo);

      // find SVN and GNSS block of satellite
      std::string svn;
      GnssBlock   block = UNKNOWN;
      UInt idAnt = transmitterInfo.findAntenna(orbitArc.times().at(0));
      if(idAnt != NULLINDEX)
      {
        svn   = transmitterInfo.antenna.at(idAnt).serial;
        block = string2block(transmitterInfo.antenna.at(idAnt).name);
      }
      if(svn.empty() || block == UNKNOWN)
        throw(Exception("satellite not found in transmitter info file"));

      StarCameraArc starCameraArc;
      if(attitudeModelChoice == "Kouba2009")
      {
        logStatus << "read yaw bias file <" << fileNameYawBias << ">" << Log::endl;
        MiscValueArc yawBiasArc = InstrumentFile::read(fileNameYawBias);
        if(!yawBiasArc.size())
          throw(Exception("yaw bias file is empty"));

        if(!fileNameYawRates.empty())
        {
          logStatus << "read yaw rates file <" << fileNameYawRates << ">" << Log::endl;
          readFileMatrix(fileNameYawRates, yawRates);
        }

        starCameraArc = orientationArcKouba2009(orbitArc, svn, block, yawBiasArc, ephemerides, eclipse, shadowBuffer, shadowThreshold);
      }
      else if(attitudeModelChoice == "glonass")
      {
        starCameraArc = orientationArcGlonassM(orbitArc, block, ephemerides, eclipse, shadowBuffer, shadowThreshold);
      }
      else if(attitudeModelChoice == "galileo")
      {
        if(block == GALILEO_1)
          starCameraArc = orientationArcGalileo1(orbitArc, ephemerides, shadowBuffer);
        else if(block == GALILEO_2)
          starCameraArc = orientationArcGalileo2(orbitArc, ephemerides, shadowBuffer);
        else
        {
          logWarning << "no attitude model implemented for "+svn+" ("+block2string(block)+"), using nominal attitude" << Log::endl;
          starCameraArc = orientationArcNominal(orbitArc, ephemerides);
        }
      }

      // recover input orbit arc structure for output starCamera file
      for(UInt i = 0; i < timesOrbitArcOriginal.size(); i++)
      {
        StarCameraArc arc(starCameraArc);
        arc.synchronize(timesOrbitArcOriginal.at(i));
        arcListOut.push_back(arc);
      }
    }

    logStatus << "write star camera file <" << fileNameOutStarCamera << ">" << Log::endl;
    InstrumentFile::write(fileNameOutStarCamera, arcListOut);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d SimulateStarCameraGnss::orientationNominal(const Vector3d &posSat, const Vector3d &posSun)
{
  try
  {
    Vector3d z = normalize(-posSat);
    Vector3d x = normalize(posSun - posSat);
    Vector3d y = normalize(crossProduct(z,x));
    return Rotary3d(crossProduct(y,z), y);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d SimulateStarCameraGnss::orientationEcom(const Vector3d &posSat, const Vector3d &posSun)
{
  try
  {
    Vector3d d = normalize(posSun - posSat);
    return Rotary3d(d, crossProduct(normalize(-posSat), d));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

StarCameraArc SimulateStarCameraGnss::orientationArcNominal(const OrbitArc &orbitArc, EphemeridesPtr ephemerides)
{
  try
  {
    StarCameraArc arc;
    for(UInt i=0; i<orbitArc.size(); i++)
    {
      StarCameraEpoch epoch;
      epoch.time   = orbitArc.at(i).time;
      epoch.rotary = orientationNominal(orbitArc.at(i).position, ephemerides->position(orbitArc.at(i).time, Ephemerides::SUN));
      arc.push_back(epoch);
    }
    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

StarCameraArc SimulateStarCameraGnss::orientationArcGlonassM(const OrbitArc &orbitArc, GnssBlock block, EphemeridesPtr ephemerides, EclipsePtr eclipse, Double shadowBuffer, Double shadowThreshold)
{
  try
  {
    // Nominal maximum hardware yaw rate
    Double yawRateHardwareMax = 0;
    if(block == GLONASS_M)
      yawRateHardwareMax = 0.25*DEG2RAD;

    std::vector<Time> times = orbitArc.times();
    std::vector<Epoch> epoch(times.size());

    // Nominal yaw and yaw rate
    for(UInt k = 0; k < times.size(); k++)
    {
      epoch.at(k).time    = times.at(k);
      epoch.at(k).pos     = orbitArc.at(k).position;
      epoch.at(k).vel     = orbitArc.at(k).velocity;
      epoch.at(k).posSun  = ephemerides->position(times.at(k), Ephemerides::SUN);
      computeYaw(epoch.at(k));
    }

    // =============================================

    // Attitude correction for noon/midnight turns
    // -------------------------------------------

    // Find noon/midnight turn starts
//    std::vector<Epoch> turnStart;
//    Epoch epochTurn;
//    for(UInt k = 0; k < times.size()-1; k++)
//      if(std::fabs(epoch.at(k).yawRateTrue) < yawRateHardwareMax && std::fabs(epoch.at(k+1).yawRateTrue) >= yawRateHardwareMax)
//        if(findNoonMidnightTurnStart(ephemerides, epoch.at(k), 1, (epoch.at(k+1).time-epoch.at(k).time).seconds(), yawRateHardwareMax, yawBiasArc, epochTurn))
//          turnStart.push_back(epochTurn);

//    // Find missing noon/midnight turn start if satellite is already potentially in noon/midnight turn when interval starts
//    if(turnStart.size() && findNoonMidnightTurnStart(ephemerides, epoch.at(0), -1, shadowBuffer, yawRateHardwareMax, yawBiasArc, epochTurn))
//        turnStart.insert(turnStart.begin(), epochTurn);

//    // Compute true yaw angle during noon/midnight turns
//    for(UInt i = 0; i < turnStart.size(); i++)
//      catchUpYaw(turnStart.at(i), epoch, yawRateHardwareMax, FALSE);

    // =============================================

    // Attitude correction for Block II/IIA and Block IIF shadow crossings
    // -------------------------------------------------------------------

    if(block == GLONASS_M)
    {
      std::vector<Epoch> shadowStart, shadowEnd;

      for(UInt k = 0; k < times.size(); k++)
        epoch.at(k).shadowFactor = eclipse ? eclipse->factor(epoch.at(k).time, epoch.at(k).pos, ephemerides) : 1.;

      // Find start and end epochs of shadow crossings
      Epoch epochShadow;
      for(UInt k = 0; k < times.size()-1; k++)
      {
        // check for shadow entry
        if(epoch.at(k).shadowFactor >= shadowThreshold && epoch.at(k+1).shadowFactor < shadowThreshold)
          if(findShadowBoundary(epoch.at(k), 1, (epoch.at(k+1).time-epoch.at(k).time).seconds(), ephemerides, eclipse, TRUE, shadowThreshold, MiscValueArc(), epochShadow))
            shadowStart.push_back(epochShadow);

        // check for shadow exit
        if(epoch.at(k).shadowFactor < shadowThreshold && epoch.at(k+1).shadowFactor >= shadowThreshold)
          if(findShadowBoundary(epoch.at(k), 1, (epoch.at(k+1).time-epoch.at(k).time).seconds(), ephemerides, eclipse, FALSE, shadowThreshold, MiscValueArc(), epochShadow))
            shadowEnd.push_back(epochShadow);
      }

      // Find missing shadow exit if satellite is still in shadow when interval ends
      if(shadowStart.size() && (!shadowEnd.size() || (shadowEnd.size() && shadowEnd.back().time < shadowStart.back().time)))
        if(findShadowBoundary(epoch.back(), 1, shadowBuffer, ephemerides, eclipse, FALSE, shadowThreshold, MiscValueArc(), epochShadow))
          shadowEnd.push_back(epochShadow);

      // Find missing shadow entry if satellite is already in shadow when interval starts
      if(shadowEnd.size() && (!shadowStart.size() || (shadowStart.size() && shadowEnd.at(0).time < shadowStart.at(0).time)))
        if(findShadowBoundary(epoch.at(0), -1, shadowBuffer, ephemerides, eclipse, TRUE, shadowThreshold, MiscValueArc(), epochShadow))
          shadowStart.insert(shadowStart.begin(), epochShadow);

      if(shadowStart.size() != shadowEnd.size())
        throw(Exception("number of shadow entries and exits not equal: " + shadowStart.size()%"%i"s + " != " + shadowEnd.size()%"%i"s));

      // Compute true yaw angle during shadow crossings
      for(UInt i = 0; i < shadowStart.size(); i++)
      {
        // Yaw rate at shadow start and during shadow crossing
        shadowStart.at(i).yawRateTrue = (shadowStart.at(i).yawRateNominal >= 0 ? yawRateHardwareMax : -yawRateHardwareMax);

        // Yaw angle during shadow crossing
        for(UInt k = 0; k < times.size(); k++)
          if(times.at(k) > shadowStart.at(i).time && times.at(k) <= shadowEnd.at(i).time)
          {
//            epoch.at(k).yawTrue     = wrapAngle(shadowStart.at(i).yawTrue + shadowStart.at(i).yawRateTrue*(times.at(k)-shadowStart.at(i).time).seconds());
//            epoch.at(k).yawRateTrue = shadowStart.at(i).yawRateTrue;
            epoch.at(k).yawTrue     = shadowEnd.at(i).yawNominal;
            epoch.at(k).yawRateTrue = 0;
          }
        catchUpYaw(shadowStart.at(i), epoch, yawRateHardwareMax, FALSE);
      }
    } // end if(block == GLONASS_M)
    else
      logWarning << "no attitude model implemented for block "+block2string(block)+", using nominal attitude" << Log::endl;

    // =============================================

    // Compute corrected orientation
    StarCameraArc starCameraArc;
    for(UInt k = 0; k < times.size(); k++)
    {
      StarCameraEpoch starCameraEpoch;
      starCameraEpoch.time   = times.at(k);
      starCameraEpoch.rotary = orientationNominal(epoch.at(k).pos, epoch.at(k).posSun) * rotaryZ(-Angle(epoch.at(k).yawTrue - epoch.at(k).yawNominal));
      starCameraArc.push_back(starCameraEpoch);
    }

    return starCameraArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

StarCameraArc SimulateStarCameraGnss::orientationArcGalileo1(const OrbitArc &orbitArc, EphemeridesPtr ephemerides, Double integrationLimit)
{
  try
  {
    auto computePsi = [](const Vector3d &S)
    {
      const Double denom = std::sqrt(1. - S.z()*S.z());
      return std::atan2(S.y()/denom, S.x()/denom);
    };

    auto computeValues = [&](const Time &time, const Vector3d &pos, const Vector3d &vel, Vector3d &S0, Rotary3d &orbitalReferenceSystem)
    {
      const Vector3d r = normalize(pos);
      const Vector3d v = normalize(vel);
      const Vector3d n = crossProduct(r, v);  // orbit normal vector

      orbitalReferenceSystem = Rotary3d(crossProduct(-r, n), -n);
      S0 = orbitalReferenceSystem.inverseRotate(normalize(ephemerides->position(time, Ephemerides::SUN) - pos));
    };

    StarCameraArc starCameraArc;

    Vector3d S0LastSwitch, S0LastEpoch;
    for(UInt idEpoch = 0; idEpoch < orbitArc.size(); idEpoch++)
    {
      Vector3d S0;
      Rotary3d orbitalReferenceSystem;
      computeValues(orbitArc.at(idEpoch).time, orbitArc.at(idEpoch).position, orbitArc.at(idEpoch).velocity, S0, orbitalReferenceSystem);
      Double psi = computePsi(S0);

      const Double sinBetaX = std::sin(15.0*DEG2RAD);
      const Double sinBetaY = std::sin(2.0*DEG2RAD);

      if(std::fabs(S0.x()) < sinBetaX && std::fabs(S0.y()) < sinBetaY)
      {
        // find S0 for last switch to modified yaw steering law
        if(idEpoch == 0) // already using modified yaw steering law at initial epoch => find last switch by extrapolating backwards using Kepler orbit
        {
          Kepler kepler(orbitArc.at(idEpoch).time, orbitArc.at(idEpoch).position, orbitArc.at(idEpoch).velocity);

          const Double integrationStep = -1;
          for(UInt i = 0; std::fabs(i*integrationStep) <= std::fabs(integrationLimit)+1; i++)
          {
            const Time time = orbitArc.at(idEpoch).time + seconds2time(i*integrationStep);
            Vector3d pos, vel, S0;
            Rotary3d orbitalReferenceSystem;
            kepler.orbit(time, pos, vel);
            computeValues(time, pos, vel, S0, orbitalReferenceSystem);

            if(std::fabs(S0.x()) > sinBetaX || std::fabs(S0.y()) > sinBetaY)
            {
              S0LastSwitch = S0;
              break;
            }
          }
        }
        else if(std::fabs(S0LastEpoch.x()) > sinBetaX || std::fabs(S0LastEpoch.y()) > sinBetaY)
          S0LastSwitch = S0;

        const Double Gamma = (S0LastSwitch.y() < 0 ? -1 : 1);
        const Double SHy = 0.5*(sinBetaY*Gamma + S0.y()) + 0.5*(sinBetaY*Gamma - S0.y()) * std::cos(PI*std::fabs(S0.x())/sinBetaX);
        const Vector3d SH(S0.x(), SHy, std::sqrt(1. - S0.x()*S0.x() - SHy*SHy)*(S0.z() < 0 ? -1 : 1));

        psi = computePsi(SH);
      }

      S0LastEpoch = S0;

      StarCameraEpoch starCameraEpoch;
      starCameraEpoch.time   = orbitArc.at(idEpoch).time;
      starCameraEpoch.rotary = orbitalReferenceSystem * rotaryZ(-Angle(psi));
      starCameraArc.push_back(starCameraEpoch);
    }

    return starCameraArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

StarCameraArc SimulateStarCameraGnss::orientationArcGalileo2(const OrbitArc &orbitArc, EphemeridesPtr ephemerides, Double integrationLimit)
{
  try
  {
    auto computeValues = [ephemerides](const Time &time, const Vector3d &pos, const Vector3d &vel, Double &beta, Double &epsilon, Double &psi, Rotary3d &orbitalReferenceSystem)
    {
      const Vector3d r = normalize(pos);
      const Vector3d v = normalize(vel);
      const Vector3d n = crossProduct(r, v);  // orbit normal vector
      const Vector3d s = normalize(ephemerides->position(time, Ephemerides::SUN));

      beta = std::acos(inner(normalize(crossProduct(v, r)), s)) - PI/2;
      epsilon = std::acos(inner(r, crossProduct(n, crossProduct(n, s))));
      if(epsilon > PI/2)
        epsilon = PI - epsilon;

      psi = std::atan2(inner(-s, n), inner(-s, crossProduct(r, n)));
      orbitalReferenceSystem = Rotary3d(crossProduct(-r, n), -n);
    };

    StarCameraArc starCameraArc;

    Double psiLastSwitch = 0;
    Time   timeLastSwitch;
    Double epsilonLastEpoch = 0;
    for(UInt idEpoch = 0; idEpoch < orbitArc.size(); idEpoch++)
    {
      Double beta, epsilon, psi;
      Rotary3d orbitalReferenceSystem;
      computeValues(orbitArc.at(idEpoch).time, orbitArc.at(idEpoch).position, orbitArc.at(idEpoch).velocity, beta, epsilon, psi, orbitalReferenceSystem);

      if(std::fabs(beta) < 4.1*DEG2RAD && std::fabs(epsilon) < 10.0*DEG2RAD)
      {
        // find time and psi at last switch to modified yaw steering law
        if(idEpoch == 0) // already using modified yaw steering law at initial epoch => find last switch by extrapolating backwards using Kepler orbit
        {
          Kepler kepler(orbitArc.at(idEpoch).time, orbitArc.at(idEpoch).position, orbitArc.at(idEpoch).velocity);

          const Double integrationStep = -1;
          for(UInt i = 0; std::fabs(i*integrationStep) <= std::fabs(integrationLimit)+1; i++)
          {
            const Time time = orbitArc.at(idEpoch).time + seconds2time(i*integrationStep);
            Vector3d pos, vel;
            Double   beta, epsilon, psi;
            Rotary3d orbitalReferenceSystem;
            kepler.orbit(time, pos, vel);
            computeValues(time, pos, vel, beta, epsilon, psi, orbitalReferenceSystem);

            if(std::fabs(epsilon) >= 10.0*DEG2RAD)
            {
              psiLastSwitch  = psi;
              timeLastSwitch = time;
              break;
            }
          }
        }
        else if(std::fabs(epsilonLastEpoch) >= 10.0*DEG2RAD)
        {
          psiLastSwitch = psi;
          timeLastSwitch = orbitArc.at(idEpoch).time;
        }

        const Double offset = (psiLastSwitch < 0 ? -PI/2 : PI/2);
        const Double tMod   = (orbitArc.at(idEpoch).time - timeLastSwitch).seconds();
        psi = offset + (psiLastSwitch - offset) * std::cos(2*PI/5656 * tMod);
      }

      epsilonLastEpoch = epsilon;

      StarCameraEpoch starCameraEpoch;
      starCameraEpoch.time   = orbitArc.at(idEpoch).time;
      starCameraEpoch.rotary = orbitalReferenceSystem * rotaryZ(-Angle(psi));
      starCameraArc.push_back(starCameraEpoch);
    }

    return starCameraArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Based on [Kouba (2009), A simplified yaw-attitude model for eclipsing GPS satellites] and IGS Repro 2 updates
StarCameraArc SimulateStarCameraGnss::orientationArcKouba2009(const OrbitArc &orbitArc, const std::string &svn, GnssBlock block,
                                                           const MiscValueArc &yawBiasArc, EphemeridesPtr ephemerides, EclipsePtr eclipse, Double shadowBuffer, Double shadowThreshold) const
{
  try
  {
    std::vector<Time> times = orbitArc.times();
    std::vector<Epoch> epoch(times.size());

    // Nominal yaw and yaw rate
    for(UInt k = 0; k < times.size(); k++)
    {
      epoch.at(k).time    = times.at(k);
      epoch.at(k).pos     = orbitArc.at(k).position;
      epoch.at(k).vel     = orbitArc.at(k).velocity;
      epoch.at(k).posSun  = ephemerides->position(times.at(k), Ephemerides::SUN);
      epoch.at(k).yawBias = yawBias(times.at(k), yawBiasArc);
      computeYaw(epoch.at(k));
    }

    // =============================================

    if(blockIsIn(block, {GPS_II, GPS_IIA, GPS_IIR, GPS_IIRA, GPS_IIRB, GPS_IIRM, GPS_IIF}))
    {
      // Nominal maximum hardware yaw rate
      Double yawRateHardwareMax = 0;
      if(blockIsIn(block, {GPS_II, GPS_IIA}))
        yawRateHardwareMax = 0.12*DEG2RAD;
      if(blockIsIn(block, {GPS_IIR, GPS_IIRA, GPS_IIRB, GPS_IIRM}))
        yawRateHardwareMax = 0.20*DEG2RAD;
      if(block == GPS_IIF)
        yawRateHardwareMax = 0.11*DEG2RAD;

      // Search for estimated maximum hardware yaw rate
      if(yawRates.size())
      {
        UInt svnInt = std::atoi(svn.substr(1).c_str());
        for(UInt i = 0; i < yawRates.rows(); i++)
          if(yawRates(i,0) == svnInt)
          {
            yawRateHardwareMax = yawRates(i,1)*DEG2RAD;
            break;
          }
      }

      // =============================================

      // Attitude correction for noon/midnight turns
      // -------------------------------------------

      // Find noon/midnight turn starts
      std::vector<Epoch> turnStart;
      Epoch epochTurn;
      for(UInt k = 0; k < times.size()-1; k++)
        if(std::fabs(epoch.at(k).yawRateTrue) < yawRateHardwareMax && std::fabs(epoch.at(k+1).yawRateTrue) >= yawRateHardwareMax)
          if(findNoonMidnightTurnStart(ephemerides, epoch.at(k), 1, (epoch.at(k+1).time-epoch.at(k).time).seconds(), yawRateHardwareMax, yawBiasArc, epochTurn))
            turnStart.push_back(epochTurn);

      // Find missing noon/midnight turn start if satellite is already potentially in noon/midnight turn when interval starts
      if(turnStart.size() && findNoonMidnightTurnStart(ephemerides, epoch.at(0), -1, shadowBuffer, yawRateHardwareMax, yawBiasArc, epochTurn))
          turnStart.insert(turnStart.begin(), epochTurn);

      // Compute true yaw angle during noon/midnight turns
      for(UInt i = 0; i < turnStart.size(); i++)
        catchUpYaw(turnStart.at(i), epoch, yawRateHardwareMax, FALSE);

      // =============================================

      // Attitude correction for Block II/IIA and Block IIF shadow crossings
      // -------------------------------------------------------------------

      if(blockIsIn(block, {GPS_II, GPS_IIA, GPS_IIF}))
      {
        std::vector<Epoch> shadowStart, shadowEnd;

        for(UInt k = 0; k < times.size(); k++)
          epoch.at(k).shadowFactor = eclipse ? eclipse->factor(epoch.at(k).time, epoch.at(k).pos, ephemerides) : 1.;

        // Find start and end epochs of shadow crossings
        Epoch epochShadow;
        for(UInt k = 0; k < times.size()-1; k++)
        {
          // check for shadow entry
          if(epoch.at(k).shadowFactor >= shadowThreshold && epoch.at(k+1).shadowFactor < shadowThreshold)
            if(findShadowBoundary(epoch.at(k), 1, (epoch.at(k+1).time-epoch.at(k).time).seconds(), ephemerides, eclipse, TRUE, shadowThreshold, yawBiasArc, epochShadow))
              shadowStart.push_back(epochShadow);

          // check for shadow exit
          if(epoch.at(k).shadowFactor < shadowThreshold && epoch.at(k+1).shadowFactor >= shadowThreshold)
            if(findShadowBoundary(epoch.at(k), 1, (epoch.at(k+1).time-epoch.at(k).time).seconds(), ephemerides, eclipse, FALSE, shadowThreshold, yawBiasArc, epochShadow))
              shadowEnd.push_back(epochShadow);
        }

        // Find missing shadow exit if satellite is still in shadow when interval ends
        if(shadowStart.size() && (!shadowEnd.size() || (shadowEnd.size() && shadowEnd.back().time < shadowStart.back().time)))
          if(findShadowBoundary(epoch.back(), 1, shadowBuffer, ephemerides, eclipse, FALSE, shadowThreshold, yawBiasArc, epochShadow))
            shadowEnd.push_back(epochShadow);

        // Find missing shadow entry if satellite is already in shadow when interval starts
        if(shadowEnd.size() && (!shadowStart.size() || (shadowStart.size() && shadowEnd.at(0).time < shadowStart.at(0).time)))
          if(findShadowBoundary(epoch.at(0), -1, shadowBuffer, ephemerides, eclipse, TRUE, shadowThreshold, yawBiasArc, epochShadow))
            shadowStart.insert(shadowStart.begin(), epochShadow);

        if(shadowStart.size() != shadowEnd.size())
          throw(Exception("number of shadow entries and exits not equal: " + shadowStart.size()%"%i"s + " != " + shadowEnd.size()%"%i"s));

        // Compute true yaw angle during shadow crossings
        for(UInt i = 0; i < shadowStart.size(); i++)
        {
          // Yaw rate at shadow start and during shadow crossing
          if(blockIsIn(block, {GPS_II, GPS_IIA}))
            shadowStart.at(i).yawRateTrue = (yawBias(shadowStart.at(i).time, yawBiasArc) >= 0 ? yawRateHardwareMax : -yawRateHardwareMax);
          if(block == GPS_IIF)
            shadowStart.at(i).yawRateTrue = wrapAngle(shadowEnd.at(i).yawTrue-shadowStart.at(i).yawTrue) / (shadowEnd.at(i).time-shadowStart.at(i).time).seconds();

          // Yaw angle during shadow crossing
          for(UInt k = 0; k < times.size(); k++)
            if(times.at(k) > shadowStart.at(i).time && times.at(k) <= shadowEnd.at(i).time)
            {
              epoch.at(k).yawTrue     = wrapAngle(shadowStart.at(i).yawTrue + shadowStart.at(i).yawRateTrue*(times.at(k)-shadowStart.at(i).time).seconds());
              epoch.at(k).yawRateTrue = shadowStart.at(i).yawRateTrue;
            }

          // Yaw angle at shadow exit
          shadowEnd.at(i).yawTrue = wrapAngle(shadowStart.at(i).yawTrue + shadowStart.at(i).yawRateTrue*(shadowEnd.at(i).time-shadowStart.at(i).time).seconds());

          // Block II/IIA post-shadow recovery maneuver
          if(blockIsIn(block, {GPS_II, GPS_IIA}))
            catchUpYaw(shadowEnd.at(i), epoch, yawRateHardwareMax, TRUE);
        }
      } // end if(block == GPS_II || block == GPS_IIA || block == GPS_IIF)
    }
    else
      logWarning << "no attitude model implemented for "+svn+" ("+block2string(block)+"), using nominal attitude" << Log::endl;

    // =============================================

    // Compute corrected orientation
    StarCameraArc starCameraArc;
    for(UInt k = 0; k < times.size(); k++)
    {
      StarCameraEpoch starCameraEpoch;
      starCameraEpoch.time   = times.at(k);
      starCameraEpoch.rotary = orientationNominal(epoch.at(k).pos, epoch.at(k).posSun) * rotaryZ(-Angle(epoch.at(k).yawTrue - epoch.at(k).yawNominal));
      starCameraArc.push_back(starCameraEpoch);
    }

    return starCameraArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Based on [Kouba (2009), A simplified yaw-attitude model for eclipsing GPS satellites]
void SimulateStarCameraGnss::computeYaw(Epoch &epoch)
{
  try
  {
    // Orbit angle rate
    Double muRate = epoch.vel.r()/epoch.pos.r();

    // Argument of latitude of satellite and sun
    Vector3d z  = normalize(crossProduct(epoch.pos, epoch.vel));
    Vector3d x  = normalize(crossProduct(Vector3d(0,0,1), z));
    Vector3d y  = crossProduct(z, x);
    Double   u  = atan2(inner(epoch.pos,y), inner(epoch.pos,x));
    Double   u0 = atan2(inner(epoch.posSun,y), inner(epoch.posSun,x));

    // Orbit angle (mu, [0..360°], 0° = midnight) and sun angle (beta, [-90..90°], 0° = sun in sat orbit plane)
    epoch.mu   = std::fmod(2*PI + (u - u0 + PI), 2*PI);
    epoch.beta = std::acos(inner(-z, epoch.posSun)/epoch.posSun.r()) - PI/2;

    // Noon turn reversal for small positive beta angles (if yaw bias > 0, II/IIA) or small negative beta angles (if yaw bias < 0, IIF)
    Double betaBias = 0;
    if(std::fabs(epoch.beta) < std::fabs(epoch.yawBias)*DEG2RAD && epoch.yawBias*epoch.beta > 0)
      betaBias = -epoch.yawBias*DEG2RAD;

    // Lambda function for yaw and yaw rate computation
    auto compute = [&] (Double beta, Double &yaw, Double &yawRate)
    {
      // Nominal yaw angle and yaw rate
      yaw     = std::atan2(-std::tan(beta), std::sin(epoch.mu));
      yawRate = muRate * std::tan(beta)*std::cos(epoch.mu) / (std::sin(epoch.mu)*std::sin(epoch.mu) + std::tan(beta)*std::tan(beta));
    };

    compute(epoch.beta,          epoch.yawNominal, epoch.yawRateNominal);
    compute(epoch.beta+betaBias, epoch.yawTrue,    epoch.yawRateTrue);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool SimulateStarCameraGnss::catchUpYaw(const Epoch epochStart, std::vector<Epoch> &epochs, Double yawRateHardwareMax, Bool allowYawReversal)
{
  try
  {
    std::vector<Double> yawDiff(epochs.size());
    for(UInt k = 1; k < epochs.size(); k++)
    {
      if(epochs.at(k).time <= epochStart.time)
        continue;

      Double yawReversalFactor = (yawDiff.at(k-1) < 0 ? -1. : 1.);
      Double yawRate = allowYawReversal ? yawReversalFactor*yawRateHardwareMax : epochStart.yawRateTrue;
      Double yawNew = wrapAngle(epochStart.yawTrue + yawRate*(epochs.at(k).time-epochStart.time).seconds());
      Double yawDiffNew = wrapAngle(epochs.at(k).yawTrue - yawNew);

      // Stop if true yaw angle catches up with nominal yaw angle (or is numerically identical)
      if((yawDiff.at(k-1) > 0 && yawDiffNew < 0) || (yawDiff.at(k-1) < 0 && yawDiffNew > 0) || std::fabs(yawDiffNew) < 1e-5*DEG2RAD)
        return TRUE;

      yawDiff.at(k)            = yawDiffNew;
      epochs.at(k).yawTrue     = yawNew;
      epochs.at(k).yawRateTrue = yawRate;
    }

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SimulateStarCameraGnss::interpolateVelocity(OrbitArc &arc, UInt interpolationDegree)
{
  try
  {
    UInt idx = 0;
    for(UInt idEpoch=0; idEpoch<arc.size(); idEpoch++)
    {
      // find optimal interval
      // ---------------------
      while((idx+interpolationDegree < arc.size()) && (arc.at(idx+interpolationDegree).time < arc.at(idEpoch).time))
        idx++;
      if(idx+interpolationDegree >= arc.size())
        break;

      UInt   idxOpt   = MAX_UINT;
      Double deltaOpt = 1e99;
      while((idx+interpolationDegree < arc.size()) && (arc.at(idx).time <= arc.at(idEpoch).time))
      {
        // interpolation point should be in the mid of the interval
        // => search minimum of the difference of the time before and after the interpolation point
        const Double delta = std::fabs(((arc.at(idx+interpolationDegree).time-arc.at(idEpoch).time)-(arc.at(idEpoch).time-arc.at(idx).time)).seconds());
        if(delta <= deltaOpt)
        {
          idxOpt   = idx;
          deltaOpt = delta;
        }
        idx++;
      }
      idx = idxOpt;

      // polynomial interpolation
      // ------------------------
      Matrix A(interpolationDegree+1, interpolationDegree+1);
      for(UInt k=0; k<interpolationDegree+1; k++)
      {
        const Double factor = (arc.at(idx+k).time-arc.at(idEpoch).time).seconds();
        A(0,k) = 1.0;
        for(UInt n=1; n<=interpolationDegree; n++)
          A(n,k) = factor * A(n-1,k);
      }
      Matrix coeff(interpolationDegree+1, 1);
      coeff(1, 0) = 1.; // velocity
      solveInPlace(A, coeff);

      arc.at(idEpoch).velocity = Vector3d();
      for(UInt k=0; k<coeff.rows(); k++)
        arc.at(idEpoch).velocity += coeff(k,0) * arc.at(idx+k).position;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool SimulateStarCameraGnss::findShadowBoundary(const Epoch &epochStart, Double integrationStep, Double integrationLimit,
                                                EphemeridesPtr ephemerides, EclipsePtr eclipse, Bool findShadowEntry, Double shadowFactorThreshold, const MiscValueArc &yawBiasArc, Epoch &epoch)
{
  try
  {
    if(shadowFactorThreshold < 0 || shadowFactorThreshold > 1)
      throw(Exception("shadow factor threshold must be between 0 (shadow) and 1 (sun): " + shadowFactorThreshold%"%f"s));

    if(integrationStep == 0)
      throw(Exception("integration step can't be zero"));

    // shadow threshold crossing check (see below) has to be reversed for backwards integration
    if(integrationStep < 0)
      findShadowEntry = !findShadowEntry;

    epoch = Epoch();
    Kepler kepler(epochStart.time, epochStart.pos, epochStart.vel);

    for(UInt i = 0; std::fabs(i*integrationStep) <= std::fabs(integrationLimit)+1; i++)
    {
      epoch.time         = epochStart.time + seconds2time(i*integrationStep);
      kepler.orbit(epoch.time, epoch.pos, epoch.vel);
      epoch.posSun       = ephemerides->position(epoch.time, Ephemerides::SUN);
      epoch.shadowFactor = eclipse ? eclipse->factor(epoch.time, epoch.pos, ephemerides) : 1.;

      if(( findShadowEntry && epoch.shadowFactor <= shadowFactorThreshold) ||
         (!findShadowEntry && epoch.shadowFactor >= shadowFactorThreshold))
      {
        epoch.yawBias = yawBiasArc.size() ? yawBias(epoch.time, yawBiasArc) : 0;
        computeYaw(epoch);
        return TRUE;
      }
    }

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool SimulateStarCameraGnss::findNoonMidnightTurnStart(EphemeridesPtr ephemerides, const Epoch &epochStart, Double integrationStep, Double integrationLimit,
                                                       Double yawRateHardwareMax, const MiscValueArc &yawBiasArc, Epoch &epoch)
{
  try
  {
    if(integrationStep == 0)
      throw(Exception("integration step cannot be zero"));

    Double yawRatePrevious = 0;
    Kepler kepler(epochStart.time, epochStart.pos, epochStart.vel);
    std::vector<Epoch> epochs;

    for(UInt i = 0; std::fabs(i*integrationStep) <= std::fabs(integrationLimit)+1; i++)
    {
      epoch = Epoch();
      epoch.time   = epochStart.time + seconds2time(i*integrationStep);
      epoch.posSun = ephemerides->position(epoch.time, Ephemerides::SUN);
      kepler.orbit(epoch.time, epoch.pos, epoch.vel);
      epoch.yawBias = yawBias(epoch.time, yawBiasArc);
      computeYaw(epoch);
      epochs.push_back(epoch);

      if((integrationStep > 0 && std::fabs(epoch.yawRateTrue) >= yawRateHardwareMax && std::fabs(yawRatePrevious) <  yawRateHardwareMax) || // forward integration
         (integrationStep < 0 && std::fabs(epoch.yawRateTrue) <  yawRateHardwareMax && std::fabs(yawRatePrevious) >= yawRateHardwareMax))   // backward integration
      {
        epoch.yawRateTrue = (epoch.yawRateTrue >= 0 ? yawRateHardwareMax : -yawRateHardwareMax);

        // ignore turn start if the turn maneuver ends before epochStart
        if(integrationStep < 0)
        {
          std::reverse(epochs.begin(), epochs.end());
          if(catchUpYaw(epoch, epochs, yawRateHardwareMax, FALSE))
            return FALSE;
        }

        return TRUE;
      }

      yawRatePrevious = epoch.yawRateTrue;
    }

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SimulateStarCameraGnss::yawBias(const Time &time, const MiscValueArc &yawBiasArc)
{
  try
  {
    std::vector<Time> times = yawBiasArc.times();
    for(UInt i = times.size(); i --> 0; )
      if(time >= times.at(i))
        return yawBiasArc.at(i).value;

    throw(Exception("no yaw bias found for time: " + time.dateTimeStr()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SimulateStarCameraGnss::wrapAngle(Double angle)
{
  try
  {
    while(angle >  PI) angle -= 2*PI;
    while(angle < -PI) angle += 2*PI;

    return angle;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
