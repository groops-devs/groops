/***********************************************/
/**
* @file gnssRinexNavigation2OrbitClock.cpp
*
* @brief Convert RINEX navigation file (e.g. broadcast ephemeris) to orbit and clock files.
*
* @author Sebastian Strasser
* @date 2017-02-28
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Evaluates orbit and clock parameters from \href{https://files.igs.org/pub/data/format/rinex305.pdf}{RINEX navigation file} \config{inputfileRinex}
at epochs given by \configClass{timeSeries}{timeSeriesType} and writes them to \configFile{outputfileOrbit}{instrument} and
\configFile{outputfileClock}{instrument}, respectively.

Orbits are rotated from TRF (as broadcasted) to CRF via \configClass{earthRotation}{earthRotationType},
but system-specific TRFs (WGS84, PZ-90, etc.) are not aligned to a common TRF.

See also \program{OrbitAddVelocityAndAcceleration}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Convert RINEX navigation file (e.g. broadcast ephemeris) to orbit and clock files.
* @ingroup programsConversionGroup */
class GnssRinexNavigation2OrbitClock
{
  class Satellite
  {
  public:
    std::vector<Time>   epoch;
    std::vector<Vector> clockParam;
    std::vector<Matrix> orbitParam;
  };

  Double rinexVersion;
  Char   system;
  Vector ionAlpha, ionBeta, deltaUTC;
  Double leapSeconds;
  std::map<GnssType, Satellite> satellites;
  std::vector<Time> times;

  void readHeader(InFile &file, UInt lineCount=MAX_UINT);
  void readData(InFile &file);
  Bool getLine(InFile &file, std::string &line, std::string &label) const;
  Bool testLabel(const std::string &labelInLine, const std::string &label, Bool optional=TRUE) const;
  OrbitEpoch rungeKutta4(const Time &time, const OrbitEpoch &refEpoch, const Vector3d &sunMoonAcceleration, EarthRotationPtr earthRotation) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssRinexNavigation2OrbitClock, SINGLEPROCESS, "Convert RINEX navigation file (e.g. broadcast ephemeris) to orbit and clock files.", Conversion, Gnss, Instrument)

/***********************************************/

void GnssRinexNavigation2OrbitClock::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outNameOrbit, outNameClock, inNameRinex;
    TimeSeriesPtr timeSeriesPtr;
    EarthRotationPtr earthRotation;
    std::vector<std::string> usePrn;

    readConfig(config, "outputfileOrbit", outNameOrbit,  Config::OPTIONAL, "", "PRN is appended to file name");
    readConfig(config, "outputfileClock", outNameClock,  Config::OPTIONAL, "", "PRN is appended to file name");
    readConfig(config, "inputfileRinex",  inNameRinex,   Config::MUSTSET,  "", "RINEX navigation file");
    readConfig(config, "timeSeries",      timeSeriesPtr, Config::MUSTSET,  "", "orbit and clock evaluation epochs");
    readConfig(config, "earthRotation",   earthRotation, Config::MUSTSET,  "", "for rotation from TRF to CRF");
    readConfig(config, "usePrn",          usePrn,        Config::OPTIONAL, "", "only export these PRNs instead of all");
    if(isCreateSchema(config)) return;

    times = timeSeriesPtr->times();

    // read file
    logStatus << "read RINEX file <" << inNameRinex << ">" << Log::endl;
    InFile file(inNameRinex);
    readHeader(file);
    readData(file);
    file.close();

    std::set<GnssType> unsupportedSatellites;
    if(!outNameOrbit.empty())
    {
      logStatus << "write orbit files <" << outNameOrbit.appendBaseName(".***") << ">" << Log::endl;

      std::map<GnssType, Double> GM {{GnssType::GPS,     3.986005e14},
                                     {GnssType::GALILEO, 3.986004418e14},
                                     {GnssType::BDS,     3.986004418e14},
                                     {GnssType::QZSS,    3.986005e14}};
      std::map<GnssType, Double> omega_e {{GnssType::GPS,     7.2921151467e-5},
                                          {GnssType::GALILEO, 7.2921151467e-5},
                                          {GnssType::BDS,     7.292115e-5},
                                          {GnssType::QZSS,    7.2921151467e-5}};
      std::map<GnssType, Time> refTime {{GnssType::GPS,     date2time(1980, 1, 6)},
                                        {GnssType::GALILEO, date2time(1980, 1, 6)},
                                        {GnssType::BDS,     date2time(2006, 1, 1, 0, 0, 14)},
                                        {GnssType::QZSS,    date2time(1980, 1, 6)}};

      for(const auto &satellite : satellites)
      {
        if(satellite.first == GnssType::SBAS)
        {
          unsupportedSatellites.insert(satellite.first);
          continue;
        }

        if(!satellite.second.epoch.size() || (usePrn.size() && std::find(usePrn.begin(), usePrn.end(), satellite.first.str().substr(3,3)) == usePrn.end()))
          continue;

        // fill arc with epochs
        UInt idx = 0;
        OrbitArc arc;
        for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
        {
          while(idx+1 < satellite.second.epoch.size() &&
                std::fabs((times.at(idEpoch)-satellite.second.epoch.at(idx+1)).seconds()) < std::fabs((times.at(idEpoch)-satellite.second.epoch.at(idx)).seconds()))
            idx++;

          OrbitEpoch epoch;
          epoch.time = times.at(idEpoch);
          const_MatrixSliceRef orbitParam = satellite.second.orbitParam.at(idx);

          if(satellite.first == GnssType::GLONASS)
          {
            OrbitEpoch refEpoch;
            refEpoch.time = timeUTC2GPS(satellite.second.epoch.at(idx));
            refEpoch.position = Vector3d(orbitParam(0,0), orbitParam(1,0), orbitParam(2,0))*1e3;
            refEpoch.velocity = Vector3d(orbitParam(0,1), orbitParam(1,1), orbitParam(2,1))*1e3;
            refEpoch.acceleration = Vector3d(orbitParam(0,2), orbitParam(1,2), orbitParam(2,2))*1e3;
            Rotary3d crf2trf  = earthRotation->rotaryMatrix(refEpoch.time);
            Vector3d omega    = earthRotation->rotaryAxis(refEpoch.time);
            refEpoch.position = crf2trf.inverseRotate(refEpoch.position); // TRF -> CRF
            refEpoch.velocity = crf2trf.inverseRotate(refEpoch.velocity) + crossProduct(omega, refEpoch.position); // TRF -> CRF

            // Runge-Kutta-4 integration with 60 second intermediate steps
            const Double dt = (epoch.time - refEpoch.time).seconds();
            const Double integrationStep = 60;
            OrbitEpoch intermediateEpoch = refEpoch;
            for(UInt i = 1; i*integrationStep < std::fabs(dt); i++)
              intermediateEpoch = rungeKutta4(refEpoch.time+seconds2time(i*integrationStep*dt/std::fabs(dt)), intermediateEpoch, refEpoch.acceleration, earthRotation);
            epoch = rungeKutta4(epoch.time, intermediateEpoch, refEpoch.acceleration, earthRotation);
            // NOTE: integration can lead to errors up to ~10 m after 15 minutes, unclear if that's the accuracy limit or there's an issue somewhere in the code
          }
          else // all systems using GPS-like ephemerides
          {
            const Double toe_week = orbitParam(4,2);
            const Double toe_sec  = orbitParam(2,0);
            const Double M0       = orbitParam(0,3);
            const Double a        = std::pow(orbitParam(1,3), 2);
            const Double delta_n  = orbitParam(0,2);
            const Double e        = orbitParam(1,1);
            const Double omega    = orbitParam(3,2);
            const Double cuc      = orbitParam(1,0);
            const Double cus      = orbitParam(1,2);
            const Double cic      = orbitParam(2,1);
            const Double cis      = orbitParam(2,3);
            const Double crc      = orbitParam(3,1);
            const Double crs      = orbitParam(0,1);
            const Double i0       = orbitParam(3,0);
            const Double iDot     = orbitParam(4,0);
            const Double Omega0   = orbitParam(2,2);
            const Double OmegaDot = orbitParam(3,3);

            // time difference with respect to time of ephemerides (t_oe)
            const Double toe = toe_week * 604800 + toe_sec;
            Double dt = (epoch.time - refTime[satellite.first & GnssType::SYSTEM]).seconds() - toe;
            while(dt > 302400)
              dt -= 604800;
            while(dt < -302400)
              dt += 604800;

            // mean, eccentric, and true anomaly [rad]
            const Double M = M0 + (std::sqrt(GM[satellite.first & GnssType::SYSTEM]/std::pow(a, 3)) + delta_n) * dt;
            Double E = M;
            for(UInt i = 0; i < 10; i++)
              E = M + e * std::sin(E);
            const Double nu = std::atan2(std::sqrt(1.-std::pow(e, 2)) * std::sin(E), std::cos(E) - e);

            // argument of latitude, radial distance, and inclination
            const Double u = omega + nu + cuc * std::cos(2*(omega+nu)) + cus * std::sin(2*(omega+nu));
            const Double r = a * (1.-e*std::cos(E)) + crc * std::cos(2*(omega+nu)) + crs * std::sin(2*(omega+nu));
            const Double i = i0 + iDot*dt + cic * std::cos(2*(omega+nu)) + cis * std::sin(2*(omega+nu));

            // longitude of ascending node and coordinates
            if(satellite.first == GnssType::BDS && (satellite.first.prn() <= 5 || (satellite.first.prn() >= 59 && satellite.first.prn() <= 63))) // BDS GEO (C01-C05, C59-C63) satellites
            {
              const Double lambda = Omega0 + OmegaDot*dt - omega_e[satellite.first & GnssType::SYSTEM]*toe_sec;
              epoch.position = (rotaryZ(Angle(-lambda))*rotaryX(Angle(-i))*rotaryZ(Angle(-u))).rotate(Vector3d(r, 0, 0));
              epoch.position = (rotaryZ(Angle(omega_e[satellite.first & GnssType::SYSTEM]*dt))*rotaryX(Angle(-5*DEG2RAD))).rotate(epoch.position);
            }
            else // all MEO and IGSO satellites
            {
              const Double lambda = Omega0 + (OmegaDot - omega_e[satellite.first & GnssType::SYSTEM])*dt - omega_e[satellite.first & GnssType::SYSTEM]*toe_sec;
              epoch.position = (rotaryZ(Angle(-lambda))*rotaryX(Angle(-i))*rotaryZ(Angle(-u))).rotate(Vector3d(r, 0, 0));
            }
            epoch.position = earthRotation->rotaryMatrix(epoch.time).inverseRotate(epoch.position); // TRF -> CRF
          }
          arc.push_back(epoch);
        }

        // write clock file
        std::list<Arc> arcList;
        arcList.push_back(arc);
        InstrumentFile::write(outNameOrbit.appendBaseName('.'+satellite.first.str().substr(3,3)), arcList);
      }

    }

    if(!outNameClock.empty())
    {
      logStatus << "write clock files <" << outNameClock.appendBaseName(".***") << ">" << Log::endl;
      for(const auto &satellite : satellites)
      {
        if(satellite.first == GnssType::SBAS)
        {
          unsupportedSatellites.insert(satellite.first);
          continue;
        }

        if(!satellite.second.epoch.size() || (usePrn.size() && std::find(usePrn.begin(), usePrn.end(), satellite.first.str().substr(3,3)) == usePrn.end()))
          continue;

        // fill arc with epochs
        UInt idx = 0;
        MiscValueArc arc;
        for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
        {
          while(idx+1 < satellite.second.epoch.size() && times.at(idEpoch) >= satellite.second.epoch.at(idx+1))
            idx++;

          MiscValueEpoch epoch;
          epoch.time = times.at(idEpoch);
          const Double dt = (times.at(idEpoch) - satellite.second.epoch.at(idx)).seconds();
          for(UInt i = 0; i < satellite.second.clockParam.at(idx).size(); i++)
            epoch.value += satellite.second.clockParam.at(idx)(i)*std::pow(dt, i);  // a_i*dt^i
          arc.push_back(epoch);
        }

        // write clock file
        std::list<Arc> arcList;
        arcList.push_back(arc);
        InstrumentFile::write(outNameClock.appendBaseName('.'+satellite.first.str().substr(3,3)), arcList);
      }
    }

    if(unsupportedSatellites.size())
    {
      std::stringstream ss;
      for(auto sat : unsupportedSatellites)
        ss << " " << sat.prnStr();
      logWarning << "conversion not implemented yet for these satellites:" << ss.str() << Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssRinexNavigation2OrbitClock::readHeader(InFile &file, UInt lineCount)
{
  try
  {
    std::string line, label;
    getLine(file, line, label);
    testLabel(label, "RINEX VERSION / TYPE", FALSE);
    rinexVersion    = String::toDouble(line.substr(0, 9));
    if(rinexVersion<2)
      logWarning << "old RINEX version: " << rinexVersion << Log::endl;
    if(line.at(20) != 'N' && !(rinexVersion < 3 && line.at(20) == 'G'))
      throw(Exception("File must contain Navigation Data"));
    system = rinexVersion < 3 ? (line.at(20) == 'N' ? 'G' : 'R') : line.at(40);

    for(UInt idLine=0; idLine<lineCount; idLine++)
    {
      if(!getLine(file, line, label))
        throw(Exception("error while reading RINEX header"));
      if(std::all_of(line.begin(), line.end(), isspace))
      {
        if(rinexVersion < 2)
          break;
        else
          continue;
      }
      // ====================================
      if(testLabel(label, "END OF HEADER"))
        break;
      // ====================================
      else if(testLabel(label, "PGM / RUN BY / DATE"))
      {
      }
      // ====================================
      else if(testLabel(label, "COMMENT"))
      {
      }
      // ====================================
      else if(testLabel(label, "ION ALPHA"))
      {
        ionAlpha = Vector(4);
        for(UInt i = 0; i < ionAlpha.size(); i++)
          ionAlpha(i) = String::toDouble(line.substr(2+i*12, 12));
      }
      // ====================================
      else if(testLabel(label, "ION BETA"))
      {
        ionBeta = Vector(4);
        for(UInt i = 0; i < ionBeta.size(); i++)
          ionBeta(i) = String::toDouble(line.substr(2+i*12, 12));
      }
      // ====================================
      else if(testLabel(label, "DELTA-UTC: A0,A1,T,W"))
      {
        deltaUTC = Vector(4);
        deltaUTC(0) = String::toDouble(line.substr(3, 19));
        deltaUTC(1) = String::toDouble(line.substr(22, 19));
        deltaUTC(2) = String::toInt(line.substr(41, 9));
        deltaUTC(3) = String::toInt(line.substr(50, 9));
      }
      // ====================================
      else if(testLabel(label, "LEAP SECONDS"))
      {
        leapSeconds = String::toInt(line.substr(0, 6));
      }
      // ====================================
      else if(testLabel(label, "CORR TO SYSTEM TIME"))
      {
      }
      // ====================================
      else if(testLabel(label, "TIME SYSTEM CORR"))
      {
      }
      // ====================================
      else if(testLabel(label, "IONOSPHERIC CORR"))
      {
      }
      // ====================================
      else
      {
        logWarning<<"Unknown header label:"<<Log::endl;
        logWarning<<"'"<<line<<"'"<<Log::endl;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssRinexNavigation2OrbitClock::readData(InFile &file)
{
  std::string line, label;

  try
  {
    while(getLine(file, line, label))
    {
      std::string prnStr = rinexVersion < 3 ? system+line.substr(0,2): line.substr(0,3);
      if(prnStr.at(1) == ' ') prnStr.at(1) = '0';
      GnssType prn("***" + prnStr);

      const UInt lineCount = (prn == GnssType::GLONASS || prn == GnssType::SBAS) ? 3 : (rinexVersion < 2 ? 6 : 7);

      // epoch
      Int year   = String::toInt(line.substr(rinexVersion < 3 ?  3 :  4, rinexVersion < 3 ? 2 : 4));
      Int month  = String::toInt(line.substr(rinexVersion < 3 ?  6 :  9, 2));
      Int day    = String::toInt(line.substr(rinexVersion < 3 ?  9 : 12, 2));
      Int hour   = String::toInt(line.substr(rinexVersion < 3 ? 12 : 15, 2));
      Int minute = String::toInt(line.substr(rinexVersion < 3 ? 15 : 18, 2));
      Double sec = String::toDouble(line.substr(rinexVersion < 3 ? 17 : 21, rinexVersion < 3 ? 5 : 2));
      if(rinexVersion < 3)
        year += ((year<=80) ? 2000 : 1900);
      const Time time = date2time(year, month, day, hour, minute, sec);

      if(satellites[prn].epoch.size() && time <= satellites[prn].epoch.back())
      {
        //logWarning << prn.str().substr(3,3) << ": duplicate entry at " << time.dateTimeStr() << Log::endl;
        for(UInt i = 0; i < lineCount; i++)
          getLine(file, line, label);
        continue;
      }
      satellites[prn].epoch.push_back(time);

      // clock polynomial
      Vector clockParam((prn != GnssType::GLONASS && prn != GnssType::SBAS) ? 3 : 2);
      for(UInt i = 0; i < clockParam.size(); i++)
        clockParam(i) = String::toDouble(line.substr((rinexVersion < 3 ? 22 : 23)+i*19, 19));
      satellites[prn].clockParam.push_back(clockParam);

      // orbit parameters
      Matrix orbitParam(lineCount,4);
      for(UInt i = 0; i < lineCount; i++)
      {
        getLine(file, line, label);
        for(UInt j = 0; j < 4; j++)
          orbitParam(i,j) = String::toDouble(line.substr((rinexVersion < 3 ? 3 : 4)+j*19, 19));
      }
      satellites[prn].orbitParam.push_back(orbitParam);
    }
  }
  catch(std::exception &e)
  {
    logError<<"'"<<line<<"'"<<Log::endl;
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssRinexNavigation2OrbitClock::getLine(InFile &file, std::string &line, std::string &label) const
{
  try
  {
    getline(file, line);
    if(line.back() == '\r')
      line.pop_back();
    if(line.size()<80)
      line.resize(80,' ');
    label = line.substr(60,20);
    return file.good();
  }
  catch(...)
  {
    line.clear();
    line.resize(80,' ');
    label = line.substr(60,20);
    return FALSE;
  }
}

/***********************************************/

Bool GnssRinexNavigation2OrbitClock::testLabel(const std::string &labelInLine, const std::string &label, Bool optional) const
{
  if(labelInLine.find(label)!=std::string::npos)
    return TRUE;
  if(optional)
    return FALSE;
  throw(Exception(std::string("In Line '")+labelInLine+"' label '"+label+"' expected\n"));
}

/***********************************************/

OrbitEpoch GnssRinexNavigation2OrbitClock::rungeKutta4(const Time &time, const OrbitEpoch &refEpoch, const Vector3d &sunMoonAcceleration, EarthRotationPtr earthRotation) const
{
  try
  {
    const Double dt = (time - refEpoch.time).seconds();
    const Double GM = 3.9860044e14;
    const Double C20 = -1082.63e-6;
    const Double a_e = 6378136;

    auto acceleration = [&](const OrbitEpoch &epoch)
    {
      Rotary3d crf2trf = earthRotation->rotaryMatrix(epoch.time);
      Double r = epoch.position.r();
      Vector3d ePos = crf2trf.rotate(epoch.position)/r;
      Double mu = GM/std::pow(r, 2);
      Double x = -mu*ePos.x() + 3/2*C20*mu*ePos.x()*std::pow(a_e/r, 2)*(1.-5.*std::pow(ePos.z(), 2)) + sunMoonAcceleration.x();
      Double y = -mu*ePos.y() + 3/2*C20*mu*ePos.y()*std::pow(a_e/r, 2)*(1.-5.*std::pow(ePos.z(), 2)) + sunMoonAcceleration.y();
      Double z = -mu*ePos.z() + 3/2*C20*mu*ePos.z()*std::pow(a_e/r, 2)*(3.-5.*std::pow(ePos.z(), 2)) + sunMoonAcceleration.z();
      return crf2trf.inverseRotate(Vector3d(x, y, z));
    };

    // orbit integration with Runge-Kutta-4 algorithm
    OrbitEpoch k1 = refEpoch;
    k1.acceleration = acceleration(refEpoch);

    OrbitEpoch k2 = k1;
    k2.time        += seconds2time(dt/2.);
    k2.position    += dt/2. * k1.velocity;
    k2.velocity    += dt/2. * k1.acceleration;
    k2.acceleration = acceleration(k2);

    OrbitEpoch k3 = k1;
    k3.time        += seconds2time(dt/2.);
    k3.position    += dt/2. * k2.velocity;
    k3.velocity    += dt/2. * k2.acceleration;
    k3.acceleration = acceleration(k3);

    OrbitEpoch k4 = k1;
    k4.time        += seconds2time(dt);
    k4.position    += dt * k3.velocity;
    k4.velocity    += dt * k3.acceleration;
    k4.acceleration = acceleration(k4);

    // Compute final value for this epoch
    OrbitEpoch epoch = k1;
    epoch.time        += seconds2time(dt);
    epoch.position    += (dt/6.) * (k1.velocity + 2*k2.velocity + 2*k3.velocity + k4.velocity);
    epoch.velocity    += (dt/6.) * (k1.acceleration + 2*k2.acceleration + 2*k3.acceleration + k4.acceleration);
    epoch.acceleration = acceleration(epoch);

    return epoch;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
