/***********************************************/
/**
* @file gnssRinexNavigation2OrbitClock.cpp
*
* @brief Convert RINEX navigation file (e.g. broadcast ephemeris) to orbit and clock files.
*
* @author Sebastian Strasser
* @author Patrick Dumitraschkewitz
* @date 2017-02-28
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Evaluates orbit and clock parameters from \href{https://files.igs.org/pub/data/format/rinex_4.00.pdf}{RINEX} (version 2, 3, and 4)
navigation file \config{inputfileRinex} at epochs given by \configClass{timeSeries}{timeSeriesType} and writes them to
\configFile{outputfileOrbit}{instrument} and \configFile{outputfileClock}{instrument}, respectively.

Orbits are rotated from TRF (as broadcasted) to CRF via \configClass{earthRotation}{earthRotationType},
but system-specific TRFs (WGS84, PZ-90, etc.) are not aligned to a common TRF.

If \config{messageType} is set (e.g., to LNAV, CNAV, or other types defined in the RINEX 4 standard), only navigation records of this type are used.
Otherwise, if multiple records are defined for the same epoch, the first one is used.

Also take note that all orbits are written out even satellites whose health flag would suggest otherwise.

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
  /**
   * @brief The NavData class will be based on the LNAV nav data message.
   */
  class NavData
  {
  public:
    NavData(const Double &rinexVersion,
            const GnssType prn,
            const Time &time_clock,
            Vector clockPoly,
            Matrix orbitData,
            const Time &gnssT0,
            std::string msg) :  msgType(msg), clockParam(clockPoly), tocValue(time_clock), data(orbitData), sat(prn)
    {
      parseData(rinexVersion, prn, gnssT0);
    }

    NavData(const GnssType prn, Vector clockPoly, const Time &time_clock, Matrix orbitData, std::string msg) : msgType(msg), clockParam(clockPoly),
                                                                                                               tocValue(time_clock), data(orbitData), sat(prn)
    {

    }

    const std::string msgType;


    Time   toc()      const  {return tocValue;}
    Double toe()      const  {return toeValue;}
    Double toe_week() const  {return toeWeekValue;}
    Double toe_sec()  const  {return toeSecValue;}
    Double crs()      const  {checkIfKeplerEphem(); return data(0, 1);}
    Double delta_n()  const  {checkIfKeplerEphem(); return data(0, 2);}
    Double m0()       const  {checkIfKeplerEphem(); return data(0, 3);}
    Double cuc()      const  {checkIfKeplerEphem(); return data(1, 0);}
    Double e()        const  {checkIfKeplerEphem(); return data(1, 1);}
    Double cus()      const  {checkIfKeplerEphem(); return data(1, 2);}
    Double a()        const  {checkIfKeplerEphem(); return std::pow(data(1, 3), 2.);}
    Double cic()      const  {checkIfKeplerEphem(); return data(2, 1);}
    Double omega0()   const  {checkIfKeplerEphem(); return data(2, 2);}
    Double cis()      const  {checkIfKeplerEphem(); return data(2, 3);}
    Double i0()       const  {checkIfKeplerEphem(); return data(3, 0);}
    Double crc()      const  {checkIfKeplerEphem(); return data(3, 1);}
    Double omega()    const  {checkIfKeplerEphem(); return data(3, 2);}
    Double omegaDot() const  {checkIfKeplerEphem(); return data(3, 3);}
    Double iDot()     const  {checkIfKeplerEphem(); return data(4, 0);}

    Double posX()     const  {checkIfStateEphem(); return data(0, 0);}
    Double velX()     const  {checkIfStateEphem(); return data(0, 1);}
    Double accX()     const  {checkIfStateEphem(); return data(0, 2);}
    Double posY()     const  {checkIfStateEphem(); return data(1, 0);}
    Double velY()     const  {checkIfStateEphem(); return data(1, 1);}
    Double accY()     const  {checkIfStateEphem(); return data(1, 2);}
    Double posZ()     const  {checkIfStateEphem(); return data(2, 0);}
    Double velZ()     const  {checkIfStateEphem(); return data(2, 1);}
    Double accZ()     const  {checkIfStateEphem(); return data(2, 2);}

    Vector clockParam;

  private:
    void time2GnssWeekSecond(const Time &t0, const Time &time, Double &week, Double &sec) const;


    void parseData(Double rinexVersion, GnssType prn, const Time &gnssT0);

    void checkIfKeplerEphem() const {if(sat == GnssType::GLONASS || sat == GnssType::SBAS) throw std::runtime_error("Ephemeris not kepler parametrizised! Check usage of access methods!");}
    void checkIfStateEphem() const {if(sat != GnssType::GLONASS && sat != GnssType::SBAS) throw std::runtime_error("Ephemeris are not state vectors! Check usage of access methods!");}

    Double toeValue = 0.;
    Double toeWeekValue = 0.;
    Double toeSecValue = 0.;
    Time tocValue;
    Matrix data;
    GnssType sat;
  };


  const std::map<std::string, UInt> message2LineCount = {{"LNAV", 7},
                                                         {"CNAV", 8},
                                                         {"CNV1", 9},
                                                         {"CNV2", 9},
                                                         {"CNV3", 8},
                                                         {"FNAV", 7},
                                                         {"INAV", 7},
                                                         {"D1",   7},
                                                         {"D2",   7},
                                                         {"FDMA", 4},
                                                         {"SBAS", 3}};


  // ref times in GPS time thats why BDS is +14 seconds
  const std::map<GnssType, Time> refTime {{GnssType::GPS,     date2time(1980, 1, 6)},
                                          {GnssType::GALILEO, date2time(1980, 1, 6)},
                                          {GnssType::BDS,     date2time(2006, 1, 1, 0, 0, 14)},
                                          {GnssType::QZSS,    date2time(1980, 1, 6)},
                                          {GnssType::IRNSS,   date2time(1980, 1, 6)}};

  // GM of the different reference systems
  const std::map<GnssType, Double> GM {{GnssType::GPS,     3.986005e14},
                                       {GnssType::GALILEO, 3.986004418e14},
                                       {GnssType::BDS,     3.986004418e14},
                                       {GnssType::QZSS,    3.986005e14},
                                       {GnssType::IRNSS,   3.986005e14},
                                       {GnssType::GLONASS, 3.9860044e14}};

  const std::map<GnssType, Double> omega_e {{GnssType::GPS,     7.2921151467e-5},
                                            {GnssType::GALILEO, 7.2921151467e-5},
                                            {GnssType::BDS,     7.292115e-5},
                                            {GnssType::QZSS,    7.2921151467e-5},
                                            {GnssType::IRNSS,   7.2921151467e-5}};

  // Parameters for GNSS that use state vectors instead of kepler elements
  const std::map<GnssType, Double> C20 {{GnssType::GLONASS, -1082.63e-6}};
  const std::map<GnssType, Double> a_e {{GnssType::GLONASS, 6378136}};


  const std::vector<GnssType> notImplementedSystems = {GnssType::SBAS};

  Double rinexVersion;
  Char   system;
  Vector ionAlpha, ionBeta, deltaUTC;
  Double leapSeconds;
  std::map<GnssType, std::vector<NavData>> satellites;
  std::vector<Time> times;
  std::string useMessageType;

  void readHeader(InFile &file, UInt lineCount=MAX_UINT);
  void readData(InFile &file, const std::vector<std::string> &usePrn, std::set<GnssType> &unsupported);
  Bool getLine(InFile &file, std::string &line, std::string &label) const;
  Bool testLabel(const std::string &labelInLine, const std::string &label, Bool optional=TRUE) const;
  OrbitEpoch rungeKutta4(const Time &time, const OrbitEpoch &refEpoch, const Vector3d &sunMoonAcceleration, EarthRotationPtr earthRotation, const Double GM, const Double C20, const Double a_e) const;

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

    readConfig(config, "outputfileOrbit", outNameOrbit,   Config::OPTIONAL, "", "PRN is appended to file name");
    readConfig(config, "outputfileClock", outNameClock,   Config::OPTIONAL, "", "PRN is appended to file name");
    readConfig(config, "inputfileRinex",  inNameRinex,    Config::MUSTSET,  "", "RINEX navigation file");
    readConfig(config, "timeSeries",      timeSeriesPtr,  Config::MUSTSET,  "", "orbit and clock evaluation epochs");
    readConfig(config, "earthRotation",   earthRotation,  Config::MUSTSET,  "", "for rotation from TRF to CRF");
    readConfig(config, "usePrn",          usePrn,         Config::OPTIONAL, "", "only export these PRNs instead of all");
    readConfig(config, "messageType",     useMessageType, Config::OPTIONAL, "", "(RINEX4) only use this navigation message (LNAV, CNAV, ...)");
    if(isCreateSchema(config)) return;

    times = timeSeriesPtr->times();

    // read file
    logStatus << "read RINEX file <" << inNameRinex << ">" << Log::endl;
    InFile file(inNameRinex);
    logStatus << "read RINEX Header file <" << inNameRinex << ">" << Log::endl;

    readHeader(file);
    logStatus << "read RINEX Data file <" << inNameRinex << ">" << Log::endl;
    std::set<GnssType> unsupportedSatellites;

    readData(file, usePrn, unsupportedSatellites);
    file.close();

    if(!outNameOrbit.empty())
    {
      for(const auto &satellite : satellites)
      {
        if(satellite.first.isInList(notImplementedSystems))
        {
          unsupportedSatellites.insert(satellite.first);
          continue;
        }

        if(!satellite.second.size())
          continue;
        logStatus << "writing orbit file <" << outNameOrbit.appendBaseName("."+satellite.first.prnStr()) << ">" << Log::endl;

        // fill arc with epochs
        OrbitArc arc;
        for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
        {
          UInt selection = NULLINDEX;
          Double currentDt = 0.;
          for(UInt idx = 0; idx < satellite.second.size(); idx++)
          {

            Double dt = 0.;
            if(satellite.first == GnssType::GLONASS)
              dt = (times.at(idEpoch) - satellite.second.at(idx).toc()).seconds();
            else
              dt = (times.at(idEpoch) - refTime.at(satellite.first & GnssType::SYSTEM)).seconds() - satellite.second.at(idx).toe();

            if(selection == NULLINDEX || fabs(dt) < currentDt)
            {
              selection = idx;
              currentDt = fabs(dt);
            }
          }

          OrbitEpoch epoch;
          epoch.time = times.at(idEpoch);
          const NavData &epochNavData = satellite.second.at(selection);


          if(satellite.first == GnssType::GLONASS)
          {
            OrbitEpoch refEpoch;
            refEpoch.time = timeUTC2GPS(epochNavData.toc());
            refEpoch.position     = Vector3d(epochNavData.posX(), epochNavData.posY(), epochNavData.posZ()) * 1.e3;
            refEpoch.velocity     = Vector3d(epochNavData.velX(), epochNavData.velY(), epochNavData.velZ()) * 1.e3;
            refEpoch.acceleration = Vector3d(epochNavData.accX(), epochNavData.accY(), epochNavData.accZ()) * 1.e3;
            Rotary3d crf2trf  = earthRotation->rotaryMatrix(refEpoch.time);
            Vector3d omega    = earthRotation->rotaryAxis(refEpoch.time);
            refEpoch.position = crf2trf.inverseRotate(refEpoch.position); // TRF -> CRF
            refEpoch.velocity = crf2trf.inverseRotate(refEpoch.velocity) + crossProduct(omega, refEpoch.position); // TRF -> CRF

            // Runge-Kutta-4 integration with 60 second intermediate steps
            const Double dt = (epoch.time - refEpoch.time).seconds();
            const Double integrationStep = 60;
            OrbitEpoch intermediateEpoch = refEpoch;
            for(UInt i = 1; i*integrationStep < std::fabs(dt); i++)
              intermediateEpoch = rungeKutta4(refEpoch.time+seconds2time(i*integrationStep*dt/std::fabs(dt)), intermediateEpoch, refEpoch.acceleration, earthRotation, GM.at(GnssType::GLONASS), C20.at(GnssType::GLONASS), a_e.at(GnssType::GLONASS));
            epoch = rungeKutta4(epoch.time, intermediateEpoch, refEpoch.acceleration, earthRotation, GM.at(GnssType::GLONASS), C20.at(GnssType::GLONASS), a_e.at(GnssType::GLONASS));
            // NOTE: integration can lead to errors up to ~10 m after 15 minutes, unclear if that's the accuracy limit or there's an issue somewhere in the code
          }
          else // all systems using GPS-like ephemerides
          {
            Double dt = (epoch.time - refTime.at(satellite.first & GnssType::SYSTEM)).seconds() - epochNavData.toe();
            while(dt > 302400)
              dt -= 604800;
            while(dt < -302400)
              dt += 604800;

            const Double n = std::sqrt(GM.at(satellite.first & GnssType::SYSTEM) / std::pow(epochNavData.a(), 3)) + epochNavData.delta_n();

            const Double M = std::fmod(epochNavData.m0() + n * dt, 2. * PI);

            Double E = M;
            for(UInt i = 0; i < 10; i++)
              E = M + epochNavData.e() * std::sin(E);

            E = std::fmod(E, 2. * PI);

            // ture anomaly
            const Double nu     = std::atan2(std::sqrt(1. - std::pow(epochNavData.e(), 2.)) * std::sin(E), std::cos(E) - epochNavData.e());

            //argument of latitude
            const Double du = epochNavData.cuc() * std::cos(2. * (nu + epochNavData.omega())) + epochNavData.cus() * std::sin(2. * (nu + epochNavData.omega()));
            const Double u = nu + epochNavData.omega() + du;

            // radius in orbital plane
            const Double dr = epochNavData.crc() * std::cos(2. * (nu + epochNavData.omega())) + epochNavData.crs() * std::sin(2. * (nu + epochNavData.omega()));
            const Double r = epochNavData.a() * (1. - epochNavData.e() * std::cos(E)) + dr;

            // orbit states
            const Double xOrb = r * std::cos(u);
            const Double yOrb = r * std::sin(u);

            // inclination
            const Double di = epochNavData.cic() * std::cos(2. * (nu + epochNavData.omega())) + epochNavData.cis() * std::sin(2. * (nu + epochNavData.omega()));
            const Double i = epochNavData.i0() + epochNavData.iDot() * dt + di;

            // longitude of ascending node
            const Double lambda = epochNavData.omega0()
                                 + (epochNavData.omegaDot() - omega_e.at(satellite.first & GnssType::SYSTEM)) * dt
                                 - omega_e.at(satellite.first & GnssType::SYSTEM) * epochNavData.toe_sec();


            // Position
            epoch.position.x() = xOrb * std::cos(lambda) - yOrb * std::sin(lambda) * std::cos(i);
            epoch.position.y() = xOrb * std::sin(lambda) + yOrb * std::cos(lambda) * std::cos(i);
            epoch.position.z() = yOrb * std::sin(i);

            //PD: No idea if the following is sitll relevant for geo satellites.
//            if(satellite.first == GnssType::BDS && (satellite.first.prn() <= 5 || (satellite.first.prn() >= 59 && satellite.first.prn() <= 63)))
            // BDS GEO (C01-C05, C59-C63) satellites
//            {
//              const Double lambda = Omega0 + OmegaDot*dt - omega_e[satellite.first & GnssType::SYSTEM]*toe_sec;
//              epoch.position = (rotaryZ(Angle(-lambda))*rotaryX(Angle(-i))*rotaryZ(Angle(-u))).rotate(Vector3d(r, 0, 0));
//              epoch.position = (rotaryZ(Angle(omega_e[satellite.first & GnssType::SYSTEM]*dt))*rotaryX(Angle(-5*DEG2RAD))).rotate(epoch.position);
//            }

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
      for(const auto &satellite : satellites)
      {
        if(satellite.first.isInList(notImplementedSystems))
        {
          unsupportedSatellites.insert(satellite.first);
          continue;
        }


        if(!satellite.second.size())
          continue;

        logStatus << "writing clock file <" << outNameOrbit.appendBaseName("."+satellite.first.prnStr()) << ">" << Log::endl;

        // fill arc with epochs
        MiscValueArc arc;
        for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
        {
          UInt selection = NULLINDEX;
          Double currentDt = 0.;
          for(UInt idx = 0; idx < satellite.second.size(); idx++)
          {
            Double dt = (times.at(idEpoch) - mjd2time(satellite.second.at(idx).toe())).seconds();
            if(selection == NULLINDEX || dt < currentDt)
              selection = idx;
          }

          const NavData &epochNavData = satellite.second.at(selection);

          MiscValueEpoch epoch;
          epoch.time = times.at(idEpoch);
          const Double dt = (times.at(idEpoch) - epochNavData.toc()).seconds();
          for(UInt i = 0; i < epochNavData.clockParam.size(); i++)
            epoch.value += epochNavData.clockParam(i)*std::pow(dt, i);

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
      else if(testLabel(label, "DOI"))
      {
      }
      // ====================================
      else if(testLabel(label, "LICENCE OF USE"))
      {
      }
      // ====================================
      else if(testLabel(label, "STATION INFORMATION"))
      {
      }
      // ====================================
      else if(testLabel(label, "REC # / TYPE / VERS"))
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

void GnssRinexNavigation2OrbitClock::readData(InFile &file, const std::vector<std::string> &usePrn,
                                              std::set<GnssType> &unsupported)
{
  std::string line, label;

  try
  {

    while(getLine(file, line, label))
    {
      // RINEX 4 record identifier
      std::string messageType;
      if(rinexVersion>=4)
      {
        if(line.substr(0,5) != "> EPH") // only ephemeris records are of interest
          continue;
        messageType = String::trim(line.substr(10,4));
        if(useMessageType.size() && messageType != useMessageType) // message type filter
          continue;
        getLine(file, line, label);
      }

      std::string prnStr = rinexVersion < 3 ? system+line.substr(0,2): line.substr(0,3);
      if(prnStr.at(1) == ' ') prnStr.at(1) = '0';
      GnssType prn("***" + prnStr);
      UInt lineCount = (prn == GnssType::GLONASS || prn == GnssType::SBAS) ? 3 : (rinexVersion < 2 ? 6 : 7);
      if(rinexVersion >= 4)
        lineCount = message2LineCount.at(messageType);

      if((usePrn.size() && std::find(usePrn.begin(), usePrn.end(), prn.str().substr(3,3)) == usePrn.end())
         ||  ((prn & GnssType::SYSTEM).isInList(notImplementedSystems)))
      {
        for(UInt i = 0; i < lineCount; i++)
          getLine(file, line, label);

        if((prn & GnssType::SYSTEM).isInList(notImplementedSystems))
          unsupported.insert(prn);
        continue;
      }

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

      // clock polynomial
      Vector clockParam((prn != GnssType::GLONASS && prn != GnssType::SBAS) ? 3 : 2);
      for(UInt i = 0; i < clockParam.size(); i++)
        clockParam(i) = String::toDouble(line.substr((rinexVersion < 3 ? 22 : 23)+i*19, 19));

      // orbit parameters
      Matrix orbitParam(lineCount,4);
      for(UInt i = 0; i < lineCount; i++)
      {
        getLine(file, line, label);
        for(UInt j = 0; j < 4; j++)
          orbitParam(i,j) = String::toDouble(line.substr((rinexVersion < 3 ? 3 : 4)+j*19, 19));
      }

      if(satellites.count(prn) == 0)
        satellites.insert({prn, std::vector<NavData>()});

      if(prn == GnssType::GLONASS)
        satellites.at(prn).push_back(NavData(prn, clockParam, time, orbitParam,messageType));
      else
        satellites.at(prn).push_back(NavData(rinexVersion, prn, time, clockParam, orbitParam, refTime.at(prn & GnssType::SYSTEM), messageType));
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

OrbitEpoch GnssRinexNavigation2OrbitClock::rungeKutta4(const Time &time, const OrbitEpoch &refEpoch, const Vector3d &sunMoonAcceleration, EarthRotationPtr earthRotation,
                                                       const Double GM, const Double C20, const Double a_e) const
{
  try
  {
    const Double dt = (time - refEpoch.time).seconds();


    auto acceleration = [&](const OrbitEpoch &epoch)
    {
      Rotary3d crf2trf = earthRotation->rotaryMatrix(epoch.time);
      Double r = epoch.position.r();
      Vector3d ePos = crf2trf.rotate(epoch.position)/r;
      Double mu = GM/std::pow(r, 2);
      Double x = -mu*ePos.x() + 3./2.*C20*mu*ePos.x()*std::pow(a_e/r, 2)*(1.-5.*std::pow(ePos.z(), 2)) + sunMoonAcceleration.x();
      Double y = -mu*ePos.y() + 3./2.*C20*mu*ePos.y()*std::pow(a_e/r, 2)*(1.-5.*std::pow(ePos.z(), 2)) + sunMoonAcceleration.y();
      Double z = -mu*ePos.z() + 3./2.*C20*mu*ePos.z()*std::pow(a_e/r, 2)*(3.-5.*std::pow(ePos.z(), 2)) + sunMoonAcceleration.z();
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

void GnssRinexNavigation2OrbitClock::NavData::time2GnssWeekSecond(const Time &t0, const Time &time, Double &gpsWeek, Double &gpsSecond) const
{
  try
  {
    gpsWeek = (time.mjdInt() - t0.mjdInt()) / 7;
    gpsSecond = std::fmod(time.mjd() - t0.mjd(), 7) * 86400;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssRinexNavigation2OrbitClock::NavData::parseData(Double rinexVersion, GnssType prn, const Time &gnssT0)
{
  try
  {
    if(rinexVersion < 4 && (prn != GnssType::GLONASS && prn != GnssType::SBAS))
    {
      toeSecValue = data(2, 0);
      toeWeekValue = data(4, 2);
      toeValue = data(4, 2) * 604800 + data(2, 0);
      return;
    }

    if(rinexVersion >= 4. && (prn == GnssType::GPS || prn == GnssType::BDS || prn == GnssType::IRNSS || prn == GnssType::GALILEO || prn == GnssType::QZSS)
       && (msgType.compare("LNAV") == 0 || msgType.compare("D1") == 0 || msgType.compare("D2") == 0 || msgType.compare("FNAV") == 0 ||
       msgType.compare("INAV") == 0))
    {
      toeSecValue = data(2, 0);
      toeWeekValue = data(4, 2);
      toeValue = data(4, 2) * 604800. + data(2, 0);
      return;
    }

    if(rinexVersion >= 4. && (prn == GnssType::GPS || prn == GnssType::QZSS)
       && (msgType.compare("CNAV") == 0 || msgType.compare("CNV1") == 0 || msgType.compare("CNV2") == 0))
    {
      time2GnssWeekSecond(gnssT0, tocValue, toeWeekValue, toeSecValue);
      return;
    }


    if(rinexVersion >= 4. && (prn == GnssType::BDS)
       && (msgType.compare("CNV1") == 0 || msgType.compare("CNV2") == 0 || msgType.compare("CNV3") == 0))
    {
      Double toe_sec_tmp;
      time2GnssWeekSecond(gnssT0, tocValue, toeWeekValue, toe_sec_tmp);
      toeSecValue = data(2, 0);
      return;
    }

    logWarning << prn.str() << " " << rinexVersion << " " << msgType << " has no toe computation assigned! Orbit may be faulty!" << Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
