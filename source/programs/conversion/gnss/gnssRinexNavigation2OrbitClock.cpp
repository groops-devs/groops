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

Furthermore, option is available to remove any satellite ephemeris data that has their satellite flag set to unhealthy.

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
  class NavData
  {
  public:
    NavData(const std::string &msgType, GnssType prn, const Time &timeClock, const Vector &clock, const std::vector<Vector> &data)
      : msgType(msgType), prn(prn), timeClock(timeClock), clock(clock), data(data) {}

    std::string         msgType;     // since RINEX 4
    GnssType            prn;
    Time                timeClock, timeEph; // time of clock and ephemeris in GPS time
    Vector              clock;
    std::vector<Vector> data;        // lines with 4 values

    Double crs()      const  {return data.at(0)(1);}
    Double delta_n()  const  {return data.at(0)(2);}
    Double m0()       const  {return data.at(0)(3);}
    Double cuc()      const  {return data.at(1)(0);}
    Double e()        const  {return data.at(1)(1);}
    Double cus()      const  {return data.at(1)(2);}
    Double a()        const  {return std::pow(data.at(1)(3), 2);}
    Double cic()      const  {return data.at(2)(1);}
    Double omega0()   const  {return data.at(2)(2);}
    Double cis()      const  {return data.at(2)(3);}
    Double i0()       const  {return data.at(3)(0);}
    Double crc()      const  {return data.at(3)(1);}
    Double omega()    const  {return data.at(3)(2);}
    Double omegaDot() const  {return data.at(3)(3);}
    Double iDot()     const  {return data.at(4)(0);}
  };

  void readRinex(const FileName &fileName, std::map<GnssType, std::vector<NavData>> &satellites);
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
    FileName                 outNameOrbit, outNameClock, inNameRinex;
    TimeSeriesPtr            timeSeriesPtr;
    EarthRotationPtr         earthRotation;
    std::vector<GnssType>    useType, ignoreType;
    Bool                     removeUnhealthySatellites = TRUE;

    readConfig(config, "outputfileOrbit",           outNameOrbit,              Config::OPTIONAL, "",  "PRN is appended to file name");
    readConfig(config, "outputfileClock",           outNameClock,              Config::OPTIONAL, "",  "PRN is appended to file name");
    readConfig(config, "inputfileRinex",            inNameRinex,               Config::MUSTSET,  "",  "RINEX navigation file");
    readConfig(config, "timeSeries",                timeSeriesPtr,             Config::MUSTSET,  "",  "orbit and clock evaluation epochs");
    readConfig(config, "earthRotation",             earthRotation,             Config::MUSTSET,  "",  "for rotation from TRF to CRF");
    readConfig(config, "useType",                   useType,                   Config::OPTIONAL, "",  "(e.g. ***G12) only use satellites with PRN that match any of these patterns");
    readConfig(config, "ignoreType",                ignoreType,                Config::OPTIONAL, "",  "(e.g. ***R**) ignore satellites PRN that match any of these patterns");
    readConfig(config, "removeUnhealthySatellites", removeUnhealthySatellites, Config::DEFAULT,  "1", "Remove satellite ephemeris that have their sat flags set to unhealthy");
    if(isCreateSchema(config)) return;

    std::vector<Time> times = timeSeriesPtr->times();

    // read file
    logStatus<<"read RINEX file <"<<inNameRinex<<">"<<Log::endl;
    std::map<GnssType, std::vector<NavData>> satellites;
    readRinex(inNameRinex, satellites);

    // check useType and ignoreType
    // ----------------------------
    for(auto iter=satellites.begin(); iter!=satellites.end();)
    {
      Bool use = (useType.size()==0) ? TRUE : FALSE;
      if(iter->first.isInList(useType))
        use = TRUE;
      if(iter->first.isInList(ignoreType))
        use = FALSE;

      if(!use)
        iter = satellites.erase(iter);
      else
        iter++;
    }

    // remove empty satellites
    // -----------------------
    for(auto iter=satellites.begin(); iter!=satellites.end();)
    {
      if(!iter->second.size())
        iter = satellites.erase(iter);
      else
        iter++;
    }

    // remove unsupported
    // ------------------
    const std::vector<GnssType> supportedSystems{GnssType::GPS, GnssType::GALILEO, GnssType::BDS, GnssType::QZSS, GnssType::IRNSS, GnssType::GLONASS, GnssType::SBAS};
    std::set<GnssType> unsupportedSatellites;
    for(auto iter=satellites.begin(); iter!=satellites.end();)
    {
      if(iter->first.isInList(supportedSystems))
        iter++;
      else
      {
        unsupportedSatellites.insert(iter->first);
        iter = satellites.erase(iter);
      }
    }
    if(unsupportedSatellites.size())
    {
      std::stringstream ss;
      for(auto sat : unsupportedSatellites)
        ss<<" "<<sat.prnStr();
      logWarning<<"conversion not implemented yet for these satellites:"<<ss.str()<<Log::endl;
    }

    // remove unhealthy satellites
    // ---------------------------
    if(removeUnhealthySatellites)
      for(auto iter=satellites.begin(); iter!=satellites.end();)
      {
        Bool unhealthy = FALSE;
        for(auto &navData : iter->second)
        {
          if((navData.prn == GnssType::GPS) || (navData.prn == GnssType::QZSS))
          {
            if((navData.msgType == "CNAV") || (navData.msgType == "CNV2"))
            {
              if(navData.data.at(5)(1) > 0)
                unhealthy = TRUE;
            }
            else if(static_cast<UInt>(navData.data.at(5)(1)) & 32)
              unhealthy = TRUE;
          }
          else if((navData.prn == GnssType::GLONASS) || (navData.prn == GnssType::SBAS))
          {
            if(static_cast<UInt>(navData.data.at(0)(3)) & (1+2+4))
              unhealthy = TRUE;
          }
          else if(navData.prn == GnssType::GALILEO)
          {
            UInt dataSourceBits = static_cast<UInt>(navData.data.at(4)(1)) & (1+2+4);
            UInt svH            = static_cast<UInt>(navData.data.at(5)(1)) & (1+8+64);
            if(((dataSourceBits & 1) && (svH & 1))  || // E1B INAV
               ((dataSourceBits & 2) && (svH & 8))  || // E5a FNAV
               ((dataSourceBits & 4) && (svH & 64)))   // E5b INAV
              unhealthy = TRUE;
          }
          else if(navData.prn == GnssType::BDS)
          {
            if(navData.msgType.empty() || navData.msgType == "D1" || navData.msgType == "D2")
            {
              if(static_cast<UInt>(navData.data.at(5)(1)) & 1)
                unhealthy = TRUE;
            }
            else if(navData.msgType == "CNV1" || navData.msgType == "CNV2")
            {
              if(static_cast<UInt>(navData.data.at(7)(1)) & 1)
                unhealthy = TRUE;
            }
            else if(navData.msgType == "CNV3")
            {
              if(static_cast<UInt>(navData.data.at(6)(1)) & 1)
                unhealthy = TRUE;
            }
          }
          else if(navData.prn == GnssType::IRNSS)
          {
            if(navData.data.at(5)(1) > 0)
              unhealthy = TRUE;
          }
        } // for each data rectord

        if(unhealthy)
        {
          logInfo<<"  "<<iter->first.prnStr()<<" is unhealthy: disabled"<<Log::endl;
          iter = satellites.erase(iter);
        }
        else
          iter++;
      }

    // disable satellties with inconsistent data records
    // -------------------------------------------------
    for(auto iter=satellites.begin(); iter!=satellites.end();)
    {
      Bool unhealthy = FALSE;
      for(UInt it1=0; it1<iter->second.size(); it1++)
        for(UInt it2=it1+1; it2<iter->second.size(); it2++)
        {
          auto &dat1 = iter->second.at(it1);
          auto &dat2 = iter->second.at(it2);
          if((dat1.msgType == dat2.msgType) && (dat1.timeClock == dat2.timeClock))
          {
            // check if data is same
            Bool different = FALSE;
            for(UInt i=0; i<std::min(dat1.data.size(), UInt(3)); i++)
              for(UInt k=0; k<dat1.data.at(i).size(); k++)
                if(std::fabs(dat1.data.at(i)(k)-dat2.data.at(i)(k)) > 1e-9 * std::max(std::fabs(dat1.data.at(i)(k)), std::fabs(dat2.data.at(i)(k))))
                  different = TRUE;
            if(different)
            {
              unhealthy = TRUE;
              break;
            }
          }
        } // for each data rectord

      if(unhealthy)
      {
        logWarning<<"  "<<iter->first.prnStr()<<" with inconsistent data sets: disabled"<<Log::endl;
        iter = satellites.erase(iter);
      }
      else
        iter++;
    }

    // convert times to GPS time ephemeris
    // -----------------------------------
    for(auto &satellite : satellites)
      for(auto &navData : satellite.second)
      {
        if((navData.prn == GnssType::GPS)     ||
           (navData.prn == GnssType::GALILEO) ||
           (navData.prn == GnssType::QZSS)    ||
           (navData.prn == GnssType::IRNSS))
        {
          navData.timeEph = date2time(1980, 1, 6) + seconds2time(navData.data.at(2)(0));
          if(navData.msgType == "CNAV" || navData.msgType == "CNAV2")
            navData.timeEph = navData.timeClock;
        }
        else if(navData.prn == GnssType::GLONASS)
        {
          navData.timeEph = navData.timeClock = timeUTC2GPS(navData.timeClock);
        }
        else if(navData.prn == GnssType::SBAS)
        {
          navData.timeEph = navData.timeClock;
        }
        else if(navData.prn == GnssType::BDS)
        {
          navData.timeClock -= seconds2time(14);
          navData.timeEph    = date2time(2006, 1, 1, 0, 0, 14) + seconds2time(navData.data.at(2)(0));
        }

        // Adjust week
        while((navData.timeEph-navData.timeClock).mjd() > 3.5)
          navData.timeEph -= mjd2time(7);
        while((navData.timeEph-navData.timeClock).mjd() < -3.5)
          navData.timeEph += mjd2time(7);
      }

    // compute orbit
    // -------------
    if(!outNameOrbit.empty())
    {
      logStatus<<"writing "<<satellites.size()<<" orbit files <"<<outNameOrbit.appendBaseName(".{prn}")<<">"<<Log::endl;
      for(const auto &satellite : satellites)
        if(satellite.second.size())
        {
          OrbitArc arc;
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          {
            OrbitEpoch epoch;
            epoch.time = times.at(idEpoch);

            const NavData &navData = *std::min_element(satellite.second.begin(), satellite.second.end(),
                                      [&](auto &d1, auto &d2){return std::fabs((d1.timeEph-times.at(idEpoch)).seconds()) < std::fabs((d2.timeEph-times.at(idEpoch)).seconds());});
            const Double dt = (epoch.time - navData.timeEph).seconds();

            if(satellite.first == GnssType::GLONASS)
            {
              const Double GM  = 3.9860044e14;
              const Double C20 = -1082.63e-6;
              const Double a_e = 6378136;

              // reference epoch
              Rotary3d crf2trf = earthRotation->rotaryMatrix(navData.timeEph);
              epoch.time       = navData.timeEph;
              epoch.position   = crf2trf.inverseRotate(1e3*Vector3d(navData.data.at(0)(0), navData.data.at(1)(0), navData.data.at(2)(0)));
              epoch.velocity   = crf2trf.inverseRotate(1e3*Vector3d(navData.data.at(0)(1), navData.data.at(1)(1), navData.data.at(2)(1)))
                               + crossProduct(earthRotation->rotaryAxis(navData.timeEph), epoch.position);
              Vector3d sunMoonAcceleration = 1e3 * Vector3d(navData.data.at(0)(2), navData.data.at(1)(2), navData.data.at(2)(2));

              // Runge-Kutta-4 integration with around 60 second step size
              if(navData.timeEph != times.at(idEpoch))
              {
                UInt steps = static_cast<UInt>(std::max(std::round(std::fabs(dt/60.)), 1.));
                for(UInt i=0; i<steps; i++)
                  epoch = rungeKutta4(navData.timeEph+seconds2time((i+1)*dt/steps), epoch, sunMoonAcceleration, earthRotation, GM, C20, a_e);
                // NOTE: integration can lead to errors up to ~10 m after 15 minutes, unclear if that's the accuracy limit or there's an issue somewhere in the code
              }
            }
            else if(satellite.first == GnssType::SBAS)
            {
              // ref to  RTCA DO229
              epoch.position     = Vector3d(navData.data.at(0)(0), navData.data.at(1)(0), navData.data.at(2)(0)) * 1e3;
              epoch.velocity     = Vector3d(navData.data.at(0)(1), navData.data.at(1)(1), navData.data.at(2)(1)) * 1e3;
              epoch.acceleration = Vector3d(navData.data.at(0)(2), navData.data.at(1)(2), navData.data.at(2)(2)) * 1e3;

              epoch.position += epoch.velocity*dt + epoch.acceleration*dt*dt/2.0;
              epoch.velocity += epoch.acceleration*dt;

              Rotary3d crf2trf = earthRotation->rotaryMatrix(epoch.time);
              Vector3d omega   = earthRotation->rotaryAxis(epoch.time);
              epoch.position   = crf2trf.inverseRotate(epoch.position); // TRF -> CRF
              epoch.velocity   = crf2trf.inverseRotate(epoch.velocity) + crossProduct(omega, epoch.position); // TRF -> CRF
            }
            else // all systems using GPS-like ephemerides
            {
              const std::map<GnssType, Double> GM      {{GnssType::GPS,     3.986005e14},
                                                        {GnssType::GALILEO, 3.986004418e14},
                                                        {GnssType::BDS,     3.986004418e14},
                                                        {GnssType::QZSS,    3.986005e14},
                                                        {GnssType::IRNSS,   3.986005e14}};

              const std::map<GnssType, Double> omega_e {{GnssType::GPS,     7.2921151467e-5},
                                                        {GnssType::GALILEO, 7.2921151467e-5},
                                                        {GnssType::BDS,     7.292115e-5},
                                                        {GnssType::QZSS,    7.2921151467e-5},
                                                        {GnssType::IRNSS,   7.2921151467e-5}};

              // GPS CNAV, GPS CNAV2 require different algorithm then LNAV for example
              // BDS CNAV, CNAV2, CNAV3
              // QZSS CNAV, CNAV2, CNAV3
              // sqrt(A) in rx4 is actually not sqrt(A)
              // https://www.gps.gov/technical/icwg/IS-GPS-800J.pdf
              // https://www.gps.gov/technical/icwg/IS-GPS-705J.pdf
              // http://en.beidou.gov.cn/SYSTEMS/ICD/201806/P020180608519640359959.pdf
              // http://en.beidou.gov.cn/SYSTEMS/ICD/201806/P020180608518432765621.pdf
              // https://qzss.go.jp/en/technical/download/pdf/ps-is-qzss/is-qzss-pnt-005.pdf?t=1708078968591
              // A, n and r are then required to be computed different
              Double n = std::sqrt(GM.at(satellite.first & GnssType::SYSTEM) / std::pow(navData.a(), 3)) + navData.delta_n();
              if((navData.prn == GnssType::GPS || navData.prn == GnssType::BDS || navData.prn == GnssType::QZSS) &&
                 (navData.msgType == "CNAV" || navData.msgType == "CNV1" || navData.msgType == "CNV2" || navData.msgType == "CNV3"))
                n += 0.5 * navData.data.at(4)(1) * dt;

              const Double M = std::fmod(navData.m0() + n * dt, 2.*PI);
              Double E = M;
              for(UInt i=0; i<10; i++)
                E = M + navData.e() * std::sin(E);
              E = std::fmod(E, 2.*PI);

              // true anomaly
              const Double nu = std::atan2(std::sqrt(1. - std::pow(navData.e(), 2.)) * std::sin(E), std::cos(E) - navData.e());
              //argument of latitude
              const Double du = navData.cuc() * std::cos(2. * (nu + navData.omega())) + navData.cus() * std::sin(2. * (nu + navData.omega()));
              const Double u  = nu + navData.omega() + du;
              // radius in orbital plane
              const Double dr = navData.crc() * std::cos(2. * (nu + navData.omega())) + navData.crs() * std::sin(2. * (nu + navData.omega()));

              Double r  = navData.a() * (1. - navData.e() * std::cos(E)) + dr;
              // sqrt(A) in rx4 is actually not sqrt(A)
              // A, n and r are then required to be computed different
              if((navData.prn == GnssType::GPS || navData.prn == GnssType::BDS || navData.prn == GnssType::QZSS) &&
                 (navData.msgType == "CNAV" || navData.msgType == "CNV1" || navData.msgType == "CNV2" || navData.msgType == "CNV3"))
                r = (navData.a() + navData.data.at(0)(0) * dt) * (1. - navData.e() * std::cos(E)) + dr;

              // inclination
              const Double di = navData.cic() * std::cos(2. * (nu + navData.omega())) + navData.cis() * std::sin(2. * (nu + navData.omega()));
              const Double i  = navData.i0() + navData.iDot() * dt + di;

              // longitude of ascending node
              Double lambda = navData.omega0()  + (navData.omegaDot() - omega_e.at(satellite.first & GnssType::SYSTEM)) * dt
                            - omega_e.at(satellite.first & GnssType::SYSTEM) * navData.data.at(2)(0);
              // CNAV toc = TOE. navData.data.at(2)(0) would be TOP for CNAV and CNAV2 GPS
              // https://files.igs.org/pub/data/format/rinex_4.00.pdf
              if(((navData.prn == GnssType::GPS) || (navData.prn == GnssType::QZSS)) &&
                 (navData.msgType == "CNAV" || navData.msgType == "CNV2"))
              {
                Double gpsWeek = (navData.timeEph.mjdInt() - date2time(1980, 1, 6).mjdInt()) / 7;
                Double gpsSecond = ((navData.timeEph.mjd() - date2time(1980, 1, 6).mjd()) / 7 - gpsWeek) * 7 * 24 * 3600;
                lambda = navData.omega0()  + (navData.omegaDot() - omega_e.at(satellite.first & GnssType::SYSTEM)) * dt
                       - omega_e.at(satellite.first & GnssType::SYSTEM) * gpsSecond;
              }
              // http://en.beidou.gov.cn/SYSTEMS/ICD/201902/P020190227702348791891.pdf p37
              else if((navData.prn == GnssType::BDS) && (navData.prn.prn() <= 5 || (59 <= navData.prn.prn() && navData.prn.prn() <= 63)))
                lambda = navData.omega0()  + navData.omegaDot() * dt
                       - omega_e.at(satellite.first & GnssType::SYSTEM) * navData.data.at(2)(0);

              // Position
              const Double xOrb = r * std::cos(u);
              const Double yOrb = r * std::sin(u);
              epoch.position.x() = xOrb * std::cos(lambda) - yOrb * std::sin(lambda) * std::cos(i);
              epoch.position.y() = xOrb * std::sin(lambda) + yOrb * std::cos(lambda) * std::cos(i);
              epoch.position.z() = yOrb * std::sin(i);

              // http://en.beidou.gov.cn/SYSTEMS/ICD/201902/P020190227702348791891.pdf p 37
              if((navData.prn == GnssType::BDS) && (navData.prn.prn() <= 5 || (navData.prn.prn() >= 59 && navData.prn.prn() <= 63)))
               epoch.position = (rotaryZ(Angle(omega_e.at(satellite.first & GnssType::SYSTEM)*dt))*rotaryX(Angle(-5*DEG2RAD))).rotate(epoch.position);

              epoch.position = earthRotation->rotaryMatrix(epoch.time).inverseRotate(epoch.position); // TRF -> CRF
            }
            arc.push_back(epoch);
          }

          InstrumentFile::write(outNameOrbit.appendBaseName('.'+satellite.first.prnStr()), arc);
        }
    }

    if(!outNameClock.empty())
    {
      logStatus<<"writing "<<satellites.size()<<" clock files <"<<outNameClock.appendBaseName(".{prn}")<<">"<<Log::endl;
      for(const auto &satellite : satellites)
        if(satellite.second.size())
        {
          MiscValueArc arc;
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          {
            const NavData &navData = *std::min_element(satellite.second.begin(), satellite.second.end(),
                                                      [&](auto &d1, auto &d2){return std::fabs((d1.timeClock-times.at(idEpoch)).seconds()) < std::fabs((d2.timeClock-times.at(idEpoch)).seconds());});

            MiscValueEpoch epoch;
            epoch.time = times.at(idEpoch);
            const Double dt = (epoch.time - navData.timeClock).seconds();
            for(UInt n=0; n<((navData.prn == GnssType::GLONASS || navData.prn == GnssType::SBAS) ? 2 : 3); n++)
              epoch.value += navData.clock(n)*std::pow(dt, n);
            arc.push_back(epoch);
          }

          // write clock file
          InstrumentFile::write(outNameClock.appendBaseName("."+satellite.first.prnStr()), arc);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssRinexNavigation2OrbitClock::readRinex(const FileName &fileName, std::map<GnssType, std::vector<NavData>> &satellites)
{
  std::string line, label;

  try
  {
    auto getLine = [](InFile &file, std::string &line, std::string &label)
    {
      try
      {
        std::getline(file, line);
        if(line.back() == '\r')
          line.pop_back();
        line.resize(std::max(UInt(80), line.size()), ' ');
        label = line.substr(60, 20);
        return file.good();
      }
      catch(...)
      {
        line.clear();
        line.resize(80,' ');
        label = line.substr(60, 20);
        return FALSE;
      }
    };

    auto testLabel = [](const std::string &labelInLine, const std::string &label, Bool optional=TRUE)
    {
      if(labelInLine.find(label) != std::string::npos)
        return TRUE;
      if(optional)
        return FALSE;
      throw(Exception(std::string("In Line '")+labelInLine+"' label '"+label+"' expected\n"));
    };

    InFile file(fileName);

    getLine(file, line, label);
    testLabel(label, "RINEX VERSION / TYPE", FALSE);
    const Double rinexVersion = String::toDouble(line.substr(0, 9));
    if(rinexVersion < 2)
      logWarning<<"old RINEX version: "<<rinexVersion<<Log::endl;
    if((line.at(20) != 'N') && !(rinexVersion < 3 && line.at(20) == 'G'))
      throw(Exception("File must contain Navigation Data"));
    Char system = (rinexVersion) < 3 ? (line.at(20) == 'N' ? 'G' : 'R') : line.at(40);

    // read header
    // -----------
    for(;;)
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
      if(testLabel(label, "END OF HEADER"))
        break;
    }

    // read data
    // ---------
    if(std::getline(file, line))
      for(;;)
      {
        if(line.empty())
          break;
        line.resize(std::max(UInt(80), line.size()), ' ');

        // RINEX 4 record identifier
        std::string messageType;
        if(rinexVersion >= 4)
        {
          if(line.substr(0,5) != "> EPH") // only ephemeris records are of interest
          {
            getline(file, line);
            continue;
          }
          messageType = String::trim(line.substr(10,4));
          std::getline(file, line);
          line.resize(std::max(UInt(80), line.size()), ' ');
        }

        std::string prnStr = (rinexVersion < 3) ? system+line.substr(0,2): line.substr(0,3);
        if(prnStr.at(1) == ' ') prnStr.at(1) = '0';
        GnssType prn("***"+prnStr);

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
        Vector clock(3);
        for(UInt i=0; i<clock.size(); i++)
          clock(i) = String::toDouble(line.substr((rinexVersion < 3 ? 22 : 23)+i*19, 19));

        // data
        std::vector<Vector> data;
        while(std::getline(file, line))
        {
          if(line.empty() || (line.at(0) != ' '))
            break;
          line.resize(std::max(UInt(80), line.size()), ' ');
          Vector dataLine(4);
          for(UInt k=0; k<4; k++)
            dataLine(k) = String::toDouble(line.substr(((rinexVersion < 3) ? 3 : 4)+k*19, 19));
          data.push_back(dataLine);
        }

        satellites[prn].push_back(NavData(messageType, prn, time, clock, data));
      }
  }
  catch(std::exception &e)
  {
    logError<<"'"<<line<<"'"<<Log::endl;
    GROOPS_RETHROW(e)
  }
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
