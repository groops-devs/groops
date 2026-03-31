/***********************************************/
/**
* @file tle2Orbit.cpp
*
* @brief Orbit from Two Line Elements (TLE/3LE).
*
* @author Torsten Mayer-Guerr
* @date 2023-04-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the \configFile{outputfileOrbit}{instrument}
from two-line elements (TLE/3LE)
as can be found at e.g. \url{http://celestrak.org/NORAD/elements/}.
The first satellite in the input file that matches the wildcard of \config{satelliteName} is used.
If more records with exactly the same name are found, the one with the closest reference epoch
is used for each point in the \configClass{timeSeries}{timeSeriesType}.

The program uses the Simplified General Perturbation (SGP) model. More information can
be found in the Revisiting Spacetrack Report 3 by Vallado et al. 2006.
)";

/***********************************************/

#include "programs/program.h"
#include "external/sgp4/SGP4.h"
#include "base/string.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Orbit from Two Line Elements (TLE/3LE).
* @ingroup programsGroup */
class Tle2Orbit
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Tle2Orbit, SINGLEPROCESS, "orbit from Two Line Elements (TLE/3LE)", Conversion, Instrument)

/***********************************************/

void Tle2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName         fileNameOrbit, fileNameTLE;
    std::string      wildcardName;
    TimeSeriesPtr    timeSeries;
    EarthRotationPtr earthRotation;

    readConfig(config, "outputfileOrbit", fileNameOrbit, Config::MUSTSET, "",  "");
    readConfig(config, "inputfileTLE",    fileNameTLE,   Config::MUSTSET, "",  "two line elements (TLE/3LE)");
    readConfig(config, "satelliteName",   wildcardName,  Config::DEFAULT, "*", "first name of wildcard match is used");
    readConfig(config, "timeSeries",      timeSeries,    Config::MUSTSET, "",  "output orbit at these times");
    readConfig(config, "earthRotation",   earthRotation, Config::MUSTSET, "",  "rotation to CRF");
    if(isCreateSchema(config)) return;

    logStatus<<"read TLE <"<<fileNameTLE<<">"<<Log::endl;
    std::vector<std::pair<Time, ElsetRec>> records;
    InFile file(fileNameTLE);
    const std::regex pattern = String::wildcard2regex(wildcardName);
    std::string line, satName;
    while(std::getline(file, line))
    {
      line = String::trim(line);
      std::string line1, line2;
      std::getline(file, line1);
      std::getline(file, line2);
      if((satName.empty() && std::regex_match(line, pattern)) || (line == satName))
      {
        satName = line;
        if(!(line1.size() && line2.size() && (line1.at(0) == '1') && (line2.at(0) == '2')))
          throw(Exception("parser error"));

        // first line
        const Double xpdotp = 24*60/(2.*PI);
        const Int yy  = String::toInt(line1.substr(18, 2)); // two digit year
        Time satEpoch = date2time(yy + ((yy > 56) ? 1900 : 2000), 1, 1) + mjd2time(String::toDouble(line1.substr(20, 12))-1);
        ElsetRec rec;
        rec.ndot     = String::toDouble(line1.substr(33, 10))/ (xpdotp*1440.);
        rec.nddot    = String::toDouble(line1.at(44)+"."s+line1.substr(45, 5)) * std::pow(10, String::toInt(line1.substr(50, 2)))/(xpdotp*1440.*1440.);
        rec.bstar    = String::toDouble(line1.at(53)+"."s+line1.substr(54, 5)) * std::pow(10, String::toInt(line1.substr(59, 2)));
        // second line
        rec.inclo    = String::toDouble(line2.substr( 8, 8)) * DEG2RAD;
        rec.nodeo    = String::toDouble(line2.substr(17, 8)) * DEG2RAD;
        rec.ecco     = String::toDouble("."s+line2.substr(26, 7));
        rec.argpo    = String::toDouble(line2.substr(34, 8)) * DEG2RAD;
        rec.mo       = String::toDouble(line2.substr(43, 8)) * DEG2RAD;
        rec.no_kozai = String::toDouble(line2.substr(52, 11))/xpdotp;

        rec.whichconst  = 2; // wgs72;
        rec.jdsatepoch  = satEpoch.mjdInt() + 2400000.5;
        rec.jdsatepochF = satEpoch.mjdMod();

        records.push_back(std::make_pair(satEpoch, rec));
        logInfo<<" used satellite: '"<<line<<"' "<<satEpoch.dateTimeStr()<<Log::endl;
      }
    }

    if(!records.size())
      throw(Exception("no satellite found"));

    logStatus<<"integrate orbit"<<Log::endl;
    const std::vector<Time> times = timeSeries->times();
    OrbitArc arc;

    auto iterRecord = records.begin();
    Single::forEach(times.size(), [&](UInt i)
    {
      auto iterRecordNew = std::min_element(records.begin(), records.end(),
                                            [&](auto &r1, auto &r2) {return std::fabs((times.at(i)-r1.first).mjd()) < std::fabs((times.at(i)-r2.first).mjd());});
      if(iterRecordNew != iterRecord)
        sgp4init('a', iterRecordNew->second);
      iterRecord = iterRecordNew;

      Double p[3], v[3];
      iterRecord->second.error = 0;
      sgp4(iterRecord->second, (timeGPS2UTC(times.at(i))-iterRecord->first).seconds()/60., p, v);
      if(iterRecord->second.error)
        logWarning<<times.at(i).dateTimeStr()<<" integration error"<<Log::endl;

      OrbitEpoch epoch;
      epoch.time     = times.at(i);
      epoch.position = Vector3d(1e3*p[0], 1e3*p[1], 1e3*p[2]);
      epoch.velocity = Vector3d(1e3*v[0], 1e3*v[1], 1e3*v[2]);

      // get actual precession and nutation
      Double xp, yp, sp, deltaUT, LOD, X, Y, S;
      earthRotation->earthOrientationParameter(times.at(i), xp, yp, sp, deltaUT, LOD, X, Y, S);
      const Double   r2  = X*X + Y*Y;
      const Double   E   = (r2 != 0.) ? std::atan2(Y, X) : 0.;
      const Double   D   = std::atan(std::sqrt(r2/(1-r2)));

      // rotate back precession
      const Double JC     = timeGPS2JC(times.at(i));
      const Double zetaA  = (2306.2181*JC + 0.30188*JC*JC + 0.017998*JC*JC*JC)* DEG2RAD/3600;
      const Double zA     = (2306.2181*JC + 1.09468*JC*JC + 0.018203*JC*JC*JC)* DEG2RAD/3600;
      const Rotary3d rot = rotaryZ(Angle(E)) * rotaryY(Angle(-D)) * rotaryZ(Angle(-E)) * rotaryZ(Angle(zetaA+zA));
      epoch.position     = rot.rotate(epoch.position);
      epoch.velocity     = rot.rotate(epoch.velocity);

      arc.push_back(epoch);
    });

    logStatus<<"write orbit <"<<fileNameOrbit<<">"<<Log::endl;
    InstrumentFile::write(fileNameOrbit, arc);
    Arc::printStatistics(arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
