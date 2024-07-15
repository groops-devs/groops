/***********************************************/
/**
* @file orbit2Cpf.cpp
*
* @brief write orbit positions to CPF.
*
* @author Cornelia Tieber-Hubmann
* @date 2023-03-17
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Writes groops orbits to \href{https://ilrs.gsfc.nasa.gov/data_and_products/formats/cpf.html}{CPF file}.

The coordinate system used in the CPF format is usually presented in ITRF.
The required time format for the input orbit file is GPS.
The time format of the output CPF file is given in UTC.

See also \program{Cpf2Orbit}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "base/time.h"
#include "files/filePlatform.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief write orbit positions to CPF.
* @ingroup programsConversionGroup */
class Orbit2Cpf
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Orbit2Cpf, SINGLEPROCESS, "write orbit positions to CPF", Conversion, Orbit, Instrument)


/***********************************************/

void Orbit2Cpf::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName          fileNameOut, fileNameIn;
    EarthRotationPtr  earthRotation;
    FileName          fileNameSatelliteInfo;
    UInt              subDailyEphemerisSequenceNumber;
    UInt              targetClass;

    readConfig(config, "outputfile",             fileNameOut,                     Config::MUSTSET,  "{target}_cpf_{fileDate}_{ephVerNum}{ver}.tug", "");
    readConfig(config, "inputfileOrbit",         fileNameIn,                      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileSatelliteInfo", fileNameSatelliteInfo,           Config::MUSTSET,  "{groopsDataDir}/slr/satellite/satelliteInfo/satelliteInfo.{satellite}.xml", "Platform File");
    readConfig(config, "earthRotation",          earthRotation,                   Config::MUSTSET,  "", "");
    readConfig(config, "versionNumber",          subDailyEphemerisSequenceNumber, Config::DEFAULT,  "01", "Version number of production day with zero leading fill, e.g. 01");
    readConfig(config, "targetClass",            targetClass,                     Config::DEFAULT,  "1",  "set 1 for passive retroreflector, set 0 for no retroreflector (includes debris)");
    if(isCreateSchema(config)) return;

    // parameters from platform file
    std::string sic("-1"), cospar("-1"), norad("-1"), targetName("-1");
    Platform platform;
    readFilePlatform(fileNameSatelliteInfo, platform);
    auto satId = platform.findEquipment<PlatformSatelliteIdentifier>(Time());
    if(satId)
    {
      if(!satId->sic.empty())    sic        = satId->sic;
      if(!satId->cospar.empty()) cospar     = satId->cospar;
      if(!satId->norad.empty())  norad      = satId->norad;
      if(!satId->name.empty())   targetName = satId->name;
    }

    // Set parameters for header
    Int compatibilityTIVs = 1; // integrable, geocentric ephemeris
    Int referenceFrame    = 0; // 0=geocentric true body-fixed (default), 1=geocentric space-fixed (i.e., Inertial) (True-of-Date), 2=geocentric space-fixed (Mean-of-Date J2000)
    Int rotAngleType      = 0; // 0=Not Applicable, 1=Lunar Euler angles, 2=North pole Right Ascension and Declination, and angle to prime meridian
    Int centerOfMassCorr  = 0; // 0=None applied. Prediction is for center of mass of target, 1=Applied. Prediction is for retro-reflector array
    Int targetLocation    = 1; // 1=Earth orbit, 0&2-10 see CPF pdf

    // Starting date of orbit
    // see Consolidated Laser Ranging Prediction Format Version 2.00 - Start dates and times
    logStatus<<"read orbit file <"<<fileNameIn<<">"<<Log::endl;
    OrbitArc orbit = InstrumentFile::read(fileNameIn);
    const Time timeStart = timeGPS2UTC(orbit.front().time);
    const Time timeEnd   = timeGPS2UTC(orbit.back().time);
    logInfo<<"  timeStart (UTC): "<<timeStart.dateTimeStr()<<Log::endl;

    logStatus<<"write file <"<<fileNameOut<<">"<<Log::endl;
    OutFile file(fileNameOut);
    // Write the header information to output file
    file<<"H1 CPF 2 tug "<<System::now()%"%y %m %d %H %O"s<<" "<<subDailyEphemerisSequenceNumber<<" "<<targetName<<" "<<"NotesOptional"<<std::endl;
	  file<<"H2 "<<cospar<<" "<<sic<<" "<<norad<<" "<<timeStart%"%y %m %d %H %M %S "s<<timeEnd%"%y %m %d %H %M %S "s
               <<medianSampling(orbit.times()).seconds()<<" "<<compatibilityTIVs<<" "<<targetClass<<" "<<referenceFrame<<" "<<rotAngleType<<" "
               <<centerOfMassCorr<<" "<<targetLocation<<std::endl;
	  file<<"H9"<<std::endl;

    // Write position information in ECEF reference frame to output file
    Single::forEach(orbit.size(), [&](UInt i)
    {
      const Vector3d positionCPF = earthRotation->rotaryMatrix(orbit.at(i).time).rotate(orbit.at(i).position);
      const Time     timeUTC     = timeGPS2UTC(orbit.at(i).time);
      file<<"10 0 "<<timeUTC.mjdInt()<<" "<<timeUTC.mjdMod()*86400.<<" 0 "<<positionCPF.x()%"%14.9f "s<<positionCPF.y()%"%14.9f "s<<positionCPF.z()%"%14.9f"s<<std::endl;
    });

    file<<"99"<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
