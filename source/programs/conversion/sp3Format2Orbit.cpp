/***********************************************/
/**
* @file sp3Format2Orbit.cpp
*
* @brief Read IGS orbits from SP3 format.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2019-10-25
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read IGS orbits from \href{https://files.igs.org/pub/data/format/sp3d.pdf}{SP3 format}
and write an \file{instrument file (ORBIT)}{instrument}.
The additional \config{outputfileClock} is an \file{instrument file (MISCVALUE)}{instrument}
and \config{outputfileCovariance} is an \file{instrument file (COVARIANCE3D)}{instrument}.

If \configClass{earthRotation}{earthRotationType} is provided the data are transformed
from terrestrial (TRF) to celestial reference frame (CRF).
Since SP3 orbits often use the center of Earth as a reference, a correction from center
of Earth to center of mass can be applied to the orbits by providing \configClass{gravityfield}{gravityfieldType} (e.g. ocean tides).

See also \program{Orbit2Sp3Format}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Read IGS orbits from SP3 format.
* @ingroup programsConversionGroup */
class Sp3Format2Orbit
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Sp3Format2Orbit, SINGLEPROCESS, "read IGS orbits from SP3 format", Conversion, Orbit, Covariance, Instrument)
GROOPS_RENAMED_PROGRAM(Sp3file2Orbit, Sp3Format2Orbit, date2time(2020, 8, 4))
GROOPS_RENAMED_PROGRAM(Igs2Orbit,     Sp3Format2Orbit, date2time(2020, 8, 4))

/***********************************************/

void Sp3Format2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOrbit, fileNameClock, fileNameCov;
    std::string           satId;
    EarthRotationPtr      earthRotation;
    GravityfieldPtr       gravityfield;
    std::vector<FileName> fileNamesIn;

    readConfig(config, "outputfileOrbit",      fileNameOrbit, Config::MUSTSET,  "", "");
    readConfig(config, "outputfileClock",      fileNameClock, Config::OPTIONAL, "", "");
    readConfig(config, "outputfileCovariance", fileNameCov,   Config::OPTIONAL, "", "3x3 epoch covariance");
    readConfig(config, "satelliteIdentifier",  satId,         Config::OPTIONAL, "", "e.g. L09 for GRACE A, empty: take first satellite");
    readConfig(config, "earthRotation",        earthRotation, Config::OPTIONAL, "file", "rotation from TRF to CRF");
    readConfig(config, "gravityfield",         gravityfield,  Config::DEFAULT,  R"([{"tides": {"tides": {"doodsonHarmonicTide": {"minDegree":1, "maxDegree":1}}}}])", "degree 1 fluid mantle for CM2CE correction (SP3 orbits should be in center of Earth)");
    readConfig(config, "inputfile",            fileNamesIn,   Config::MUSTSET,  "", "orbits in SP3 format");
    if(isCreateSchema(config)) return;

    // ==============================

    OrbitArc        orbit;
    MiscValueArc    clock;
    Covariance3dArc cov;
    for(FileName &filenameIn : fileNamesIn)
    {
      try
      {
        logStatus<<"read file <"<<filenameIn<<">"<<Log::endl;
        InFile file(filenameIn);
        std::string line;
        enum TimeSystem {GPS, UTC, TAI};
        TimeSystem timeSystem = GPS;
        Time time;
        Vector3d cm2ceCorrection;
        Bool positionRecord = FALSE;
        while(std::getline(file, line))
        {
          // Header
          // ------
          if(String::startsWith(line, "#") ||   // first 2 lines
             String::startsWith(line, "/*") ||  // comment lines
             String::startsWith(line, "%f") ||  // floating point base base numbers
             String::startsWith(line, "%i"))    // additional parameters
            continue;

          if(String::startsWith(line, "+"))     // satellite list and orbit accuracy lines
          {
            if(satId.empty() && String::toInt(line.substr(3, 3)) > 0)
              satId = line.substr(9, 3);
            continue;
          }

          if(String::startsWith(line, "%c"))    // file type and time system definition lines
          {
            if(line.substr(9, 3) == "GPS") timeSystem = GPS;
            else if(line.substr(9, 3) == "UTC") timeSystem = UTC;
            else if(line.substr(9, 3) == "TAI") timeSystem = TAI;
            else logWarning<<"Unknown time system ("<<line.substr(9, 3)<<"), assuming GPS time"<<Log::endl;
            std::getline(file, line); // skip second %c line
            continue;
          }

          // Epoch
          // -----
          if(String::startsWith(line, "* "))
          {
            UInt   year  = String::toInt(line.substr(3, 4));
            UInt   month = String::toInt(line.substr(8, 2));
            UInt   day   = String::toInt(line.substr(11, 2));
            UInt   hour  = String::toInt(line.substr(14, 2));
            UInt   min   = String::toInt(line.substr(17, 2));
            Double sec   = String::toDouble(line.substr(20, 11));
            time = date2time(year, month, day, hour, min, sec);
            if(timeSystem == UTC)
              time = timeUTC2GPS(time);
            else if(timeSystem == TAI)
              time -= seconds2time(DELTA_TAI_GPS);
            positionRecord = FALSE;

            const SphericalHarmonics harmonics = gravityfield->sphericalHarmonics(time, 1, 1);
            const Vector coeff = harmonics.x(); // [c00, c10, c11, s11]
            cm2ceCorrection = std::sqrt(3.) * harmonics.R() * Vector3d(coeff(2), coeff(3), coeff(1));
          }

          // Position
          // --------
          if(String::startsWith(line, "P"))
          {
            if(line.substr(1,3) != satId)
              positionRecord = FALSE;
            else
            {
              positionRecord = TRUE;
              Double x = String::toDouble(line.substr(4, 14));
              Double y = String::toDouble(line.substr(18, 14));
              Double z = String::toDouble(line.substr(32, 14));
              Double c = String::toDouble(line.substr(46, 14));

              if(x != 0. && y != 0. && z != 0.)
              {
                OrbitEpoch epoch;
                epoch.time     = time;
                epoch.position = 1e3*Vector3d(x,y,z) - cm2ceCorrection; // km -> m
                orbit.push_back(epoch);
              }
              if(c < 999999)
              {
                MiscValueEpoch epoch;
                epoch.time  = time;
                epoch.value = 1e-6*c; // microsecond -> second
                clock.push_back(epoch);
              }
            }
          }

          // Position covariance
          // -------------------
          if(String::startsWith(line, "EP") && positionRecord)
          {
            Double xx = String::toDouble(line.substr(4, 4));
            Double yy = String::toDouble(line.substr(9, 4));
            Double zz = String::toDouble(line.substr(14, 4));
            Double xy = String::toDouble(line.substr(27, 8));
            Double xz = String::toDouble(line.substr(36, 8));
            Double yz = String::toDouble(line.substr(54, 8));
            Covariance3dEpoch epochCov;
            epochCov.time = time;
            // mm -> m, correlation [1e-7] -> covariance
            epochCov.setData(Vector({std::pow(1e-3*xx,2), std::pow(1e-3*yy,2), std::pow(1e-3*zz,2), 1e-13*xy*xx*yy, 1e-13*xz*xx*zz, 1e-13*yz*yy*zz}));
            cov.push_back(epochCov);
          }

          // Velocity
          // --------
          if(String::startsWith(line, "V") && (line.substr(1,3) == satId) && positionRecord)
          {
            Double x = String::toDouble(line.substr(4, 14));
            Double y = String::toDouble(line.substr(18, 14));
            Double z = String::toDouble(line.substr(32, 14));
            if(x != 0. && y != 0. && z != 0.)
              orbit.at(orbit.size()-1).velocity = 0.1*Vector3d(x,y,z);  // dm/s -> m/s
          }

          // end of file
          // -----------
          if(String::startsWith(line, "EOF"))
            break;
        } // for(;;)
      }
      catch(std::exception &e)
      {
        logWarning<<std::endl<<e.what()<<" continue..."<<Log::endl;
        break;
      }
    } // for(inputFiles)

    if(orbit.size() == 0)
      throw(Exception("empty arc"));

    // ==============================

    // Rotation TRF -> CRF
    // -------------------
    if(earthRotation)
    {
      logStatus<<"rotation from TRF to CRF"<<Log::endl;
      UInt idxCov = 0;
      Single::forEach(orbit.size(), [&](UInt i)
      {
        const Rotary3d rotation = inverse(earthRotation->rotaryMatrix(orbit.at(i).time));
        const Vector3d omega    = earthRotation->rotaryAxis(orbit.at(i).time);
        orbit.at(i).position = rotation.rotate(orbit.at(i).position);
        if(orbit.at(i).velocity.r() > 0)
          orbit.at(i).velocity = rotation.rotate(orbit.at(i).velocity) + crossProduct(omega, orbit.at(i).position);
        if(cov.size() && !fileNameCov.empty())
        {
          while((idxCov < cov.size()) && (cov.at(idxCov).time < orbit.at(i).time))
            idxCov++;
          if((idxCov < cov.size()) && (cov.at(idxCov).time == orbit.at(i).time))
            cov.at(idxCov).covariance = rotation.rotate(cov.at(idxCov).covariance);
        }
      });
    }

    // ==============================

    // write results
    // -------------
    if(!fileNameOrbit.empty() && orbit.size())
    {
      logStatus<<"write orbit data to file <"<<fileNameOrbit<<">"<<Log::endl;
      InstrumentFile::write(fileNameOrbit, orbit);
      Arc::printStatistics(orbit);
    }
    if(!fileNameClock.empty() && clock.size())
    {
      logStatus<<"write clock data to file <"<<fileNameClock<<">"<<Log::endl;
      InstrumentFile::write(fileNameClock, clock);
    }
    if(!fileNameCov.empty() && cov.size())
    {
      logStatus<<"write covariance data to file <"<<fileNameCov<<">"<<Log::endl;
      InstrumentFile::write(fileNameCov, cov);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
