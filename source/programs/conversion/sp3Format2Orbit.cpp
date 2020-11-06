/***********************************************/
/**
* @file sp3Format2Orbit.cpp
*
* @brief Read IGS orbits from SP3 format.
*
* @author Torsten Mayer-Guerr
* @date 2019-10-25
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read IGS orbits from \href{ftp://igs.org/pub/data/format/sp3c.txt}{SP3 format}
and write an \file{instrument file (ORBIT)}{instrument}.
The additional \config{outputfileClock} is an \file{instrument file (MISCVALUE)}{instrument}
and \config{outputfileCovariance} is an \file{instrument file (COVARIANCE3D)}{instrument}.

If \configClass{earthRotation}{earthRotationType} is provided the data are transformed
from terrestrial (TRF) to celestial reference frame (CRF).

See also \program{Orbit2Sp3Format}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Read IGS orbits from SP3 format.
* @ingroup programsConversionGroup */
class Sp3Format2Orbit
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Sp3Format2Orbit, SINGLEPROCESS, "read IGS orbits from SP3 format", Conversion, Orbit, Covariance, Instrument)
GROOPS_RENAMED_PROGRAM(Sp3file2Orbit, Sp3Format2Orbit, date2time(2020, 8, 4))
GROOPS_RENAMED_PROGRAM(Igs2Orbit,     Sp3Format2Orbit, date2time(2020, 8, 4))

/***********************************************/

void Sp3Format2Orbit::run(Config &config)
{
  try
  {
    FileName              fileNameOrbit, fileNameClock, fileNameCov;
    std::string           satId;
    EarthRotationPtr      earthRotation;
    std::vector<FileName> fileNamesIn;

    readConfig(config, "outputfileOrbit",      fileNameOrbit, Config::MUSTSET,  "", "");
    readConfig(config, "outputfileClock",      fileNameClock, Config::OPTIONAL, "", "");
    readConfig(config, "outputfileCovariance", fileNameCov,   Config::OPTIONAL, "", "3x3 epoch covariance");
    readConfig(config, "satelliteIdentifier",  satId,         Config::OPTIONAL, "", "e.g. L09 for GRACE A, empty: take first satellite");
    readConfig(config, "earthRotation",        earthRotation, Config::OPTIONAL, "", "rotation from TRF to CRF");
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

        // read header
        // -----------
        std::string line;           //                  0123456789|123456789|123456789|123456789|123456789|123456789
        std::getline(file, line);   // SP3 Line 01      #cV2018 12 24 13 51  0.00000000   14307 ORBIT ITRF  FIT CNES
        std::getline(file, line);   // SP3 Line 02      ## 2033 136260.00000000    60.00000000 58476 0.5770833333333
        std::getline(file, line);   // SP3 Line 03      +    1   L39  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        if(satId.empty())
          satId = line.substr(9, 3);
        std::getline(file, line);   // SP3 Line 04      +          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        std::getline(file, line);   // SP3 Line 05      +          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        std::getline(file, line);   // SP3 Line 06      +          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        std::getline(file, line);   // SP3 Line 07      +          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        std::getline(file, line);   // SP3 Line 08      ++
        std::getline(file, line);   // SP3 Line 09      ++
        std::getline(file, line);   // SP3 Line 10      ++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        std::getline(file, line);   // SP3 Line 11      ++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        std::getline(file, line);   // SP3 Line 12      ++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        std::getline(file, line);   // SP3 Line 13      %c L  cc TAI ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
        enum TimeSystem {GPS, UTC, TAI};
        TimeSystem timeSystem = GPS;
        if(line.substr(9, 3) == "GPS") timeSystem = GPS;
        else if(line.substr(9, 3) == "UTC") timeSystem = UTC;
        else if(line.substr(9, 3) == "TAI") timeSystem = TAI;
        else logWarning<<"Unkonwn time system ("<<line.substr(9, 3)<<"), assuming GPS time"<<Log::endl;
        std::getline(file, line);   // SP3 Line 14      %c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
        std::getline(file, line);   // SP3 Line 15      %f  1.2500000  1.025000000  0.00000000000  0.000000000000000
        std::getline(file, line);   // SP3 Line 16      %f  0.0000000  0.000000000  0.00000000000  0.000000000000000
        std::getline(file, line);   // SP3 Line 17      %i    0    0    0    0      0      0      0      0         0
        std::getline(file, line);   // SP3 Line 18      %i    0    0    0    0      0      0      0      0         0
        for(UInt i=0; i<4; i++)
          std::getline(file, line); // SP3 Lines 19-22  /* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        Time time;
        Bool positionRecord = FALSE;
        for(;;)
        {
          std::string line;
          std::getline(file, line);
          std::string lineID = line.substr(0,2);
          std::stringstream ss(line.substr(1));

          // Epoch
          // -----
          if(lineID=="* ")
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
          }

          // Position
          // --------
          if(lineID.at(0)=='P')
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

              OrbitEpoch epoch;
              epoch.time     = time;
              epoch.position = 1e3*Vector3d(x,y,z); // km -> m
              orbit.push_back(epoch);
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
          if(lineID=="EP" && positionRecord)
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
          if((lineID.at(0)=='V') && (line.substr(1,3) == satId) && positionRecord)
          {
            Double x = String::toDouble(line.substr(4, 14));
            Double y = String::toDouble(line.substr(18, 14));
            Double z = String::toDouble(line.substr(32, 14));
            orbit.at(orbit.size()-1).velocity = 0.1*Vector3d(x,y,z);  // dm/s -> m/s
          }

          // end of file
          // -----------
          if(lineID=="EO")
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
      logTimerStart;
      for(UInt i=0; i<orbit.size(); i++)
      {
        logTimerLoop(i, orbit.size());
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
      }
      logTimerLoopEnd(orbit.size());
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
