/***********************************************/
/**
* @file cpf2Orbit.cpp
*
* @brief Read orbit data from SLR prediction (CPF) format
*
* @author Barbara Suesser-Rechberger
* @date 2022-11-24
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts \href{https://ilrs.gsfc.nasa.gov/data_and_products/formats/cpf.html}{CPF file}
and writes an \file{instrument file (ORBIT)}{instrument}.

The time format of the CPF file is UTC.
The coordinate system used in the CPF format is usually represented in TRF.
If \configClass{earthRotation}{earthRotationType} is provided the data are transformed
from terrestrial (TRF) to celestial reference frame (CRF).

See also \program{Orbit2Cpf}
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "classes/earthRotation/earthRotation.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read orbit from SLR prediction (CPF) format
* @ingroup programsConversionGroup */
class Cpf2Orbit
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Cpf2Orbit, SINGLEPROCESS, "Read orbit from SLR prediction (CPF) format", Conversion, Slr, Instrument)

/***********************************************/

void Cpf2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName         fileNameOrbit;
    FileName         fileNameIn;
    EarthRotationPtr earthRotation;

    readConfig(config, "outputfileOrbit", fileNameOrbit, Config::MUSTSET, "output/orbit_{satellite}_{loopNumberEphVersion:%03i}{loopNumberVersion:%02i}_{loopSource}_{loopTime:%D}.dat", "");
    readConfig(config, "earthRotation",   earthRotation, Config::MUSTSET, "", "");
    readConfig(config, "inputfile",       fileNameIn,    Config::MUSTSET, "", "SLR CPF file");
    if(isCreateSchema(config)) return;

    logStatus<<"read file <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);
    std::string line;

    std::string satelliteName;
    std::map<Time, OrbitEpoch> orbEpoch;
    Time time;

    // Get a reference to OrbitEpoch object at specified time. OrbitEpoch object
    // will be default created (all NAN) if it does not already exist.
    auto orbitEpochGetOrInsert = [&orbEpoch](const Time &time)-> OrbitEpoch&
    {
      if(orbEpoch.count(time) == 0)
      {
        orbEpoch[time] = {};
        orbEpoch[time].time = time;
        orbEpoch[time].position = Vector3d(NAN,NAN,NAN);
        orbEpoch[time].velocity = Vector3d(NAN,NAN,NAN);
        orbEpoch[time].acceleration = Vector3d(NAN,NAN,NAN);
      }
      return orbEpoch[time];
    };

    while(std::getline(file, line))
    {
      if(line.empty())
        continue;
      std::string type = String::upperCase(line.substr(0,2));
      std::stringstream ss(line.substr(2));

      if(type == "H1") // basic information 1
      {
        Int version, year, month, day, hour, ephemeridenSubSequenceNumber;
        std::string cpf, source, ephemeridenSequenceNumber;

        ss>>cpf>>version>>source>>year>>month>>day>>hour>>ephemeridenSequenceNumber>>ephemeridenSubSequenceNumber>>satelliteName;

      } else if(type == "H2") // basic information 2
      {
        Int cosparId, sic, noradId, startYear, startMonth, startDay, startHour, startMinute, startSecond, endYear, endMonth, endDay, endHour, endMinute, endSecond, time, comp, targetType, refFrame, comCorr;

        ss>>cosparId>>sic>>noradId>>startYear>>startMonth>>startDay>>startHour>>startMinute>>startSecond>>endYear>>endMonth>>endDay>>endHour>>endMinute>>endSecond>>time>>comp>>targetType>>refFrame>>comCorr;

        if(refFrame > 0)
        {
          // Other reference frame than geocentric true body-fixed is not considered at the moment.
          logStatus<<"Other reference frame than geocentric true body-fixed is not considered at the moment: satellite <"<<satelliteName<<">"<<Log::endl;
          continue;
        }

        if(comCorr > 0)
        {
          // COM correction applied, Prediction is for retro-reflector array.
          logStatus<<"Center of mass correction applied. Prediction is for retroreflector array. Data set not considered: satellite: <"<<satelliteName<<">"<<Log::endl;
        }

      } if(type == "10") // position record
      {
        Int directionFlag, mjdDay,leapSecondFlag;
        Double sec, x, y, z;

        ss>>directionFlag>>mjdDay>>sec>>leapSecondFlag>>x>>y>>z;

        if(directionFlag > 0)
        {
          // direction flag 1 or 2 is not considered in the moment
          logStatus<<"Direction flag 1 or 2 is not considered in the moment: satellite <"<<satelliteName<<">"<<Log::endl;
          continue;
        }

        // SLR predictions are given in UTC
        time = timeUTC2GPS(mjd2time(mjdDay) + seconds2time(sec));
        OrbitEpoch& orbitEpoch = orbitEpochGetOrInsert(time);
        // Set the position
        orbitEpoch.position = Vector3d(x,y,z);

      }
      else if(type == "20") // velocity record
      {
        Int directionFlag;
        Double velX, velY, velZ;

        ss>>directionFlag>>velX>>velY>>velZ;

        if(directionFlag > 0)
        {
          // direction flag 1 or 2 is not considered in the moment
          logStatus<<"Direction flag 1 or 2 is not considered: satellite <"<<satelliteName<<">"<<Log::endl;
          continue;
        }

        OrbitEpoch& orbitEpoch = orbitEpochGetOrInsert(time);
        // Set velocity at correct time
        orbitEpoch.velocity = Vector3d(velX,velY,velZ);

      } else if(type == "30")
      {
        // Corrections
        logStatus<<"Record 30: Corrections not need for near earth satellites" <<Log::endl;
      } else if(type == "99")
      {
        // File end
      }
    } // while

    OrbitArc orbArc;

    for(const auto &orb : orbEpoch)
      orbArc.push_back(orb.second);

    if(earthRotation)
    {
      // Rotation from TRF to CRF
      logStatus<<"rotation from TRF to CRF"<<Log::endl;
      Single::forEach(orbArc.size(), [&](UInt i)
      {
      const Rotary3d rotation = inverse(earthRotation->rotaryMatrix(orbArc.at(i).time));
      orbArc.at(i).position = rotation.rotate(orbArc.at(i).position);
      if(orbArc.at(i).velocity.r() > 0)
          orbArc.at(i).velocity = rotation.rotate(orbArc.at(i).velocity) + crossProduct(earthRotation->rotaryAxis(orbArc.at(i).time), orbArc.at(i).position);
      });
    }

    // ==============================
    // write results
    // -------------
    if(orbArc.size())
    {
      logStatus<<"write orbit data to file <"<<fileNameOrbit<<">"<<Log::endl;
      InstrumentFile::write(fileNameOrbit, orbArc);
      Arc::printStatistics(orbArc);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
