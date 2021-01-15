/***********************************************/
/**
* @file berneseKinematic2Orbit.cpp
*
* @brief Read kinematic orbits in Bernese format.
*
* @author Torsten Mayer-Guerr
* @date 2009-04-03
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read kinematic orbits in Bernese format.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Read kinematic orbits in Bernese format.
* @ingroup programsConversionGroup */
class BerneseKinematic2Orbit
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(BerneseKinematic2Orbit, SINGLEPROCESS, "read kinematic orbits in Bernese format", Conversion, Orbit, Covariance, Instrument)
GROOPS_RENAMED_PROGRAM(Kinematic2Orbit, BerneseKinematic2Orbit, date2time(2020, 8, 4))

/***********************************************/

void BerneseKinematic2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName orbitName, covarianceName;
    std::vector<FileName> fileNamesIn;
    EarthRotationPtr earthRotation;

    readConfig(config, "outputfileOrbit",      orbitName,      Config::MUSTSET,  "", "");
    readConfig(config, "outputfileCovariance", covarianceName, Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",        earthRotation,  Config::OPTIONAL, "", "from TRF to CRF");
    readConfig(config, "inputfile",            fileNamesIn,    Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // ==============================

    OrbitArc        orbit;
    Covariance3dArc covariance;
    for(const auto &fileName : fileNamesIn)
    {
      try
      {
        logStatus<<"read file <"<<fileName<<">"<<Log::endl;
        InFile file(fileName);

        // header
        Double rms = 1.0;
        std::string line;
        for(UInt i=0; i<6; i++)
        {
          std::getline(file, line);
          if(file.eof())
            break;
          if(line.find("RMS") != std::string::npos)
            rms = String::toDouble(line.substr(97, 10));
        }

        while(std::getline(file, line))
        {
          if(line.empty())
            break;
          std::stringstream ss(line);
          ss.exceptions(std::ios::badbit|std::ios::failbit);

          std::string name, code;
          UInt   week;
          Double seconds;
          Char   flag;

          OrbitEpoch orbitEpoch;
          ss>>name>>code>>week>>seconds>>orbitEpoch.position.x()>>orbitEpoch.position.y()>>orbitEpoch.position.z()>>flag;
          orbitEpoch.time = date2time(1980, 1, 6) + mjd2time(week*7.) + seconds2time(seconds);

          Covariance3dEpoch covEpoch;
          ss>>covEpoch.covariance.xx()>>covEpoch.covariance.yy()>>covEpoch.covariance.zz();
          ss>>covEpoch.covariance.xy()>>covEpoch.covariance.xz()>>covEpoch.covariance.yz();
          covEpoch.time        = orbitEpoch.time;
          covEpoch.covariance *= (rms*rms);

          if(flag != 'K') //&&(flag!='G')) //&&(flag!='S')
            continue;

          orbit.push_back(orbitEpoch);
          covariance.push_back(covEpoch);
        }
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    } // for(idFile)

    // ==============================

    // Rotation TRF -> CRF
    // -------------------
    if(earthRotation)
    {
      logStatus<<"rotation from TRF to CRF"<<Log::endl;
      Single::forEach(orbit.size(), [&](UInt i)
      {
        const Rotary3d rotation = inverse(earthRotation->rotaryMatrix(orbit.at(i).time));
        orbit.at(i).position = rotation.rotate(orbit.at(i).position);
        if(orbit.at(i).velocity.r() > 0)
          orbit.at(i).velocity = rotation.rotate(orbit.at(i).velocity) + crossProduct(earthRotation->rotaryAxis(orbit.at(i).time), orbit.at(i).position);
        if(covariance.size())
          covariance.at(i).covariance = rotation.rotate(covariance.at(i).covariance);
      });
    }

    // ==============================

    // Write data
    // ----------
    if(!orbitName.empty())
    {
      logStatus<<"write orbit data to file <"<<orbitName<<">"<<Log::endl;
      InstrumentFile::write(orbitName, orbit);
      Arc::printStatistics(orbit);
    }
    if(!covarianceName.empty())
    {
      logStatus<<"write covariance data to file <"<<covarianceName<<">"<<Log::endl;
      InstrumentFile::write(covarianceName, covariance);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
