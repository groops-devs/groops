/***********************************************/
/**
* @file groopsAscii2Orbit.cpp
*
* @brief Read Orbits given in groops kinematic orbit ASCII format.
*
* @author Norbert Zehentner
* @date 2015-01-23
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read Orbits given in groops kinematic orbit ASCII format with covariance information.

See also \program{Orbit2GroopsAscii}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Read Orbits given in groops kinematic orbit ASCII format.
* @ingroup programsConversionGroup */
class GroopsAscii2Orbit
{
  void readFile(const FileName &fileName, OrbitArc &arcOrb, Covariance3dArc &arcCov);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GroopsAscii2Orbit, SINGLEPROCESS, "read Orbits given in groosp kinematic ASCII format", Conversion, Orbit, Covariance, Instrument)
GROOPS_RENAMED_PROGRAM(AsciiKinematic2OrbitCovariance, GroopsAscii2Orbit, date2time(2020, 6, 14))

/***********************************************/

void GroopsAscii2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOrbit, fileNameCovariance;
    std::vector<FileName> fileNamesIn;
    EarthRotationPtr      earthRotation;

    readConfig(config, "outputfileOrbit",      fileNameOrbit,      Config::OPTIONAL,  "", "");
    readConfig(config, "outputfileCovariance", fileNameCovariance, Config::OPTIONAL,  "", "");
    readConfig(config, "earthRotation",        earthRotation,      Config::MUSTSET,   "", "");
    readConfig(config, "inputfile",            fileNamesIn,        Config::MUSTSET,   "", "");
    if(isCreateSchema(config)) return;

    // ===============================================

    logStatus<<"read input files"<<Log::endl;
    OrbitArc        orbit;
    Covariance3dArc cov;
    for(const auto &fileName : fileNamesIn)
    {
      try
      {
        logStatus<<"read file <"<<fileName<<">"<<Log::endl;
        InFile file(fileName);

        // skip header
        std::string line;
        std::getline(file, line);
        std::getline(file, line);

        while(std::getline(file, line))
        {
          if(line.empty())
            break;
          std::stringstream ss(line);
          ss.exceptions(std::ios::badbit | std::ios::failbit);

          LongDouble mjd;
          Vector3d   pos;
          Vector     x(6);

          ss>>mjd>>pos.x()>>pos.y()>>pos.z()>>x(0)>>x(1)>>x(2)>>x(3)>>x(4)>>x(5);

          OrbitEpoch        epochOrb;
          Covariance3dEpoch epochCov;
          epochOrb.time = epochCov.time = mjd2time(mjd);
          epochOrb.position = pos;
          epochCov.setData(x);

          orbit.push_back(epochOrb);
          cov.push_back(epochCov);
        }
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    } // for(fileName)

    // Rotation TRF -> CRF
    // -------------------
    if(earthRotation)
    {
      logStatus<<"rotation from TRF to CRF"<<Log::endl;
      Single::forEach(orbit.size(), [&](UInt i)
      {
        const Rotary3d rotation = inverse(earthRotation->rotaryMatrix(orbit.at(i).time));
        orbit.at(i).position = rotation.rotate(orbit.at(i).position);
        cov.at(i).covariance = rotation.rotate(cov.at(i).covariance);
      });
    }

    if(!fileNameOrbit.empty())
    {
      logStatus<<"write orbit data to file <"<<fileNameOrbit<<">"<<Log::endl;
      InstrumentFile::write(fileNameOrbit, orbit);
      Arc::printStatistics(orbit);
    }

    if(!fileNameCovariance.empty())
    {
      logStatus<<"write covariance data to file <"<<fileNameCovariance<<">"<<Log::endl;
      InstrumentFile::write(fileNameCovariance, cov);
      Arc::printStatistics(cov);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
