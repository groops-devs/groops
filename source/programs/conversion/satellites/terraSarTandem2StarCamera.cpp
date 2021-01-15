/***********************************************/
/**
* @file terraSarTandem2StarCamera.cpp
*
* @brief read TerraSar-X or Tandem-X star camera data.
*
* @author Norbert Zehentner
* @date 2014-03-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in TerraSar-X or Tandem-X star camera data given in the special format.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read TerraSar-X or Tandem-X star camera data.
* @ingroup programsConversionGroup */
class TerraSarTandem2StarCamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(TerraSarTandem2StarCamera, SINGLEPROCESS, "read TerraSar-X or Tandem-X star camera data", Conversion, Instrument)
GROOPS_RENAMED_PROGRAM(Tsxtdx2Starcamera, TerraSarTandem2StarCamera, date2time(2020, 8, 4))

/***********************************************/

void TerraSarTandem2StarCamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              outNameStar;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileStarCamera", outNameStar, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",            fileNameIn,  Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    StarCameraArc arc;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      try
      {
        logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
        InFile file(fileNameIn.at(idFile));

        std::string line;
        while(std::getline(file, line))
        {
          if(line.empty())
            continue;
          std::stringstream ss(line);
          ss.exceptions(std::ios::badbit | std::ios::failbit);

          Vector q(4);
          Double second;
          ss>>second>>q(1)>>q(2)>>q(3)>>q(0);

          if(std::fabs(norm(q)-1) > 1e-5)
          {
            logWarning<<"strange norm = "<<norm(q)<<Log::endl;
            continue;
          }

          StarCameraEpoch epoch;
          epoch.time   = mjd2time(44244.0) + seconds2time(second);
          epoch.rotary = Rotary3d(q);
          arc.push_back(epoch);
        }
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    }

    logStatus<<"write star camera data to file <"<<outNameStar<<">"<<Log::endl;
    Arc::printStatistics(arc);
    InstrumentFile::write(outNameStar, arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
