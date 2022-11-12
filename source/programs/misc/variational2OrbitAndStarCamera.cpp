/***********************************************/
/**
* @file variational2OrbitAndStarCamera.cpp
*
* @brief Extracts orbit and star camera data from variational file.
*
* @author Sebastian Strasser
* @author Andreas Kvas
* @date 2016-07-13
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Extracts the reference \configFile{outputfileOrbit}{instrument}, \configFile{outputfileStarCamera}{instrument},
and \configFile{outputfileEarthRotation}{instrument} from \configFile{inputfileVariational}{variationalEquation}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileVariationalEquation.h"

/***** CLASS ***********************************/

/** @brief Extracts orbit and star camera data from variational file.
* @ingroup programsGroup */
class Variational2OrbitAndStarCamera
{
 public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Variational2OrbitAndStarCamera, SINGLEPROCESS, "Extracts orbit and star camera data from variational file.", Misc, VariationalEquation, Instrument)
GROOPS_RENAMED_PROGRAM(Variational2Orbit, Variational2OrbitAndStarCamera, date2time(2022, 10, 7))

/***********************************************/

void Variational2OrbitAndStarCamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutOrbit, fileNameOutStarCamera, fileNameOutEarthRotation, fileNameInVariational;

    renameDeprecatedConfig(config, "outputFileOrbit",      "outputfileOrbit",      date2time(2020, 7, 15));
    renameDeprecatedConfig(config, "inputFileVariational", "inputfileVariational", date2time(2020, 7, 15));

    readConfig(config, "outputfileOrbit",         fileNameOutOrbit,         Config::OPTIONAL, "", "output orbit (instrument) file");
    readConfig(config, "outputfileStarCamera",    fileNameOutStarCamera,    Config::OPTIONAL, "", "output satellite attidude as star camera (instrument) file");
    readConfig(config, "outputfileEarthRotation", fileNameOutEarthRotation, Config::OPTIONAL, "", "output Earth rotation as star camera (instrument) file");
    readConfig(config, "inputfileVariational",    fileNameInVariational,    Config::MUSTSET,  "", "input variational file");
    if(isCreateSchema(config)) return;

    if(fileNameOutOrbit.empty() && fileNameOutStarCamera.empty() && fileNameOutEarthRotation.empty())
      return;

    FileVariationalEquation fileVariational(fileNameInVariational);
    std::list<Arc> arcListOrbit, arcListStarCamera, arcListEarthRotation;
    for(UInt arcNo = 0; arcNo < fileVariational.arcCount(); arcNo++)
    {
      VariationalEquationArc arc = fileVariational.readArc(arcNo);
      arcListOrbit.push_back(arc.orbitArc());

      StarCameraArc starCamera;
      for(UInt i=0; i<arc.times.size(); i++)
      {
        StarCameraEpoch e;
        e.time   = arc.times.at(i);
        e.rotary = arc.rotSat.at(i);
        starCamera.push_back(e);
      }
      arcListStarCamera.push_back(starCamera);

      StarCameraArc earthRotation;
      for(UInt i=0; i<arc.times.size(); i++)
      {
        StarCameraEpoch e;
        e.time   = arc.times.at(i);
        e.rotary = arc.rotEarth.at(i);
        earthRotation.push_back(e);
      }
      arcListEarthRotation.push_back(earthRotation);
    }

    if(!fileNameOutOrbit.empty())
    {
      logStatus<<"write orbit to file <"<<fileNameOutOrbit<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutOrbit, arcListOrbit);
    }
    if(!fileNameOutStarCamera.empty())
    {
      logStatus<<"write star camera to file <"<<fileNameOutStarCamera<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutStarCamera, arcListStarCamera);
    }
    if(!fileNameOutEarthRotation.empty())
    {
      logStatus<<"write Earth rotation to file <"<<fileNameOutEarthRotation<<">"<< Log::endl;
      InstrumentFile::write(fileNameOutEarthRotation, arcListEarthRotation);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
