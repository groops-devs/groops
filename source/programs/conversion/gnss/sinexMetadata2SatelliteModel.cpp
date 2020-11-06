/***********************************************/
/**
* @file sinexMetadata2SatelliteModel.cpp
*
* @brief Create satellite model from IGS SINEX metadata file.
*
* TODO: Consider time-variability in mass and transmit power.
*
* @author Sebastian Strasser
* @date 2018-09-17
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create \configFile{outputfileSatelliteModel}{satelliteModel} from \href{http://mgex.igs.org/IGS_MGEX_Metadata.php}{IGS SINEX metadata format}.

If \configFile{inputfileSatelliteModel}{satelliteModel} is provided it is used as a basis and values are updated from the metadata file.

See also \program{SatelliteModelCreate}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileSatelliteModel.h"
#include "inputOutput/fileSinex.h"

/***** CLASS ***********************************/

/** @brief Create satellite model from IGS SINEX metadata file.
* @ingroup programsConversionGroup */
class SinexMetadata2SatelliteModel
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(SinexMetadata2SatelliteModel, SINGLEPROCESS, "Create satellite model from IGS SINEX metadata file.", Conversion, Gnss)

/***********************************************/

void SinexMetadata2SatelliteModel::run(Config &config)
{
  try
  {
    FileName inNameSinexMetadata, inNameSatelliteModel, outNameSatelliteModel;
    std::string svn;

    readConfig(config, "outputfileSatelliteModel", outNameSatelliteModel, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileSinexMetadata",   inNameSinexMetadata,   Config::MUSTSET,  "", "IGS SINEX metadata file");
    readConfig(config, "inputfileSatelliteModel",  inNameSatelliteModel,  Config::OPTIONAL, "", "base satellite model");
    readConfig(config, "svn",                      svn,                   Config::MUSTSET,  "", "e.g. G040, R736, E204, C211");
    if(isCreateSchema(config)) return;

    Sinex sinex(inNameSinexMetadata);
    std::vector<std::string> massLines  = sinex.getBlock<Sinex::SinexText>("SATELLITE/MASS")->lines();
    std::vector<std::string> powerLines = sinex.getBlock<Sinex::SinexText>("SATELLITE/TX_POWER")->lines();

    SatelliteModelPtr satellite(new SatelliteModel);
    if(!inNameSatelliteModel.empty())
      readFileSatelliteModel(inNameSatelliteModel, satellite);

    satellite->satelliteName = svn;

    auto massIter  = std::find_if(massLines.rbegin(),  massLines.rend(),  [&](const std::string &line){ return (line.size() >= 45 && line.substr(1,4) == svn); });
    if(massIter != massLines.rend())
      satellite->mass = String::toDouble(massIter->substr(36, 9));
    else
      logWarning << "mass missing for " + svn << Log::endl;

    auto powerIter = std::find_if(powerLines.rbegin(), powerLines.rend(), [&](const std::string &line){ return (line.size() >= 40 && line.substr(1,4) == svn); });
    if(powerIter != powerLines.rend())
    {
      auto iter = std::remove_if(satellite->modules.begin(), satellite->modules.end(), [&](SatelliteModelModulePtr module){ return module->type() == SatelliteModelModule::ANTENNATHRUST; });
      satellite->modules.erase(iter, satellite->modules.end());
      SatelliteModelModuleAntennaThrust *module = new SatelliteModelModuleAntennaThrust;
      module->thrust = Vector3d(0, 0, String::toInt(powerIter->substr(36, 4)));
      satellite->modules.push_back(SatelliteModelModulePtr(module));
    }
    else
      logWarning << "transmit power missing for " + svn << Log::endl;

    if(!satellite->surfaces.size())
    {
      // Box-wing surfaces
      satellite->surfaces.resize(8);
      satellite->surfaces.at(0).normal = Vector3d( 0, 0, 1); //  Z (TOWARDS THE EARTH)
      satellite->surfaces.at(1).normal = Vector3d( 0, 0,-1); // -Z
      satellite->surfaces.at(2).normal = Vector3d( 0, 1, 0); //  Y (ALONG SOLAR PANELS BEAMS)
      satellite->surfaces.at(3).normal = Vector3d( 0,-1, 0); // -Y
      satellite->surfaces.at(4).normal = Vector3d( 1, 0, 0); //  X (ALONG BUS DIRECTION ALWAYS ILLUMINATED BY THE SUN)
      satellite->surfaces.at(5).normal = Vector3d(-1, 0, 0); // -X
      satellite->surfaces.at(6).normal = Vector3d( 0, 0, 1); // SOLAR PANELS
      satellite->surfaces.at(7).normal = Vector3d( 0, 0,-1);

      // Rotation of solar panels
      SatelliteModelModuleSolarPanel *module = new SatelliteModelModuleSolarPanel;
      satellite->modules.push_back(SatelliteModelModulePtr(module));
      module->rotationAxis = Vector3d(0,1,0);
      module->normal       = Vector3d(0,0,1);
      module->indexSurface.push_back(6);
      module->indexSurface.push_back(7);
    }

    logStatus<<"write satellite model file <"<<outNameSatelliteModel<<">"<<Log::endl;
    writeFileSatelliteModel(outNameSatelliteModel, satellite);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
