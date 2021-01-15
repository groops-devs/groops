/***********************************************/
/**
* @file satelliteModelCreate.cpp
*
* @brief Create satellite macro model.
*
* @author Torsten Mayer-Guerr
* @author Sandro Krauss
* @date 2015-27-05
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program creates a satellite macro model for the estimation of non-gravitational accelerations acting on a satellite.
Mandatory input values are the \config{satelliteName}, \config{mass}, \config{coefficientDrag} and information
about the satellite \config{surfaces}. For low Earth orbiting satellites, like GRACE for instance, a good guess
for the drag coefficient could be 2.3. Apart from that, it is latter on possible to estimate a more precise variable drag coefficient
(e.g. \configClass{miscAccelerations:atmosphericDrag}{miscAccelerationsType:atmosphericDrag}), which will override this initial guess.
Concerning the satellite surfaces an external file must be imported which must contain information about each single
 satellite plate in terms of plate \config{area}, the associated plate normal and re-radiation properties
(reflexion, diffusion and absorption) properties in the visible and IR part. Examplarily, a description of
the macro model for GRACE can be found under:
\url{https://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/docs/ProdSpecDoc_v4.6.pdf}
Additionally, it is possible to add further information like antennaThrust, solar panel, temporal mass changes and
massInstrument using the modules option.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"

/***** CLASS ***********************************/

/** @brief Create satellite macro model.
* @ingroup programsGroup */
class SatelliteModelCreate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SatelliteModelCreate, SINGLEPROCESS, "Create satellite macro model.", Misc)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, SatelliteModelModulePtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    std::string choice;
    if(!readConfigChoice(config, name, choice, mustSet, defaultValue, annotation))
      return FALSE;

    // Antenna thrust
    // --------------
    if(readConfigChoiceElement(config, "antennaThrust", choice, ""))
    {
      SatelliteModelModuleAntennaThrust *modul = new SatelliteModelModuleAntennaThrust;
      var = SatelliteModelModulePtr(modul);

      readConfig(config, "thrustX", modul->thrust.x(), Config::DEFAULT,  "0", "");
      readConfig(config, "thrustY", modul->thrust.y(), Config::DEFAULT,  "0", "");
      readConfig(config, "thrustZ", modul->thrust.z(), Config::DEFAULT,  "0", "");
    }

    // Solar Panel
    // -----------
    if(readConfigChoiceElement(config, "solarPanel", choice, ""))
    {
      SatelliteModelModuleSolarPanel *modul = new SatelliteModelModuleSolarPanel;
      var = SatelliteModelModulePtr(modul);

      readConfig(config, "rotationAxisX", modul->rotationAxis.x(), Config::DEFAULT,  "0", "");
      readConfig(config, "rotationAxisY", modul->rotationAxis.y(), Config::DEFAULT,  "0", "");
      readConfig(config, "rotationAxisZ", modul->rotationAxis.z(), Config::DEFAULT,  "0", "");
      readConfig(config, "normalX",       modul->normal.x(),       Config::DEFAULT,  "0", "Direction to sun");
      readConfig(config, "normalY",       modul->normal.y(),       Config::DEFAULT,  "0", "Direction to sun");
      readConfig(config, "normalZ",       modul->normal.z(),       Config::DEFAULT,  "0", "Direction to sun");
      readConfig(config, "indexSurface",  modul->indexSurface,     Config::MUSTSET,   "", "index of solar panel surfaces");
    }

    // Satellite mass change
    // ---------------------
    if(readConfigChoiceElement(config, "massChange", choice, ""))
    {
      SatelliteModelModuleMassChange *modul = new SatelliteModelModuleMassChange;
      var = SatelliteModelModulePtr(modul);

      readConfig(config, "time", modul->times, Config::MUSTSET,  "", "");
      readConfig(config, "mass", modul->mass,  Config::MUSTSET,  "", "");
    }

    // Satellite mass change
    // ---------------------
    if(readConfigChoiceElement(config, "massInstrument", choice, ""))
    {
      std::vector<FileName> fileName;
      readConfig(config, "inputfileInstrument", fileName, Config::MUSTSET,  "", "");

      if(!isCreateSchema(config))
      {
        MassArc massArc;
        for(UInt i=0; i<fileName.size(); i++)
          massArc.append(InstrumentFile::read(fileName.at(i)));
        massArc.sort();

        SatelliteModelModuleMassChange *modul = new SatelliteModelModuleMassChange;
        var = SatelliteModelModulePtr(modul);
        for(UInt i=0; i<massArc.size(); i++)
          if(massArc.at(i).massThr > 0.0)
          {
            modul->times.push_back(massArc.at(i).time);
            modul->mass.push_back(massArc.at(i).massThr);
          }
      }
    }

    endChoice(config);

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, std::vector<SatelliteModel::Surface> &surfaces, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    FileName fileNameIn;
    ExpressionVariablePtr typeExpr, areaExpr;
    ExpressionVariablePtr normalXExpr, normalYExpr, normalZExpr;
    ExpressionVariablePtr absorptionVisibleExpr, diffusionVisibleExpr, reflexionVisibleExpr;
    ExpressionVariablePtr absorptionInfraredExpr, diffusionInfraredExpr, reflexionInfraredExpr;
    ExpressionVariablePtr hasThermalReemissionExpr;

    readConfig(config, "inputfile",            fileNameIn,               Config::MUSTSET, "",       "each line must contain one surface element");
    readConfig(config, "type",                 typeExpr,                 Config::MUSTSET, "0",      "0: plate, 1: sphere, 2: cylinder");
    readConfig(config, "area",                 areaExpr,                 Config::MUSTSET, "data0",  "[m**2]");
    readConfig(config, "normalX",              normalXExpr,              Config::MUSTSET, "data1",  "");
    readConfig(config, "normalY",              normalYExpr,              Config::MUSTSET, "data2",  "");
    readConfig(config, "normalZ",              normalZExpr,              Config::MUSTSET, "data3",  "");
    readConfig(config, "reflexionVisible",     reflexionVisibleExpr,     Config::MUSTSET, "data4",  "");
    readConfig(config, "diffusionVisible",     diffusionVisibleExpr,     Config::MUSTSET, "data5",  "");
    readConfig(config, "absorptionVisible",    absorptionVisibleExpr,    Config::MUSTSET, "data6",  "");
    readConfig(config, "reflexionInfrared",    reflexionInfraredExpr,    Config::MUSTSET, "data7",  "");
    readConfig(config, "diffusionInfrared",    diffusionInfraredExpr,    Config::MUSTSET, "data8",  "");
    readConfig(config, "absorptionInfrared",   absorptionInfraredExpr,   Config::MUSTSET, "data9",  "");
    readConfig(config, "hasThermalReemission", hasThermalReemissionExpr, Config::MUSTSET, "data10", "(0 = no, 1 = yes) spontaneous reemission of absorbed radiation");
    endSequence(config);
    if(isCreateSchema(config))
      return TRUE;

    // read data
    // ---------
    logStatus<<"read input from <"<<fileNameIn<<">"<<Log::endl;
    Matrix data;
    readFileMatrix(fileNameIn, data);

    // calculate output
    // -----------------
    auto varList = config.getVarList();
    std::set<std::string> usedVariables;
    typeExpr                ->usedVariables(varList, usedVariables);
    areaExpr                ->usedVariables(varList, usedVariables);
    normalXExpr             ->usedVariables(varList, usedVariables);
    normalYExpr             ->usedVariables(varList, usedVariables);
    normalZExpr             ->usedVariables(varList, usedVariables);
    absorptionVisibleExpr   ->usedVariables(varList, usedVariables);
    diffusionVisibleExpr    ->usedVariables(varList, usedVariables);
    reflexionVisibleExpr    ->usedVariables(varList, usedVariables);
    absorptionInfraredExpr  ->usedVariables(varList, usedVariables);
    diffusionInfraredExpr   ->usedVariables(varList, usedVariables);
    reflexionInfraredExpr   ->usedVariables(varList, usedVariables);
    hasThermalReemissionExpr->usedVariables(varList, usedVariables);
    addDataVariables(data, varList, usedVariables);
    surfaces.resize(data.rows());
    for(UInt i=0; i<data.rows(); i++)
    {
      evaluateDataVariables(data, i, varList);

      // evaluate expression
      surfaces.at(i).type                 = static_cast<SatelliteModel::Surface::Type>(typeExpr->evaluate(varList));
      surfaces.at(i).area                 = areaExpr              ->evaluate(varList);
      surfaces.at(i).normal.x()           = normalXExpr           ->evaluate(varList);
      surfaces.at(i).normal.y()           = normalYExpr           ->evaluate(varList);
      surfaces.at(i).normal.z()           = normalZExpr           ->evaluate(varList);
      surfaces.at(i).absorptionVisible    = absorptionVisibleExpr ->evaluate(varList);
      surfaces.at(i).diffusionVisible     = diffusionVisibleExpr  ->evaluate(varList);
      surfaces.at(i).reflexionVisible     = reflexionVisibleExpr  ->evaluate(varList);
      surfaces.at(i).absorptionInfrared   = absorptionInfraredExpr->evaluate(varList);
      surfaces.at(i).diffusionInfrared    = diffusionInfraredExpr ->evaluate(varList);
      surfaces.at(i).reflexionInfrared    = reflexionInfraredExpr ->evaluate(varList);
      surfaces.at(i).hasThermalReemission = static_cast<Bool>(hasThermalReemissionExpr->evaluate(varList));

      surfaces.at(i).normal.normalize();
    }

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, SatelliteModelPtr &satellite, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    satellite = SatelliteModelPtr(new SatelliteModel);

    renameDeprecatedConfig(config, "modul", "module", date2time(2020, 8, 20));

    readConfig(config, "satelliteName",   satellite->satelliteName,   Config::MUSTSET,  "", "");
    readConfig(config, "mass",            satellite->mass,            Config::MUSTSET,  "", "");
    readConfig(config, "coefficientDrag", satellite->coefficientDrag, Config::MUSTSET,  "", "");
    readConfig(config, "surfaces",        satellite->surfaces,        Config::MUSTSET,  "", "");
    readConfig(config, "module",          satellite->modules,         Config::OPTIONAL, "", "");
    endSequence(config);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SatelliteModelCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut;
    std::vector<SatelliteModelPtr> satellites;

    readConfig(config, "outputfileSatelliteModel", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "satellite",                satellites,  Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // write file
    // ----------
    logStatus<<"write satellite model file <"<<fileNameOut<<">"<<Log::endl;
    writeFileSatelliteModel(fileNameOut, satellites);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
