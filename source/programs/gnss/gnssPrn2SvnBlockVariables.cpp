/***********************************************/
/**
* @file gnssPrn2SvnBlockVariables.cpp
*
* @brief Create variable of SVN and block for a GNSS PRN at a specific time.
*
* @author Torsten Mayer-Guerr
* @date 2017-03-14
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create \reference{variables}{general.parser} containing SVN and block for a GNSS PRN at a specific \config{time}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGnssStationInfo.h"

/***** CLASS ***********************************/

/** @brief Create variable of SVN and block for a GNSS PRN at a specific time.
* @ingroup programsGroup */
class GnssPrn2SvnBlockVariables
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssPrn2SvnBlockVariables, PARALLEL, "Create variable of SVN and block for a GNSS PRN at a specific time.", Gnss)

/***********************************************/

void GnssPrn2SvnBlockVariables::run(Config &config)
{
  try
  {
    std::string nameSVN, nameBlock;
    FileName    fileNameTransmitterInfo;
    Time        time;

    readConfig(config, "variableSVN",              nameSVN,                 Config::OPTIONAL,  "svn",   "name of the SVN variable");
    readConfig(config, "variableBlock",            nameBlock,               Config::OPTIONAL,  "block", "name of the satellites block variable");
    readConfig(config, "inputfileTransmitterInfo", fileNameTransmitterInfo, Config::MUSTSET,   "{groopsDataDir}/gnss/transmitterGPS/transmitterInfo/igs/igs14/transmitterInfo_igs14.{prn}.xml", "used for GNSS PRN-to-SVN/model relation");
    readConfig(config, "time",                     time,                    Config::MUSTSET,   "", "used for GNSS PRN-to-SVN/model relation");
    if(isCreateSchema(config)) return;

    logStatus<<"read transmitter info file <"<<fileNameTransmitterInfo<<">"<<Log::endl;
    GnssStationInfo transmitterInfo;
    readFileGnssStationInfo(fileNameTransmitterInfo, transmitterInfo);
    const UInt idAntenna = transmitterInfo.findAntenna(time);
    if(idAntenna == NULLINDEX)
      throw(Exception("satellite SVN not found in transmitter info"));

    if(!nameSVN.empty())
      addVariable(nameSVN, transmitterInfo.antenna.at(idAntenna).serial, config.getVarList());

    if(!nameBlock.empty())
      addVariable(nameBlock, transmitterInfo.antenna.at(idAntenna).name, config.getVarList());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
