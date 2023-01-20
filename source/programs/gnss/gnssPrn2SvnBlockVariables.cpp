/***********************************************/
/**
* @file gnssPrn2SvnBlockVariables.cpp
*
* @brief Create variables containing SVN and block based on a transmitter info file of a GNSS satellite/PRN and a specified time.
*
* @author Torsten Mayer-Guerr
* @date 2017-03-14
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create \reference{variables}{general.parser} containing SVN and block based on an
\configFile{inputfileTransmitterInfo}{platform} of a GNSS satellite/PRN and
a specified \config{time}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief Create variables containing SVN and block based on a transmitter info file of a GNSS satellite/PRN and a specified time.
* @ingroup programsGroup */
class GnssPrn2SvnBlockVariables
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssPrn2SvnBlockVariables, PARALLEL, "Create variables containing SVN and block based on a transmitter info file of a GNSS satellite/PRN and a specified time.", Gnss)

/***********************************************/

void GnssPrn2SvnBlockVariables::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    std::string nameSVN, nameBlock;
    FileName    fileNameTransmitterInfo;
    Time        time;

    readConfig(config, "variableSVN",              nameSVN,                 Config::OPTIONAL,  "svn",   "name of the SVN variable");
    readConfig(config, "variableBlock",            nameBlock,               Config::OPTIONAL,  "block", "name of the satellites block variable");
    readConfig(config, "inputfileTransmitterInfo", fileNameTransmitterInfo, Config::MUSTSET,   "{groopsDataDir}/gnss/transmitter/transmitterInfo/igs/igs14/transmitterInfo_igs14.{prn}.xml", "used for GNSS PRN-to-SVN/model relation");
    readConfig(config, "time",                     time,                    Config::MUSTSET,   "",      "used for GNSS PRN-to-SVN/model relation");
    if(isCreateSchema(config)) return;

    logStatus<<"read transmitter info file <"<<fileNameTransmitterInfo<<">"<<Log::endl;
    Platform transmitterInfo;
    readFilePlatform(fileNameTransmitterInfo, transmitterInfo);
    auto antenna = transmitterInfo.findEquipment<PlatformGnssAntenna>(time);
    if(!antenna)
      throw(Exception("satellite SVN not found in transmitter info"));

    if(!nameSVN.empty())
      addVariable(nameSVN, antenna->serial, config.getVarList());

    if(!nameBlock.empty())
      addVariable(nameBlock, antenna->name, config.getVarList());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
