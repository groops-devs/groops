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
\configFile{inputfileTransmitterInfo}{gnssStationInfo} of a GNSS satellite/PRN and
a specified \config{time}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGnssStationInfo.h"

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
    readConfig(config, "inputfileTransmitterInfo", fileNameTransmitterInfo, Config::MUSTSET,   "",      "used for GNSS PRN-to-SVN/model relation");
    readConfig(config, "time",                     time,                    Config::MUSTSET,   "",      "used for GNSS PRN-to-SVN/model relation");
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
