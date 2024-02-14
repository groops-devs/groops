/***********************************************/
/**
* @file gnssPrn2SvnBlockVariables.cpp
*
* @brief DEPRECATED. This program will no longer work from the next release! See documentation for help.
*
* @author Torsten Mayer-Guerr
* @date 2017-03-14
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
DEPRECATED. This program will no longer work from the next release!

Setup up a \configClass{loop:platformEquipment}{loopType:platformEquipment} instead with
\begin{itemize}
  \item \configFile{inputfilePlatform}{platform}: the old \config{inputfileTransmitterInfo}
  \item \config{equipmentType}         = \verb|gnssAntenna|
  \item \config{variableLoopName}      = \verb|block|
  \item \config{variableLoopSerial}    = \verb|svn|
  \item \config{variableLoopTimeStart} = \verb|svnTimeStart|
  \item \config{variableLoopTimeEnd}   = \verb|svnTimeEnd|
  \item \configClass{condition:expression}{conditionType:expression}
  \begin{itemize}
    \item \config{expression} = \verb|(svnTimeStart <= time) && (time < svnTimeEnd)|
  \end{itemize}
\end{itemize}
Attribute this loop to programs, which uses the variables.
)";

/***********************************************/

#include "programs/program.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief DEPRECATED. This program no longer works! See documentation for help.
* @ingroup programsGroup */
class GnssPrn2SvnBlockVariables
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssPrn2SvnBlockVariables, PARALLEL, "This program no longer works! See documentation for help.", Deprecated)

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

    logWarningOnce<<"This program will no longer work from the next release!!!!! See documentation for help."<<Log::endl;
    logWarningOnce<<"  Setup up a loop->platformEquipment instead with"<<Log::endl;
    logWarningOnce<<"    - inputfilePlatform     = inputfileTransmitterInfo"<<Log::endl;
    logWarningOnce<<"    - equipmentType         = gnssAntenna"<<Log::endl;
    logWarningOnce<<"    - variableLoopName      = block"<<Log::endl;
    logWarningOnce<<"    - variableLoopSerial    = svn"<<Log::endl;
    logWarningOnce<<"    - variableLoopTimeStart = svnTimeStart"<<Log::endl;
    logWarningOnce<<"    - variableLoopTimeEnd   = svnTimeEnd"<<Log::endl;
    logWarningOnce<<"    - condition->expression"<<Log::endl;
    logWarningOnce<<"      -- expression = (svnTimeStart <= time) && (time < svnTimeEnd)"<<Log::endl;
    logWarningOnce<<"  Attribute the loop to programs, which uses the variables."<<Log::endl;

    logStatus<<"read transmitter info file <"<<fileNameTransmitterInfo<<">"<<Log::endl;
    Platform transmitterInfo;
    readFilePlatform(fileNameTransmitterInfo, transmitterInfo);
    auto antenna = transmitterInfo.findEquipment<PlatformGnssAntenna>(time);
    if(!antenna)
      throw(Exception(fileNameTransmitterInfo.str()+" contains no satellite at "+time.dateTimeStr()));

    if(!nameSVN.empty())
      config.getVarList().setVariable(nameSVN, antenna->serial);

    if(!nameBlock.empty())
      config.getVarList().setVariable(nameBlock, antenna->name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
