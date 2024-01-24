/***********************************************/
/**
* @file gnssStationIdShort2Long.cpp
*
* @brief Create variable containing long (9 char) site ID based on short (4 char) site ID.
*
* @author André Hauschild
* @date 2023-04-26
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create \reference{variables}{general.parser} containing long (9 char) site ID based on short (4 char) site ID.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileStringTable.h"

#include <map>
#include <string>
#include <fstream>

/***** CLASS ***********************************/

/** @brief Create variable containing long (9 char) site ID based on short (4 char) site ID.
* @ingroup programsGroup */
class GnssStationIdShort2Long
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssStationIdShort2Long, PARALLEL, "Create variable containing long (9 char) site ID based on short (4 char) site ID.", Gnss)

/***********************************************/

void GnssStationIdShort2Long::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {

     std::string site4char;
    std::string site9char;
    FileName    fileNameSiteID;

    readConfig(config, "variableSite9char", site9char,      Config::OPTIONAL, "site9c",   "9 character site ID");
    readConfig(config, "site4char",         site4char,      Config::MUSTSET,  "site4c",   "4 character site ID");
    readConfig(config, "inputfileSiteID",   fileNameSiteID, Config::MUSTSET,  "{groopsDataDir}/SiteId4char9char.txt", "used for site ID conversion");
    if(isCreateSchema(config)) return;

    logStatus<<"read site ID conversion file <"<<fileNameSiteID<<">"<<Log::endl;

    std::vector<std::vector<std::string>> siteIDs;

    readFileStringTable(fileNameSiteID, siteIDs);

    auto iter = std::find_if(siteIDs.begin(), siteIDs.end(),
                             [&](const auto &x){return x.size() && (x.at(0) == site4char);});
    if(iter == siteIDs.end())
      throw(Exception(fileNameSiteID.str()+" contains no site ID <"+site4char+">"));

    config.getVarList().setVariable(site9char, iter->at(1));

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
