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

    std::string line;
    std::map<std::string,std::string> siteIDs;

    std::ifstream inp(fileNameSiteID.c_str());
    if (inp.fail()) {
      throw(Exception(fileNameSiteID.str()+" cannot open input file!"));
    };
    while(!inp.eof()){

      getline(inp,line);
      if (line.size()<14) continue;

      std::string site4c = line.substr(0, 4);
      std::string site9c = line.substr(5,14);

      siteIDs[site4c] = site9c;

      //logStatus <<"  found <" << site4c << "> -> <" << site9c << ">" << Log::endl;

    };
    inp.close();

    if (siteIDs.find(site4char)==siteIDs.end())
      throw(Exception(fileNameSiteID.str()+" contains no site ID <"+site4char+">"));

    if(!site9char.empty())
      addVariable(site9char, siteIDs.find(site4char)->second, config.getVarList());

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
