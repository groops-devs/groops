/***********************************************/
/**
* @file gnssReceiverDefinitionCreate.cpp
*
* @brief Create GNSS receiver definition file.
**
* @author Sebastian Strasser
* @date 2018-08-28
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create a \file{GNSS receiver definition file}{gnssReceiverDefinition}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGnssReceiverDefinition.h"

/***** CLASS ***********************************/

/** @brief Create GNSS receiver definition file.
* @ingroup programsGroup */
class GnssReceiverDefinitionCreate
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssReceiverDefinitionCreate, SINGLEPROCESS, "Create GNSS receiver definition file.", Gnss)

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, GnssReceiverDefinitionPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    var = GnssReceiverDefinitionPtr(new GnssReceiverDefinition);

    readConfig(config, "name",     var->name,    Config::MUSTSET,  "",  "");
    readConfig(config, "serial",   var->serial,  Config::OPTIONAL, "",  "");
    readConfig(config, "version",  var->version, Config::OPTIONAL, "",  "");
    readConfig(config, "comment",  var->comment, Config::OPTIONAL, "",  "");
    readConfig(config, "gnssType", var->types,   Config::OPTIONAL, "",  "");
    endSequence(config);

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiverDefinitionCreate::run(Config &config)
{
  try
  {
    FileName fileNameOut;
    std::vector<GnssReceiverDefinitionPtr> receivers;

    readConfig(config, "outputfileGnssReceiverDefinition", fileNameOut, Config::MUSTSET, "", "");
    readConfig(config, "receiverDefinition",               receivers,   Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"write GNSS receiver definition <"<<fileNameOut<<">"<<Log::endl;
    writeFileGnssReceiverDefinition(fileNameOut, receivers);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
