/***********************************************/
/**
* @file fileCreateDirectories.cpp
*
* @brief Creates the directory and parent directories as needed.
*
* @author Torsten Mayer-Guerr
* @date 2018-05-12
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Creates the directory and parent directories as needed.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief Creates the directory and parent directories as needed.
* @ingroup programsGroup */
class FileCreateDirectories
{
 public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(FileCreateDirectories, SINGLEPROCESS, "Creates the directory and parent directories as needed.", System)

/***********************************************/

void FileCreateDirectories::run(Config &config)
{
  try
  {
    std::vector<FileName> fileNames;

    readConfig(config, "directory", fileNames, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    for(UInt i=0; i<fileNames.size(); i++)
    {
      logStatus<<"Create directory <"<<fileNames.at(i)<<">"<<Log::endl;
      if(!System::createDirectories(fileNames.at(i)))
        throw(Exception("Cannot create directory <"+fileNames.at(i).str()+">"));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
