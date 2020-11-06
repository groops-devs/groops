/***********************************************/
/**
* @file fileRemove.cpp
*
* @brief Remove files or directories.
*
* @author Torsten Mayer-Guerr
* @date 2018-05-12
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Remove files or directories.
Deletes also the content recursivley if one of \config{files} is a directory.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief Remove files or directories.
* @ingroup programsGroup */
class FileRemove
{
 public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(FileRemove, SINGLEPROCESS, "Remove files or directories.", System)

/***********************************************/

void FileRemove::run(Config &config)
{
  try
  {
    std::vector<FileName> fileNames;

    readConfig(config, "files", fileNames, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    for(UInt i=0; i<fileNames.size(); i++)
      if(!fileNames.at(i).empty())
      {
        logStatus<<"Remove file <"<<fileNames.at(i)<<">"<<Log::endl;
        if(!System::remove(fileNames.at(i)))
          logWarning<<"Unable to remove file <"<<fileNames.at(i)<<">."<<Log::endl;
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
