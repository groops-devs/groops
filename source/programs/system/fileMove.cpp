/***********************************************/
/**
* @file fileMove.cpp
*
* @brief Move/rename file or directory.
*
* @author Torsten Mayer-Guerr
* @date 2023-05-31
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Move/rename file or directory. If the \config{outputfile} is an existing directory
the \config{inputfile} is moved into it.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief Move/rename file or directory.
* @ingroup programsGroup */
class FileMove
{
 public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(FileMove, SINGLEPROCESS, "Move/rename file or directory.", System)

/***********************************************/

void FileMove::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOld, fileNameNew;

    readConfig(config, "outputfile", fileNameNew, Config::MUSTSET, "", "target name or directory for the move/rename");
    readConfig(config, "inputfile",  fileNameOld, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"Move file <"<<fileNameOld<<"> to <"<<fileNameNew<<">"<<Log::endl;
    if(!System::move(fileNameOld, fileNameNew))
      logWarning<<"Unable to move/rename file <"<<fileNameOld<<">"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
