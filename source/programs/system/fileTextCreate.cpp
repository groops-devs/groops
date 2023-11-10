/***********************************************/
/**
* @file fileTextCreate.cpp
*
* @brief Create text file.
*
* @author Torsten Mayer-Guerr
* @date 2023-05-31
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create text \config{outputfile} containing \config{line}s.
This program can be a powerful tool,
if the \config{line} is repeated with a \configClass{loop}{loopType}
together with the \reference{text parser}{general.parser:text}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"

/***** CLASS ***********************************/

/** @brief Create text file.
* @ingroup programsGroup */
class FileTextCreate
{
 public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(FileTextCreate, SINGLEPROCESS, "Create text file.", System)

/***********************************************/

void FileTextCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileName;
    std::vector<std::string> lines;

    readConfig(config, "outputfile", fileName, Config::MUSTSET,  "", "");
    readConfig(config, "line",       lines,    Config::OPTIONAL, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"Create <"<<fileName<<"> with "<<lines.size()<<" lines"<<Log::endl;
    OutFile file(fileName);
    for(auto &line : lines)
      file<<line<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
