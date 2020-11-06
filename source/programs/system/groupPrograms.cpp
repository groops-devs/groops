/***********************************************/
/**
* @file groupPrograms.cpp
*
* @brief Runs programs in group. Sole purpose is to structure GROOPS files.
*
* @author Sebastian Strasser
* @date 2017-09-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Runs programs in a group. The sole purpose is to structure GROOPS config files.
)";

/***********************************************/

#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief  Runs programs in group. Sole purpose is to structure GROOPS files.
* @ingroup programsGroup */
class GroupPrograms
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GroupPrograms, PARALLEL, "Runs programs in group.", System)
GROOPS_RENAMED_PROGRAM(GroupProgramme, GroupPrograms, date2time(2020, 6, 3))

/***********************************************/

void GroupPrograms::run(Config &config)
{
  try
  {
    renameDeprecatedConfig(config, "programme", "program", date2time(2020, 6, 3));

    if(isCreateSchema(config))
    {
      config.xselement("program", "programType", Config::DEFAULT,  Config::UNBOUNDED, "", "");
      return;
    }

    programRun(config);
    programRemove(config);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
