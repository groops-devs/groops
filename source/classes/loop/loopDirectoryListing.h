/***********************************************/
/**
* @file loopDirectoryListing.h
*
* @brief Loop over files of a directory.
*
* @author Torsten Mayer.Guerr
* @date 2017-02-10
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPDIRECTORYLISTING__
#define __GROOPS_LOOPDIRECTORYLISTING__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopDirectoryListing = R"(
\subsection{DirectoryListing}\label{loopType:directoryListing}
Loop over files of a directory.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "parallel/parallel.h"
#include "inputOutput/logging.h"
#include "inputOutput/system.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over files of a directory.
* @ingroup LoopGroup
* @see Loop */
class LoopDirectoryListing : public Loop
{
  std::string           nameFile, nameIndex, nameCount;
  std::vector<FileName> files;
  UInt                  index;

public:
  LoopDirectoryListing(Config &config);

  UInt count() const override {return files.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopDirectoryListing::LoopDirectoryListing(Config &config)
{
  try
  {
    FileName    directory;
    std::string pattern;
    Bool        isRegularExpression;

    readConfig(config, "directory",           directory,           Config::MUSTSET,  "",         "directory");
    readConfig(config, "pattern",             pattern,             Config::DEFAULT,  "*",        "wildcard pattern");
    readConfig(config, "isRegularExpression", isRegularExpression, Config::DEFAULT,  "0",        "pattern is a regular expression");
    readConfig(config, "variableLoopFile",    nameFile,            Config::OPTIONAL, "loopFile", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex",   nameIndex,           Config::OPTIONAL, "",         "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",   nameCount,           Config::OPTIONAL, "",         "variable with total number of iterations");
    if(isCreateSchema(config)) return;

    files = System::directoryListing(directory, isRegularExpression ? std::regex(pattern) : String::wildcard2regex(pattern));
    std::sort(files.begin(), files.end());
    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopDirectoryListing::iteration(VariableList &varList)
{
  if(index >= count())
    return FALSE;

  if(!nameFile.empty())  varList.setVariable(nameFile,  files.at(index).str());
  if(!nameIndex.empty()) varList.setVariable(nameIndex, index);
  if(!nameCount.empty()) varList.setVariable(nameCount, count());

  index++;
  return TRUE;
}

/***********************************************/

#endif
