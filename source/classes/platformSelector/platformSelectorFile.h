/***********************************************/
/**
* @file platformSelectorFile.h
*
* @brief Select platforms from file list.
* @see PlatformSelector
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_PLATFORMSELECTORFILE__
#define __GROOPS_PLATFORMSELECTORFILE__

// Latex documentation
#ifdef DOCSTRING_PlatformSelector
static const char *docstringPlatformSelectorFile = R"(
\subsection{File}\label{platformSelectorType:file}
Select receivers/transmitters from each row of
\configFile{inputfileStringTable}{stringTable}.
Additional columns in a row represent alternatives
if previous names are not available (e.g. without observation file).
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileStringTable.h"
#include "classes/platformSelector/platformSelector.h"

/***** CLASS ***********************************/

/** @brief Select platforms from file list.
* @ingroup platformSelectorGroup
* @see PlatformSelector */
class PlatformSelectorFile : public PlatformSelectorBase
{
  std::vector<std::vector<std::string>> names;

public:
  PlatformSelectorFile(Config &config);
  void select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const override;
};

/***********************************************/

inline PlatformSelectorFile::PlatformSelectorFile(Config &config)
{
  try
  {
    FileName fileNameList;
    readConfig(config, "inputfileStringTable", fileNameList, Config::MUSTSET,  "",  "list of names with alternatives");
    readConfig(config, "exclude",              exclude,      Config::DEFAULT,  "0", "deselect first matching platforms");
    if(isCreateSchema(config)) return;

    readFileStringTable(fileNameList, names);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void PlatformSelectorFile::select(const Time &/*timeStart*/, const Time &/*timeEnd*/, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const
{
  try
  {
    for(UInt i=0; i<names.size(); i++)
      for(UInt k=0; k<names.at(i).size(); k++) // alternatives
      {
        const UInt id = std::distance(platforms.begin(), std::find_if(platforms.begin(), platforms.end(),
                                                                      [&](auto &t){return t && (t->name == names.at(i).at(k));}));
        if((id < platforms.size()) && (selected.at(id) == exclude))
        {
          selected.at(id) = !exclude;
          break; // skip alternative stations
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
