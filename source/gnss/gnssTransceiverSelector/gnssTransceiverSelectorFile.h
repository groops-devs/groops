/***********************************************/
/**
* @file gnssTransceiverSelectorFile.h
*
* @brief Select transceivers from file list.
* @see GnssTransceiverSelector
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSCEIVERSELECTORFILE__
#define __GROOPS_GNSSTRANSCEIVERSELECTORFILE__

// Latex documentation
#ifdef DOCSTRING_GnssTransceiverSelector
static const char *docstringGnssTransceiverSelectorFile = R"(
\subsection{File}\label{gnssTransceiverSelectorType:file}
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
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"

/***** CLASS ***********************************/

/** @brief Select transceivers from file list.
* @ingroup gnssTransceiverSelectorGroup
* @see GnssTransceiverSelector */
class GnssTransceiverSelectorFile : public GnssTransceiverSelectorBase
{
  std::vector<std::vector<std::string>> names;

public:
  GnssTransceiverSelectorFile(Config &config);
  void select(const std::vector<GnssTransceiverPtr> &transceivers, std::vector<Byte> &selected) const override;
};

/***********************************************/

inline GnssTransceiverSelectorFile::GnssTransceiverSelectorFile(Config &config)
{
  try
  {
    FileName fileNameList;
    readConfig(config, "inputfileStringTable", fileNameList, Config::MUSTSET, "", "list of names with alternatives");
    if(isCreateSchema(config)) return;

    readFileStringTable(fileNameList, names);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssTransceiverSelectorFile::select(const std::vector<GnssTransceiverPtr> &transceivers, std::vector<Byte> &selected) const
{
  try
  {
    for(UInt i=0; i<names.size(); i++)
      for(UInt k=0; k<names.at(i).size(); k++) // alternatives
      {
        const UInt id = std::distance(transceivers.begin(), std::find_if(transceivers.begin(), transceivers.end(),
                                                                         [&](auto &t){return t->name() == names.at(i).at(k);}));
        if(id < transceivers.size())
        {
          selected.at(id) = transceivers.at(id)->useable();
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
