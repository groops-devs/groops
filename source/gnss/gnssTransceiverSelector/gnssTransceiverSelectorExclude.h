/***********************************************/
/**
* @file gnssTransceiverSelectorExclude.h
*
* @brief Exclude transceivers.
* @see GnssTransceiverSelector
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSCEIVERSELECTOREXCLUDE__
#define __GROOPS_GNSSTRANSCEIVERSELECTOREXCLUDE__

// Latex documentation
#ifdef DOCSTRING_GnssTransceiverSelector
static const char *docstringGnssTransceiverSelectorExclude = R"(
\subsection{Exclude}\label{gnssTransceiverSelectorType:exclude}
Select all receivers/transmitters except
\configClass{selector}{gnssTransceiverSelectorType}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"

/***** CLASS ***********************************/

/** @brief Exclude transceivers.
* @ingroup gnssTransceiverSelectorGroup
* @see GnssTransceiverSelector */
class GnssTransceiverSelectorExclude : public GnssTransceiverSelectorBase
{
  GnssTransceiverSelectorPtr selector;

public:
  GnssTransceiverSelectorExclude(Config &config);
  void select(const std::vector<GnssTransceiverPtr> &transceivers, std::vector<Byte> &selected) const override;
  Bool exclude() const override {return TRUE;}
};

/***********************************************/

inline GnssTransceiverSelectorExclude::GnssTransceiverSelectorExclude(Config &config)
{
  try
  {
    readConfig(config, "selector", selector, Config::MUSTSET, "", "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssTransceiverSelectorExclude::select(const std::vector<GnssTransceiverPtr> &transceivers, std::vector<Byte> &selected) const
{
  try
  {
    std::vector<Byte> excl = selector->select(transceivers);
    for(UInt i=0; i<excl.size(); i++)
      if(excl.at(i))
        selected.at(i) = FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
