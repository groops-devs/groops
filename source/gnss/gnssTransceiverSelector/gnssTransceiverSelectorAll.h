/***********************************************/
/**
* @file gnssTransceiverSelectorAll.h
*
* @brief All available stations.
* @see GnssTransceiverSelector
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSCEIVERSELECTORALL__
#define __GROOPS_GNSSTRANSCEIVERSELECTORALL__

// Latex documentation
#ifdef DOCSTRING_GnssTransceiverSelector
static const char *docstringGnssTransceiverSelectorAll = R"(
\subsection{All}\label{gnssTransceiverSelectorType:all}
Select all receivers/transmitters.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"

/***** CLASS ***********************************/

/** @brief All available stations.
* @ingroup gnssTransceiverSelectorGroup
* @see GnssTransceiverSelector */
class GnssTransceiverSelectorAll : public GnssTransceiverSelectorBase
{
public:
  GnssTransceiverSelectorAll(Config &/*config*/) {}

  void select(const std::vector<GnssTransceiverPtr> &transceivers, std::vector<Byte> &selected) const override;
};

/***********************************************/

inline void GnssTransceiverSelectorAll::select(const std::vector<GnssTransceiverPtr> &/*transceivers*/, std::vector<Byte> &selected) const
{
  std::fill(selected.begin(), selected.end(), TRUE);
}

/***********************************************/

#endif
