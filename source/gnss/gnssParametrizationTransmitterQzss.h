/***********************************************/
/**
* @file gnssParametrizationTransmitterQzss.h
*
* @brief QZSS satellites (transmitter).
*
* @author Sebastian Strasser
* @date 2021-02-08
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERQZSS__
#define __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERQZSS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationTransmitter
static const char *docstringGnssParametrizationTransmitterQzss = R"(
\subsection{QZSS}\label{gnssParametrizationTransmitterType:qzss}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssParametrizationTransmitter.h"

/***** CLASS ***********************************/

/** @brief QZSS satellites (transmitter).
* @ingroup gnssParametrizationTransmitterGroup
* @see GnssParametrizationTransmitter */
class GnssParametrizationTransmitterQzss : public GnssParametrizationTransmitter
{
public:
  GnssParametrizationTransmitterQzss(Config &config);
  virtual ~GnssParametrizationTransmitterQzss() {}
  std::string name()   const override {return "qzss";}
  GnssType    system() const override {return GnssType::QZSS;}
};

/***********************************************/

#endif /* __GROOPS___ */
