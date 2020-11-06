/***********************************************/
/**
* @file gnssParametrizationTransmitterBeidou.h
*
* @brief Beidou satellites (transmitter).
*
* @author Sebastian Strasser
* @date 2018-09-18
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERBEIDOU__
#define __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERBEIDOU__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationTransmitter
static const char *docstringGnssParametrizationTransmitterBeidou = R"(
\subsection{BeiDou}\label{gnssParametrizationTransmitterType:beidou}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssParametrizationTransmitter.h"

/***** CLASS ***********************************/

/** @brief Beidou satellites (transmitter).
* @ingroup gnssParametrizationTransmitterGroup
* @see GnssParametrizationTransmitter */
class GnssParametrizationTransmitterBeidou : public GnssParametrizationTransmitter
{
public:
  GnssParametrizationTransmitterBeidou(Config &config);
  virtual ~GnssParametrizationTransmitterBeidou() {}
  std::string name()   const override {return "beidou";}
  GnssType    system() const override {return GnssType::BDS;}
};

/***********************************************/

#endif /* __GROOPS___ */
