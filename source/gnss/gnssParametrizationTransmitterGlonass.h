/***********************************************/
/**
* @file gnssParametrizationTransmitterGlonass.h
*
* @brief GLONASS satellites (transmitter).
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2018-08-08
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERGLONASS__
#define __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERGLONASS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationTransmitter
static const char *docstringGnssParametrizationTransmitterGlonass = R"(
\subsection{GLONASS}\label{gnssParametrizationTransmitterType:glonass}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssParametrizationTransmitter.h"

/***** CLASS ***********************************/

/** @brief GLONASS satellites (transmitter).
* @ingroup gnssParametrizationTransmitterGroup
* @see GnssParametrizationTransmitter */
class GnssParametrizationTransmitterGlonass : public GnssParametrizationTransmitter
{
public:
  GnssParametrizationTransmitterGlonass(Config &config);
  virtual ~GnssParametrizationTransmitterGlonass() {}
  std::string name()   const override {return "glonass";}
  GnssType    system() const override {return GnssType::GLONASS;}
};

/***********************************************/

#endif /* __GROOPS___ */
