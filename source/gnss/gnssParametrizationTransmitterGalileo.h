/***********************************************/
/**
* @file gnssParametrizationTransmitterGalileo.h
*
* @brief Galileo satellites (transmitter).
*
* @author Sebastian Strasser
* @date 2018-09-18
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERGALILEO__
#define __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERGALILEO__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationTransmitter
static const char *docstringGnssParametrizationTransmitterGalileo = R"(
\subsection{Galileo}\label{gnssParametrizationTransmitterType:galileo}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssParametrizationTransmitter.h"

/***** CLASS ***********************************/

/** @brief Galileo satellites (transmitter).
* @ingroup gnssParametrizationTransmitterGroup
* @see GnssParametrizationTransmitter */
class GnssParametrizationTransmitterGalileo : public GnssParametrizationTransmitter
{
public:
  GnssParametrizationTransmitterGalileo(Config &config);
  virtual ~GnssParametrizationTransmitterGalileo() {}
  std::string name()   const override {return "galileo";}
  GnssType    system() const override {return GnssType::GALILEO;}
};

/***********************************************/

#endif /* __GROOPS___ */
