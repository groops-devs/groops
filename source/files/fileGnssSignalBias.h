/***********************************************/
/**
* @file fileGnssSignalBias.h
*
* @brief Code/Phase biases.
*
* @author Torsten Mayer-Guerr
* @date 2013-08-11
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSSIGNALBIAS__
#define __GROOPS_GNSSSIGNALBIAS__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_GnssSignalBias
static const char *docstringGnssSignalBias = R"(
Signal biases of GNSS transmitters or receivers for different \configClass{gnssType}{gnssType}.

\begin{verbatim}
groops gnssSignalBias version=20200123
          5 # number of signals
# type   bias [m]
# ===============================
 C1CG06 -1.752461109688110974e-01
 C1WG06  4.005884595055994590e-02
 C2WG06  6.597469378913034532e-02
 L1*G06 -2.736169875580296909e-02
 L2*G06  3.422596762686257871e-02
 \end{verbatim}

See also \program{GnssProcessing}, \program{GnssSimulateReceiver}, \program{GnssSignalBias2Matrix}, \program{GnssSignalBias2SinexBias}.
)";
#endif

/***********************************************/

#include "base/gnssType.h"
#include "inputOutput/fileName.h"
#include "inputOutput/archive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_GNSSSIGNALBIAS_TYPE = "gnssSignalBias";

/***** TYPES ***********************************/

class GnssSignalBias;
typedef std::shared_ptr<GnssSignalBias> GnssSignalBiasPtr;

/***** CLASS ***********************************/

/** @brief Code/Phase biases. */
class GnssSignalBias
{
  public:
  std::vector<GnssType> type;
  std::vector<Double>   bias;

  Vector compute(const std::vector<GnssType> &type) const;
  void applyTecBias(Double tecBias);
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const GnssSignalBias &x);
template<> void load(InArchive  &ar, GnssSignalBias &x);

/** @brief Write into a GnssSignalBias file. */
void writeFileGnssSignalBias(const FileName &fileName, const GnssSignalBias &x);

/** @brief Read from a GnssSignalBias file. */
void readFileGnssSignalBias(const FileName &fileName, GnssSignalBias &x);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
