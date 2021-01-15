/***********************************************/
/**
* @file gnssParametrizationAntenna.h
*
* @brief GNSS antenna center variations.
*
* @author Torsten Mayer-Guerr
* @date 2019-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONANTENNA__
#define __GROOPS_GNSSPARAMETRIZATIONANTENNA__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationAntenna
static const char *docstringGnssParametrizationAntenna = R"(
\section{GnssParametrizationAntenna}\label{gnssParametrizationAntennaType}

Antenna center offsets or variations for satellite or ground station antennas can be estimated by defining
the parametrization via \configClass{antennaCenterVariations}{parametrizationGnssAntennaType}.
The defined parametrization is set up for each type listed in \configClass{patternTypes}{gnssType}, so,
for example, GPS phase center offsets can be estimated separately for \verb|L1*G|, \verb|L2*G|, and \verb|L5*G|.

If \config{groupAntennas}=\verb|yes| parameters are set up per antenna type instead of per individual antenna.
This enables the estimation of, for example, one set of antenna center variations per GPS satellite block.

See also \program{GnssProcessing}.
)";
#endif

/***********************************************/

#include "gnss/gnss.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssStationInfo;
class GnssParametrizationAntenna;
typedef std::shared_ptr<GnssParametrizationAntenna> GnssParametrizationAntennaPtr;

/***** CLASS ***********************************/

/** @brief GNSS antenna center variations.
* An Instance of this class can be created by @ref readConfig. */
class GnssParametrizationAntenna : public Gnss::Parametrization
{
  std::vector<GnssType>              typesPattern;
  Bool                               addNonMatchingTypes;
  Bool                               ignoreSerial;

  Bool                               deviceIsReceiver;
  std::vector<UInt>                  device2antenna;

  std::vector<std::string>           antennaNames;
  std::vector<std::vector<Gnss::ParameterIndex>> indexParameter; // for each antenna and pattern
  std::vector<std::vector<GnssType>> types;          // for each antenna and pattern
  std::vector<std::vector<Bool>>     usedTypes;      // for each antenna and pattern
  ParametrizationGnssAntennaPtr      parametrization;

public:
  /// Constructor.
  GnssParametrizationAntenna(Config &config, const std::string &name);

  /// Destructor.
  virtual ~GnssParametrizationAntenna() {}

  // must be called from all nodes
  void addAntennas(UInt idDevice, Bool deviceIsReceiver, const GnssStationInfo &stationInfo, const std::vector<GnssType> &types);

  // Realization of Gnss::Parametrization
  // ------------------------------------
  void   initParameter(Gnss::NormalEquationInfo &normalEquationInfo) override;
  Bool   isDesignMatrix(const Gnss::NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch) const override;
  void   designMatrix(const Gnss::NormalEquationInfo &normalEquationInfo, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const override;

  /** @brief creates an derived instance of this class. */
  static GnssParametrizationAntennaPtr create(Config &config, const std::string &name) {return std::make_shared<GnssParametrizationAntenna>(config, name);}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssParametrizationAntenna.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssParametrizationAntenna */
template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationAntennaPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
