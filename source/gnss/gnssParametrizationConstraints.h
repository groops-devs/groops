/***********************************************/
/**
* @file gnssParametrizationConstraints.h
*
* @brief GNSS parameter constraints.
*
* @author Torsten Mayer-Guerr
* @date 2019-05-28
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTS__
#define __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationConstraints
static const char *docstringGnssParametrizationConstraints = R"(
\section{GnssParametrizationConstraints}\label{gnssParametrizationConstraintsType}

Constraints for selected \configClass{parameters}{parameterSelectorType} can be added to the normal equation system
in case there are no direct constraint options in the respective config element (e.g. \config{noNetRotationSigma} in
\configClass{gnssParametrizationReceiver:stationNetwork}{gnssParametrizationReceiverType:stationNetwork}).

For details on possible constraints see \configClass{gnssParametrizationReceiver:stationNetwork}{gnssParametrizationReceiverType:stationNetwork}
(tropospheric parameters), \configClass{gnssParametrizationTransmitter}{gnssParametrizationTransmitterType} (antenna center offsets),
or \configClass{gnssParametrizationEarthRotation}{gnssParametrizationEarthRotationType} (dUT1).

See also \program{GnssProcessing}.
)";
#endif

/***********************************************/

#include "classes/parameterSelector/parameterSelector.h"
#include "gnss/gnss.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssParametrizationConstraints;
typedef std::shared_ptr<GnssParametrizationConstraints> GnssParametrizationConstraintsPtr;

/***** CLASS ***********************************/

/** @brief GNSS parameter constraints.
* An Instance of this class can be created by @ref readConfig. */
class GnssParametrizationConstraints : public Gnss::Parametrization
{
  ParameterSelectorPtr parameterSelector;
  Double      sigma;
  Double      bias;
  Bool        relativeToApriori;
  std::string comment;

public:
  GnssParametrizationConstraints(Config &config, const std::string &name);
  virtual ~GnssParametrizationConstraints() {}

  void observationEquation(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;

  /** @brief creates an derived instance of this class. */
  static GnssParametrizationConstraintsPtr create(Config &config, const std::string &name) {return std::make_shared<GnssParametrizationConstraints>(config, name);}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssParametrizationConstraints.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssParametrizationConstraints */
template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationConstraintsPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
