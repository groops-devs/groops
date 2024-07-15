/***********************************************/
/**
* @file slrStationGenerator.h
*
* @brief Provides a list of stations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRRECEIVERGENERATOR__
#define __GROOPS_SLRRECEIVERGENERATOR__

// Latex documentation
#ifdef DOCSTRING_SlrStationGenerator
static const char *docstringSlrStationGenerator = R"(
\section{SlrStationGenerator}\label{slrStationGeneratorType}
Definition and basic information of SLR ground stations.

See also \program{SlrProcessing}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "slr/slrStation.h"
#include "slr/slrSatellite.h"

/**
* @defgroup slrStationGeneratorGroup SlrStationGenerator
* @brief Provides a list of stations.
* @ingroup classesGroup
* The interface is given by @ref SlrStationGenerator.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Slr;
class SlrStationGenerator;
class SlrStationGeneratorBase;
typedef std::shared_ptr<SlrStationGenerator> SlrStationGeneratorPtr;

/***** CLASS ***********************************/

/** @brief Provides a list of stations.
* An Instance of this class can be created by @ref readConfig. */
class SlrStationGenerator
{
  std::vector<SlrStationGeneratorBase*> base;

public:
  /** @brief Constructor from config. */
  SlrStationGenerator(Config &config, const std::string &name);

  /// Destructor.
 ~SlrStationGenerator();

  /** @brief Iniatialize and returns a vector of stations. */
  std::vector<SlrStationPtr> stations(const std::vector<Time> &times, const std::vector<SlrSatellitePtr> &satellites,
                                      EarthRotationPtr earthRotation);

  /** @brief preprocess the observations of stations. */
  void preprocessing(Slr *slr);

//   /** @brief simulate the observations of stations. */
//   void simulation(NoiseGeneratorPtr noiseObs, Slr *slr);

  /** @brief creates an derived instance of this class. */
  static SlrStationGeneratorPtr create(Config &config, const std::string &name) {return SlrStationGeneratorPtr(new SlrStationGenerator(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class SlrStationGenerator.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates SlrStationGenerator */
template<> Bool readConfig(Config &config, const std::string &name, SlrStationGeneratorPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class SlrStationGeneratorBase
{
public:
  virtual ~SlrStationGeneratorBase() {}

  virtual void init(const std::vector<Time> &times, const std::vector<SlrSatellitePtr> &satellites,
                    EarthRotationPtr earthRotation, std::vector<SlrStationPtr> &stations) = 0;

//   virtual void preprocessing(Slr *slr, Parallel::CommunicatorPtr comm) = 0;
//
//   virtual void simulation(NoiseGeneratorPtr noiseObs, Slr *slr, Parallel::CommunicatorPtr comm) = 0;

  static void printPreprocessingInfos(const std::string &header, const std::vector<SlrStationPtr> &stations,
                                      Bool disabledOnly);
};

/***********************************************/

#endif /* __GROOPS___ */
