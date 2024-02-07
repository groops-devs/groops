/***********************************************/
/**
* @file noiseGenerator.h
*
* @brief Generate different types of noise.
*
* @author Matthias Ellmer
* @date 2013-09-18
*
*/
/***********************************************/

#ifndef __GROOPS_NOISEGENERATOR__
#define __GROOPS_NOISEGENERATOR__

// Latex documentation
#ifdef DOCSTRING_NoiseGenerator
static const char *docstringNoiseGenerator = R"(
\section{NoiseGenerator}\label{noiseGeneratorType}
This class implements the generation of different types of noise.
It provides a generic interface that can be implemented by different
types of generators. The characteristics of the generated noise
is determined by the generators. See the appropriate documentation
for more information.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup noiseGeneratorGroup NoiseGenerator
* @brief Generate different types of noise.
* @ingroup classesGroup
* The interface is given by @ref NoiseGenerator. */
/// @{

/***** TYPES ***********************************/

class NoiseGenerator;
class NoiseGeneratorBase;
typedef std::shared_ptr<NoiseGenerator> NoiseGeneratorPtr;

/***** CLASS ***********************************/

/** @brief Generate different types of noise.
* An Instance of this class can be created by @ref readConfig. */
class NoiseGenerator
{
public:
/// Constructor
NoiseGenerator(Config &config, const std::string &name);

/// Destructor.
~NoiseGenerator();

/** @brief Generalized Matrix containing one or more columns of noise.
* The type of noise, as well as details like correlations etc. are up to
* the specific implementation.
* @param samples number of noise samples to be returned (rows).
* @param series  number of noise series to be returned (columns)*/
Matrix noise(UInt samples, UInt series) const;

/** @brief Realization of noise vector.
* The type of noise, as well as details like correlations etc. are up to
* the specific implementation.
* @param samples number of noise samples to be returned. */
Vector noise(UInt samples) const {return noise(samples, 1);}

/** @brief Covariance function of generated noise
 * @param length of the covariance function.
 * @param sampling of the time series in seconds.
 * @return @p length x 2 Matrix with time lag in first column and autocovariance in second */
Matrix covarianceFunction(UInt length, Double sampling = 1) const;

/** @brief creates an derived instance of this class. */
static NoiseGeneratorPtr create(Config &config, const std::string &name) {return NoiseGeneratorPtr(new NoiseGenerator(config, name));}

private:
  std::vector<NoiseGeneratorBase*> noiseGenerator;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class NoiseGenerator.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a noiseGenerator with zero noise is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] noiseGenerator Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates NoiseGenerator */
template<> Bool readConfig(Config &config, const std::string &name, NoiseGeneratorPtr &noiseGenerator, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class NoiseGeneratorBase
{
public:
  virtual ~NoiseGeneratorBase() {}
  virtual Matrix noise(UInt samples, UInt series) = 0;
  virtual Vector covarianceFunction(UInt length, Double sampling = 1) = 0;
};

/***********************************************/

#endif
