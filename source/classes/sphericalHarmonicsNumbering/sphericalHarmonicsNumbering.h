/***********************************************/
/**
* @file sphericalHarmonicsNumbering.h
*
* @brief Numbering schema of spherical harmonics coefficients.
*
* @author Torsten Mayer-Guerr
* @date 2009-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICSNUMBERING__
#define __GROOPS_SPHERICALHARMONICSNUMBERING__

// Latex documentation
#ifdef DOCSTRING_SphericalHarmonicsNumbering
static const char *docstringSphericalHarmonicsNumbering = R"(
\section{SphericalHarmonicsNumbering}\label{sphericalHarmonicsNumberingType}
This class organizes the numbering scheme of spherical harmonics coefficients
in a parameter vector (e.g \program{Gravityfield2SphericalHarmonicsVector} and the design matrix of
\configClass{parametrizationGravity:sphericalHarmoncis}{parametrizationGravityType:sphericalHarmonics}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup sphericalHarmonicsNumberingGroup SphericalHarmonicsNumbering
* @brief Numbering schema of spherical harmonics coefficients.
* @ingroup classesGroup
* The interface is given by @ref SphericalHarmonicsNumbering. */
/// @{

/***** TYPES ***********************************/

class SphericalHarmonicsNumbering;
typedef std::shared_ptr<SphericalHarmonicsNumbering> SphericalHarmonicsNumberingPtr;

/***** CLASS ***********************************/

/** @brief Numbering schema of spherical harmonics coefficients.
* An Instance of this class can be created by @ref readConfig. */
class SphericalHarmonicsNumbering
{
public:
  /// Destructor.
  virtual ~SphericalHarmonicsNumbering() {}

  /** @brief Total number of coefficients. */
  virtual UInt parameterCount(UInt maxDegree, UInt minDegree) const = 0;

  /** @brief Numbering scheme to interpret values in x-vector.
  * @param maxDegree input
  * @param minDegree input
  * @param[out] degree degree n of the value in x(idx) at index idx.
  * @param[out] order  order m of the value in x(idx) at index idx.
  * @param[out] cs     is the value in x(idx) at index idx Cnm (cs[idx]==0) or Snm (cs[idx]==1). */
  void numbering(UInt maxDegree, UInt minDegree, std::vector<UInt> &degree, std::vector<UInt> &order, std::vector<UInt> &cs) const;

  /** @brief Numbering scheme to find coefficient in x-vector.
  * Coefficients not included in the x vector are indicated with NULLINDEX.
  * @param maxDegree input
  * @param minDegree input
  * @param[out] Cnm index in x vector for Cnm(n,m).
  * @param[out] Snm index in x vector for Cnm(n,m). */
  virtual void numbering(UInt maxDegree, UInt minDegree, std::vector<std::vector<UInt>> &Cnm, std::vector<std::vector<UInt>> &Snm) const = 0;

  /** @brief Creates an derived instance of this class. */
  static SphericalHarmonicsNumberingPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class SphericalHarmonicsNumbering.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a numbering remains untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] numbering Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates SphericalHarmonicsNumbering */
template<> Bool readConfig(Config &config, const std::string &name, SphericalHarmonicsNumberingPtr &numbering,
                           Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif
