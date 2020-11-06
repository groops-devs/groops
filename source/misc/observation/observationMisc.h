/***********************************************/
/**
* @file observationMisc.h
*
* @brief Right hand sides.
*
* @author Torsten Mayer-Guerr
* @date 2008-07-28
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONMISC__
#define __GROOPS_OBSERVATIONMISC__

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "classes/forces/forces.h"
#include "classes/tides/tides.h"
#include "classes/gravityfield/gravityfield.h"

/** @addtogroup miscGroup */
/// @{

/***** CLASS ***********************************/

// Latex documentation
#ifdef DOCSTRING_PodRightSide
static const char *docstringPodRightSide = R"(
\section{PodRightSide}\label{podRightSideType}
Observation vector for precise orbit data (POD) of \configClass{observation}{observationType}
equations in a least squares adjustment. The observations are reduced by the effect of
\configFile{inputfileAccelerometer}{instrument} and \configClass{forces}{forcesType}
(observed minus computed).
)";
#endif

class PodRightSide;
typedef std::shared_ptr<PodRightSide> PodRightSidePtr;

/** @brief Input for observation vectors (precise orbit data, POD). */
class PodRightSide
{
public:
  InstrumentFilePtr orbitFile;
  InstrumentFilePtr accelerometerFile;
  ForcesPtr         forces;

  /// Constructor
  PodRightSide(Config &config, const std::string &name);

  /** @brief creates an derived instance of this class. */
  static PodRightSidePtr create(Config &config, const std::string &name) {return PodRightSidePtr(new PodRightSide(config, name));}
};

/** @brief Creates a class PodRightSide.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PodRightSide */
template<> Bool readConfig(Config &config, const std::string &name, PodRightSidePtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***** CLASS ***********************************/

// Latex documentation
#ifdef DOCSTRING_SstRightSide
static const char *docstringSstRightSide = R"(
\section{SstRightSide}\label{sstRightSideType}
Observation vector for GRACE like data (satellite-tracking and precise orbit data (POD))
of \configClass{observation}{observationType} equations in a least squares adjustment.
The observations are reduced by the effect of \configFile{inputfileAccelerometer}{instrument}
and \configClass{forces}{forcesType} (observed minus computed).
)";
#endif

class SstRightSide;
typedef std::shared_ptr<SstRightSide> SstRightSidePtr;

/** @brief Input for observation vectors (GRACE POD and KBR data). */
class SstRightSide
{
public:
  std::vector<InstrumentFilePtr> sstFile;
  InstrumentFilePtr              orbit1File, orbit2File;
  InstrumentFilePtr              accelerometer1File, accelerometer2File;
  ForcesPtr                      forces;

  /// Constructor
  SstRightSide(Config &config, const std::string &name);

  /** @brief creates an derived instance of this class. */
  static SstRightSidePtr create(Config &config, const std::string &name) {return SstRightSidePtr(new SstRightSide(config, name));}
};

/** @brief Creates a class SstRightSide.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates SstRightSide */
template<> Bool readConfig(Config &config, const std::string &name, SstRightSidePtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***** CLASS ***********************************/

// Latex documentation
#ifdef DOCSTRING_SggRightSide
static const char *docstringSggRightSide = R"(
\section{SggRightSide}\label{sggRightSideType}
Observation vector for gradiometer data (satellite gravity gradiometry, SGG)
of \configClass{observation}{observationType} equations in a least squares adjustment.
The observations are reduced by an \configFile{inputfileReferenceGradiometer}{instrument},
the effect of \configClass{referencefield}{gravityfieldType}, and \configClass{tides}{tidesType}
(observed minus computed).

The reference gradiometer data can be precomputed with \program{SimulateGradiometer}.
)";
#endif

class SggRightSide;
typedef std::shared_ptr<SggRightSide> SggRightSidePtr;

/** @brief Input for observation vectors (GOCE SGG). */
class SggRightSide
{
public:
  InstrumentFilePtr              gradiometerFile;
  std::vector<InstrumentFilePtr> referenceFile;
  GravityfieldPtr                referencefield;
  TidesPtr                       tides;

  /// Constructor
  SggRightSide(Config &config, const std::string &name);

  /** @brief creates an derived instance of this class. */
  static SggRightSidePtr create(Config &config, const std::string &name) {return SggRightSidePtr(new SggRightSide(config, name));}
};

/** @brief Creates a class SggRightSide.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates SggRightSide */
template<> Bool readConfig(Config &config, const std::string &name, SggRightSidePtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS_OBSERVATIONMISC__ */
