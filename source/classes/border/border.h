/***********************************************/
/**
* @file border.h
*
* @brief Borders of an area on sphere/ellipsoid.
*
* @author Torsten Mayer-Guerr
* @date 2004-10-28
*
*/
/***********************************************/

#ifndef __GROOPS_BORDER__
#define __GROOPS_BORDER__

// Latex documentation
#ifdef DOCSTRING_Border
static const char *docstringBorder = R"(
\section{Border}\label{borderType}
With this class you can select one or more region on the surface of the Earth.
In every instance of Border you can choose whether the specific region is excluded
from the overall result with the switch \config{exclude}.
To determine whether a specific point will be used furthermore the following algorithm will be applied:
In a first step all points are selected if first border excludes points otherwise all points excluded.
When every point will be tested for each instance of border from top to bottom.
If the point is not in the selected region nothing happens.
Otherwise it will included or excluded depending on the switch \config{exclude}.

First Example: The border excludes all continental areas.
The result are points on the oceans only.

Second Example: First border describes the continent north america. The next borders
excludes the great lakes and the last border describes Washington island.
In this configuration points are selected if they are inside north america
but not in the area of the great lakes. But if the point is on Washington island
it will be included again.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup borderGroup Border
* @brief Borders of an area on sphere/ellipsoid.
* @ingroup classesGroup
* The interface is given by @ref Border.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Border;
class BorderBase;
typedef std::shared_ptr<Border> BorderPtr;

/***** CLASS ***********************************/

/** @brief Borders of an area on sphere/ellipsoid.
* An instance of this class can be created with @ref readConfig. */
class Border
{
public:
  /// Constructor.
  Border(Config &config, const std::string &name);

  /// Destructor.
  ~Border();

  /** @brief Is this point inside the area?
  * @param lambda longitude (-PI,PI]
  * @param phi    latitude [PI/2,PI/2] */
  Bool isInnerPoint(Angle lambda, Angle phi) const;

  /** @brief Is this point inside the area?
  * @param point point
  * @param ellipsoid longitude, latitude of point relates to this ellipsoid. */
  Bool isInnerPoint(const Vector3d &point, const Ellipsoid &ellipsoid=Ellipsoid()) const;

  /** @brief creates an derived instance of this class. */
  static BorderPtr create(Config &config, const std::string &name) {return BorderPtr(new Border(config, name));}

private:
  std::vector<BorderBase*> border;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Border.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and a global border is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] border Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Border */
template<> Bool readConfig(Config &config, const std::string &name, BorderPtr &border, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class BorderBase
{
public:
  virtual ~BorderBase() {}

  virtual Bool isInnerPoint(Angle lambda, Angle phi) const = 0;
  virtual Bool isExclude() const = 0;
};

/***********************************************/

#endif /* __GROOPS_BORDER__ */
