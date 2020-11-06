/***********************************************/
/**
* @file plotMapProjection.h
*
* @brief Map projections used for plotting.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2015-10-23
*
*/
/***********************************************/

#ifndef __GROOPS_PLOTMAPROJECTION__
#define __GROOPS_PLOTMAPROJECTION__

// Latex documentation
#ifdef DOCSTRING_PlotMapProjection
static const char *docstringPlotMapProjection = R"(
\section{PlotMapProjection}\label{plotMapProjectionType}
Selects the underlying projection of \program{PlotMap}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/** @addtogroup plotGroup */
/// @{

/***** TYPES ***********************************/

class PlotMapProjection;
typedef std::shared_ptr<PlotMapProjection> PlotMapProjectionPtr;

/***** CLASS ***********************************/

/** @brief Layers for maps.
* An Instance of this class can be created by @ref readConfig. */
class PlotMapProjection
{
public:
  virtual ~PlotMapProjection() {}

  virtual std::string scriptEntry(Double width, Double height) const = 0;

  virtual Double aspectRatio() const {return 1;}

  /** @brief creates an derived instance of this class. */
  static PlotMapProjectionPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates a class PlotMapProjection.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a plotMapProjection is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] plotMapProjection Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlotMapProjection */
template<> Bool readConfig(Config &config, const std::string &name, PlotMapProjectionPtr &plotMapProjection, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif
