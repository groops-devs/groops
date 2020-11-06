/***********************************************/
/**
* @file plotColorbar.h
*
* @brief Colorbar.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2020-05-03
*
*/
/***********************************************/

#ifndef __GROOPS_PLOTCOLORBAR__
#define __GROOPS_PLOTCOLORBAR__

// Latex documentation
#ifdef DOCSTRING_PlotColorbar
static const char *docstringPlotColorbar = R"(
\section{PlotColorbar}\label{plotColorbarType}
A colorbar as used in \program{PlotMap}, \program{PlotMatrix}, \program{PlotSphericalHarmonicsTriangle}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/** @addtogroup plotGroup */
/// @{

/***** TYPES ***********************************/

class PlotColorbar;
typedef std::shared_ptr<PlotColorbar> PlotColorbarPtr;

/***** CLASS ***********************************/

/** @brief A Colorbar for GMT plots.
* An Instance of this class can be created by @ref readConfig. */
class PlotColorbar
{
  Double      vMin, vMax;
  std::string colorTable;
  Double      annotation;
  std::string label, unit;
  Bool        isLog;
  Bool        triangleLeft, triangleRight;
  Bool        illuminate;
  Bool        vertical;
  Double      margin;
  Double      length;
  Bool        isPlot;

public:
  /** @brief Constructor. */
  PlotColorbar(Config &config, const std::string &name);

  void setAutoInterval(Double minAuto, Double maxAuto);

  /** @brief Get minimum value for this axis. */
  Double getMin() const {return vMin;}

  /** @brief Get maximum value for this axis. */
  Double getMax() const {return vMax;}

  Bool isLogarithmic() const {return isLog;}

  std::string scriptColorTable() const;
  std::string scriptEntry(Double width, Double height, Double marginX, Double marginY) const;

  /** @brief creates an derived instance of this class. */
  static PlotColorbarPtr create(Config &config, const std::string &name) {return PlotColorbarPtr(new PlotColorbar(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates a class PlotColorbar.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a plotColorbar is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] plotColorbar Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlotColorbar */
template<> Bool readConfig(Config &config, const std::string &name, PlotColorbarPtr &plotColorbar, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif
