/***********************************************/
/**
* @file plotLegend.h
*
* @brief Legend.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2020-05-03
*
*/
/***********************************************/

#ifndef __GROOPS_PLOTLEGEND__
#define __GROOPS_PLOTLEGEND__

// Latex documentation
#ifdef DOCSTRING_PlotLegend
static const char *docstringPlotLegend = R"(
\section{PlotLegend}\label{plotLegendType}
Plot a legend of the descriptions provided in
\configClass{plotGraphLayer}{plotGraphLayerType} in \program{PlotGraph}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "plot/plotMisc.h"

/** @addtogroup plotGroup */
/// @{

/***** TYPES ***********************************/

class PlotLegend;
typedef std::shared_ptr<PlotLegend> PlotLegendPtr;

/***** CLASS ***********************************/

/** @brief A Legend for GMT plots.
* An Instance of this class can be created by @ref readConfig. */
class PlotLegend
{
  Double       positionX, positionY;
  UInt         ncolumns;
  Double       width, height;
  PlotLinePtr  edgeLine;
  PlotColorPtr textColor, fillColor;
  std::string  anchor;
  Bool         _hasEntries;

public:
  /** @brief Constructor. */
  PlotLegend(Config &config, const std::string &name);

  std::string scriptEntry() const;

  void writeDataFile(const FileName &workingDirectory, const std::vector<std::string> &legendText);

  /** @brief creates an derived instance of this class. */
  static PlotLegendPtr create(Config &config, const std::string &name) {return PlotLegendPtr(new PlotLegend(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates a class PlotLegend.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a plotLegend is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] plotLegend Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlotLegend */
template<> Bool readConfig(Config &config, const std::string &name, PlotLegendPtr &plotLegend, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif
