/***********************************************/
/**
* @file plotMisc.h
*
* @brief Miscellaneous classes for plots.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2010-01-13
*
*/
/***********************************************/

#ifndef __GROOPS_PLOTMISC__
#define __GROOPS_PLOTMISC__

#include "base/import.h"
#include "config/config.h"

/** @addtogroup plotGroup */
/// @{

/***** CLASS ***********************************/

// Latex documentation
#ifdef DOCSTRING_PlotColor
static const char *docstringPlotColor = R"(
\section{PlotColor}\label{plotColorType}
Selects a color.
Used in \program{PlotDegreeAmplitudes}, \program{PlotGraph}, \program{PlotMap},
\program{PlotMatrix}, \program{PlotSphericalHarmonicsTriangle}.

)";
#endif

class PlotColor;
typedef std::shared_ptr<PlotColor> PlotColorPtr;

/** @brief Color options for GMT plots.
* An Instance of this class can be created by @ref readConfig. */
class PlotColor
{
  std::string colorStr; //!< color string needed by GMT

public:
  /// Constructor.
  PlotColor(Config &config, const std::string &name);

  std::string str() const {return colorStr;}

  /** @brief creates an derived instance of this class. */
  static PlotColorPtr create(Config &config, const std::string &name) {return PlotColorPtr(new PlotColor(config, name));}
};

/** @brief Creates a class PlotColor.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a plotColor is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] plotColor Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlotColor */
template<> Bool readConfig(Config &config, const std::string &name, PlotColorPtr &plotColor, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***** CLASS ***********************************/

// Latex documentation
#ifdef DOCSTRING_PlotLine
static const char *docstringPlotLine = R"(
\section{PlotLine}\label{plotLineType}
Defines the line style to be plotted.
)";
#endif

class PlotLine;
typedef std::shared_ptr<PlotLine> PlotLinePtr;

/** @brief Line style for plots
* An Instance of this class can be created by @ref readConfig. */
class PlotLine
{
protected:
  std::string  style;
  Double       width;
  PlotColorPtr color;

public:
  virtual ~PlotLine() {}

  virtual std::string str() const;

  /** @brief creates an derived instance of this class. */
  static PlotLinePtr create(Config &config, const std::string &name);
};

/** @brief Creates a class PlotLine.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a PlotLine is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] PlotLine Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlotLine */
template<> Bool readConfig(Config &config, const std::string &name, PlotLinePtr &PlotLine, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***** CLASS ***********************************/

// Latex documentation
#ifdef DOCSTRING_PlotSymbol
static const char *docstringPlotSymbol = R"(
\section{PlotSymbol}\label{plotSymbolType}
Plots a symbol as used e.g. in \configClass{plotGraphLayer:linesAndPoints}{plotGraphLayerType:linesAndPoints}
or \configClass{plotMapLayer:points}{plotMapLayerType:points}.
)";
#endif

class PlotSymbol;
typedef std::shared_ptr<PlotSymbol> PlotSymbolPtr;

/** @brief Symbol options for GMT plots.
* An Instance of this class can be created by @ref readConfig. */
class PlotSymbol
{
  Bool        colorFromValue_;
  std::string scriptLine;
  std::string legendLine;

public:
  /// Constructor.
  PlotSymbol(Config &config, const std::string &name);

  /** @brief Return whether symbol color is determined externally. */
  Bool requiresColorBar() const {return colorFromValue_;}

  /** @brief Script entry. */
  std::string str() const {return scriptLine;}

  std::string legendStr() const {return legendLine;}

  /** @brief creates an derived instance of this class. */
  static PlotSymbolPtr create(Config &config, const std::string &name) {return PlotSymbolPtr(new PlotSymbol(config, name));}
};

/** @brief Creates a class PlotSymbol.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a plotSymbol is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] plotSymbol Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlotSymbol */
template<> Bool readConfig(Config &config, const std::string &name, PlotSymbolPtr &plotSymbol, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***** CLASS ***********************************/

/** @brief Basic options for GMT plots. */
class PlotBasics
{
public:
  std::string nameProgram;
  FileName    fileNamePlot;
  FileName    baseDirectory, workingDirectory;
  Double      width, height;
  Bool        drawGridOnTop;
  std::string title;
  Double      marginTitle;
  UInt        titleSize;
  Bool        transparent;
  UInt        dpi;
  Bool        removeFiles;
  Bool        view;
  std::string optionsString;

  void read(Config &config, const std::string &nameProgram_, const FileName &fileNamePlot_, const std::string &title_,
            const std::string &defaultWidth, const std::string &defaultHeight, const std::vector<std::string> &optionsDefault=std::vector<std::string>());

  FileName    fileNameScript() const;
  std::string scriptHeader() const;
  std::string scriptTrailer() const;
  Bool        runScript() const;

  /** @brief Axis anntotation in GMT plots. */
  static std::string axisTicks(Bool isLog, Double vmin, Double vmax, Double annotation, Double frame, Double grid, const std::string &unit="", const std::string &label="");

  static std::string scriptSetVariable(const std::string &var, const std::string &value);
  static std::string scriptVariable(const std::string &var);
  static std::string scriptError2Null();

  static UInt gmtVersion();
};

/***********************************************/

/// @}

/***********************************************/

#endif /* __GROOPS_PLOTMISC__ */
