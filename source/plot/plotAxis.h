/***********************************************/
/**
* @file plotAxis.h
*
* @brief Axis ticks, limits and labels.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2016-07-23
*
*/
/***********************************************/

#ifndef __GROOPS_PLOTAXIS__
#define __GROOPS_PLOTAXIS__

// Latex documentation
#ifdef DOCSTRING_PlotAxis
static const char *docstringPlotAxis = R"(
\section{PlotAxis}\label{plotAxisType}
Defines the style of the axes of \program{PlotGraph}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/** @addtogroup plotGroup */
/// @{

/***** TYPES ***********************************/

class PlotAxis;
typedef std::shared_ptr<PlotAxis> PlotAxisPtr;

/***** CLASS ***********************************/

/** @brief Axis ticks, limits and labels.
* An Instance of this class can be created by @ref readConfig. */
class PlotAxis
{
protected:
  Double       vMin;
  Double       vMax;
  PlotColorPtr color;
  PlotLinePtr  gridLine;
  Double       margin;          // margin between label and ticks
  Bool         changeDirection; // flag whether to reverse the axis

public:
  virtual ~PlotAxis() {}

  virtual void setAutoInterval(Double minAuto, Double maxAuto) = 0;

  /** @brief Get minimum value for this axis. */
  virtual Double getMin() const {return vMin;}

  /** @brief Get maximum value for this axis. */
  virtual Double getMax() const {return vMax;}

  /** @brief return margin of this axis
  * @return margin height/width in cm */
  virtual Double getMargin() const {return margin;}

  /** @brief direction of the axis (+1 or -1) */
  virtual int direction() const {return changeDirection ? -1 : +1;}

  /** @brief return if axis is in logarithmic spacing */
  virtual Bool isLogarithmic() const {return FALSE;}

  /** @brief GMT axis modifier (i.e. logarithmic, normal or time axis)
  * @return modifier string representation of axis type */
  virtual std::string axisModifier() const = 0;

  /** @brief Create temporary files for axis plots. */
  virtual void writeDataFile(const FileName &/*workingDirectory*/, const std::string &/*axis*/) {}

  /** @brief GMT script line for axis object
  * @return line GMT string representation of axis object */
  virtual std::string scriptEntry(const std::string &axis, Bool withGrid=TRUE) const = 0;

  /** @brief creates an derived instance of this class. */
  static PlotAxisPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates a class PlotAxis.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a plotAxis is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] plotAxis Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlotAxis */
template<> Bool readConfig(Config &config, const std::string &name, PlotAxisPtr &plotAxis, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
