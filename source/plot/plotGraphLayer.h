/***********************************************/
/**
* @file plotGraphLayer.h
*
* @brief Lines, points in 2d plots.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2016-07-23
*
*/
/***********************************************/

#ifndef __GROOPS_PLOTGRAPHLAYER__
#define __GROOPS_PLOTGRAPHLAYER__

// Latex documentation
#ifdef DOCSTRING_PlotGraphLayer
static const char *docstringPlotGraphLayer = R"(
\section{PlotGraphLayer}\label{plotGraphLayerType}
Defines the content of an xy-plot of \program{PlotGraph}.
Multiple layers are are plotted sequentially. With \config{plotOnSecondAxis}
the alternative y-axis on the right hand side can be selected if provided.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/** @addtogroup plotGroup */
/// @{

/***** TYPES ***********************************/

class PlotGraphLayer;
typedef std::shared_ptr<PlotGraphLayer> PlotGraphLayerPtr;

/***** CLASS ***********************************/

/** @brief Lines, points in 2d plots.
* An Instance of this class can be created by @ref readConfig. */
class PlotGraphLayer
{
protected:
  FileName dataFileName;
  Matrix   data;
  Bool     onSecondAxis;

  virtual Double bufferX() const {return 0;}
  virtual Double bufferY() const {return 0;}

public:
  virtual ~PlotGraphLayer() {}

  Bool drawOnSecondAxis() const {return onSecondAxis;}
  virtual Bool requiresColorBar() const {return FALSE;}
  virtual void getIntervalX(Bool isLogarithmic, Double &minX, Double &maxX) const;
  virtual void getIntervalY(Bool isLogarithmic, Double minX, Double maxX, Double &minY, Double &maxY) const;
  virtual void getIntervalZ(Bool isLogarithmic, Double minX, Double maxX, Double minY, Double maxY, Double &minZ, Double &maxZ) const;
  virtual void writeDataFile(const FileName &workingDirectory, UInt idxLayer, Double minX, Double maxX, Double minY, Double maxY);
  virtual std::string scriptEntry() const = 0;
  virtual std::string legendEntry() const {return std::string();}

  /** @brief creates an derived instance of this class. */
  static PlotGraphLayerPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates a class PlotGraphLayer.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a plotGraphLayer is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] plotGraphLayer Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlotGraphLayer */
template<> Bool readConfig(Config &config, const std::string &name, PlotGraphLayerPtr &plotGraphLayer, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
