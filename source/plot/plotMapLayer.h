/***********************************************/
/**
* @file plotMapLayer.h
*
* @brief plot layers of maps.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2015-10-23
*
*/
/***********************************************/

#ifndef __GROOPS_PLOTMAPLAYER__
#define __GROOPS_PLOTMAPLAYER__

#include "base/import.h"
#include "config/config.h"

/** @addtogroup plotGroup */
/// @{

/***** CLASS ***********************************/

// Latex documentation
#ifdef DOCSTRING_PlotMapLayer
static const char *docstringPlotMapLayer = R"(
\section{PlotMapLayer}\label{plotMapLayerType}
Defines the content of a map of \program{PlotMap}. Multiple layers are are plotted sequentially.
)";
#endif

/***** TYPES ***********************************/

class PlotMapLayer;
typedef std::shared_ptr<PlotMapLayer> PlotMapLayerPtr;

/***** CLASS ***********************************/

/** @brief Layers for maps.
* An Instance of this class can be created by @ref readConfig. */
class PlotMapLayer
{
protected:
  FileName              dataFileName;
  std::vector<Vector3d> points;
  std::vector<Double>   areas;
  Matrix                data;
  Angle                 bufferLon, bufferLat;

public:
  virtual ~PlotMapLayer() {}

  virtual Bool requiresColorBar() const {return FALSE;}
  virtual void boundary(const Ellipsoid &ellipsoid, Angle &minL, Angle &maxL, Angle &minB, Angle &maxB) const;
  virtual void getIntervalZ(Bool isLogarithmic, Double &minZ, Double &maxZ) const;
  virtual void writeDataFile(const Ellipsoid &ellipsoid, const FileName &workingDirectory, UInt idxLayer);
  virtual std::string scriptStatisticsInfo(UInt fontSize, Double width, const FileName &workingDirectory, UInt idxLayer) const;
  virtual std::string scriptEntry() const = 0;
  virtual std::string legendEntry(const FileName &workingDirectory, UInt idxLayer) const;

  /** @brief creates an derived instance of this class. */
  static PlotMapLayerPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates a class PlotMapLayer.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a plotMapLayer is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] plotMapLayer Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlotMapLayer */
template<> Bool readConfig(Config &config, const std::string &name, PlotMapLayerPtr &plotMapLayer, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/


#endif /* __GROOPS__ */
