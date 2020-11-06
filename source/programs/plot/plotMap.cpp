/***********************************************/
/**
* @file plotMap.cpp
*
* @brief Plot maps.
*
* @author Enrico Kurtenbach
* @author Torsten Mayer-Guerr
* @author Christian Pock
* @author Andreas Kvas
* @date 2015-10-20
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Generates a map using the GMT Generic Mapping Tools (\url{https://www.generic-mapping-tools.org}).
A variety of image file formats are supported (e.g. png, jpg, eps) determined by the extension of \config{outputfile}.

The base map is defined by a \configClass{projection}{plotMapProjectionType} of an ellipsoid (\config{R}, \config{inverseFlattening}).
The content of the map itself is defined by one or more \configClass{layer}{plotMapLayerType}s.

The plot programs create a temporary directory in the path of \config{outputfile}, writes all needed data into it,
generates a batch/shell script with the GMT commands, execute it, and remove the temporary directory.
With setting \config{options:removeFiles}=false the last step is skipped and it is possible to adjust the plot manually
to specific publication needs. Individual GMT settings are adjusted with \config{options:options}="\verb|FORMAT=value|",
see \url{https://docs.generic-mapping-tools.org/latest/gmt.conf.html}.

See also: \program{PlotDegreeAmplitudes}, \program{PlotGraph}, \program{PlotMatrix}, \program{PlotSphericalHarmonicsTriangle}.

\fig{!hb}{0.8}{plotMap}{fig:plotMap}{A Robinson projection with griddedData (geoid), coast, polygon (amazon), and points (IGS stations) layer.}
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "plot/plotColorbar.h"
#include "plot/plotMapLayer.h"
#include "plot/plotMapProjection.h"
#include "plot/plotMisc.h"

/***********************************************/

/** @brief Plot maps.
* @ingroup programsGroup */
class PlotMap
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(PlotMap, SINGLEPROCESS, "Plot maps", Plot, Grid)

/***********************************************/

void PlotMap::run(Config &config)
{
  try
  {
    FileName                     fileNamePlot;
    std::string                  title;
    Bool                         statisticInfos;
    std::vector<PlotMapLayerPtr> layer;
    Double                       a, f;
    Angle                        minL(NAN_EXPR), maxL(NAN_EXPR), minB(NAN_EXPR), maxB(NAN_EXPR);
    Angle                        annotation, frame, grid;
    PlotColorbarPtr              colorbar;
    PlotMapProjectionPtr         projection;
    PlotBasics                   plotBasics;

    renameDeprecatedConfig(config, "annotation", "majorTickSpacing", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "frame",      "minorTickSpacing", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "grid",       "gridLineSpacing",  date2time(2020, 4, 23));

    readConfig(config, "outputfile",        fileNamePlot,   Config::MUSTSET,  "",  "*.png, *.jpg, *.eps, ...");
    readConfig(config, "title",             title,          Config::OPTIONAL, "",  "");
    readConfig(config, "statisticInfos",    statisticInfos, Config::DEFAULT,  "0", "");
    readConfig(config, "layer",             layer,          Config::MUSTSET,  "",  "");
    readConfig(config, "R",                 a,              Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening", f,              Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    readConfig(config, "minLambda",         minL,           Config::OPTIONAL, "-180", "min. longitude (default: compute from input data)");
    readConfig(config, "maxLambda",         maxL,           Config::OPTIONAL, "180",  "max. longitude (default: compute from input data)");
    readConfig(config, "minPhi",            minB,           Config::OPTIONAL, "-90",  "min. latitude (default: compute from input data)");
    readConfig(config, "maxPhi",            maxB,           Config::OPTIONAL, "90",   "max. latitude (default: compute from input data)");
    readConfig(config, "majorTickSpacing",  annotation,     Config::DEFAULT,  "0",    "boundary annotation");
    readConfig(config, "minorTickSpacing",  frame,          Config::DEFAULT,  "0",    "frame tick spacing");
    readConfig(config, "gridLineSpacing",   grid,           Config::DEFAULT,  "0",    "gridline spacing");
    readConfig(config, "colorbar",          colorbar,       Config::OPTIONAL, "1",    "");
    readConfig(config, "projection",        projection,     Config::MUSTSET,  "",     "map projection");
    plotBasics.read(config, "PlotMap", fileNamePlot, title, "12", "");
    if(isCreateSchema(config)) return;

    // ===========================================================

    // calculate boundary
    // ------------------
    Ellipsoid ellipsoid(a, f);
    Angle minL2(+180.*DEG2RAD);
    Angle maxL2(-180.*DEG2RAD);
    Angle minB2(+ 90.*DEG2RAD);
    Angle maxB2(- 90.*DEG2RAD);
    for(UInt i=0; i<layer.size(); i++)
      layer.at(i)->boundary(ellipsoid, minL2, maxL2, minB2, maxB2);
    if(std::isnan(minL)) minL = minL2;
    if(std::isnan(maxL)) maxL = maxL2;
    if(std::isnan(minB)) minB = minB2;
    if(std::isnan(maxB)) maxB = maxB2;

    if(minL == Angle(+180.*DEG2RAD) && (maxL == Angle(-180.*DEG2RAD)))
      std::swap(minL, maxL);
    if(minB > maxB)
      std::swap(minB, maxB);

    // colorbar needed?
    // ----------------
    Bool needColorBar = FALSE;
    for(UInt i=0; i<layer.size(); i++)
      needColorBar = needColorBar || layer.at(i)->requiresColorBar();
    if(needColorBar && !colorbar)
      throw(Exception("Layer requires a color bar"));

    // calculate min and max
    // ---------------------
    if(needColorBar)
    {
      Double vmin =  1e99;
      Double vmax = -1e99;
      for(UInt i=0; i<layer.size(); i++)
        layer.at(i)->getIntervalZ(colorbar->isLogarithmic(), vmin, vmax);
      colorbar->setAutoInterval(vmin, vmax);
    }

    plotBasics.optionsString += " PROJ_ELLIPSOID="+a%"%f"s + f%"/%f"s;
    if(std::isnan(plotBasics.height))
    {
      plotBasics.height = projection->aspectRatio() * plotBasics.width;
      if(projection->aspectRatio() == 0)
        plotBasics.height = (maxB-minB)/(maxL-minL) * plotBasics.width;
    }

    // create data files
    // -----------------
    logStatus<<"create temporary data files"<<Log::endl;
    for(UInt i=0; i<layer.size(); i++)
      layer.at(i)->writeDataFile(ellipsoid, plotBasics.workingDirectory, i);

    // calculate annotation
    // --------------------
    Double margin = (annotation!=0) ? 0.5 : 0.0;

    // create scriptfile
    // -----------------
    logStatus<<"create scriptfile"<<Log::endl;
    {
      OutFile file(plotBasics.fileNameScript());
      file<<plotBasics.scriptHeader();
      if(statisticInfos)
        for(UInt i=0; i<layer.size(); i++)
          file<<layer.at(i)->scriptStatisticsInfo(plotBasics.titleSize, plotBasics.width, plotBasics.workingDirectory, i);
      if(colorbar)
        file<<colorbar->scriptColorTable();
      file<<std::endl;
      file<<PlotBasics::scriptSetVariable("region", RAD2DEG*minL%"%g/"s+RAD2DEG*maxL%"%g/"s+RAD2DEG*minB%"%g/"s+RAD2DEG*maxB%"%g"s)<<std::endl;
      file<<PlotBasics::scriptSetVariable("projection", projection->scriptEntry(plotBasics.width, plotBasics.height).substr(2))<<std::endl;
      file<<"gmt psbasemap -Y"<<-plotBasics.height-margin-plotBasics.marginTitle<<"c "<<projection->scriptEntry(plotBasics.width, plotBasics.height);
      file<<" -R"<<PlotBasics::scriptVariable("region")<<" -Bg"<<grid*RAD2DEG<<"a"<<annotation*RAD2DEG<<"f"<<frame*RAD2DEG<<"WSNE -O -K >> groopsPlot.ps"<<std::endl;
      for(UInt i=0; i<layer.size(); i++)
        file<<layer.at(i)->scriptEntry();
      file<<"gmt psbasemap -J -R -BWSNE -Ba"<<annotation*RAD2DEG<<"f"<<frame*RAD2DEG<<"g"<<(plotBasics.drawGridOnTop ? grid*RAD2DEG : 0.)<<" -O -K >> groopsPlot.ps"<<std::endl;
      if(colorbar)
        file<<colorbar->scriptEntry(plotBasics.width, plotBasics.height, margin, margin);
      for(UInt i=0; i<layer.size(); i++)
        file<<layer.at(i)->legendEntry(plotBasics.workingDirectory, i);
      file<<plotBasics.scriptTrailer();
    }

    // run scriptfile
    // --------------
    logStatus<<"run scriptfile"<<Log::endl;
    plotBasics.runScript();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
