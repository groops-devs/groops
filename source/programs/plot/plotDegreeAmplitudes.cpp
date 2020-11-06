/***********************************************/
/**
* @file plotDegreeAmplitudes.cpp
*
* @brief Plot Degree Amplitudes.
*
* @author Enrico Kurtenbach
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2015-10-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Plot degree amplitudes of potential coefficients computed by \program{Gravityfield2DegreeAmplitudes}
or \program{PotentialCoefficients2DegreeAmplitudes} using the GMT Generic Mapping Tools
(\url{https://www.generic-mapping-tools.org}).
A variety of image file formats are supported (e.g. png, jpg, eps) determined by the extension of \config{outputfile}.
This is a convenience program with meaningful default values. The same plots can be generated with the more general \program{PlotGraph}.

\fig{!hb}{0.8}{plotDegreeAmplitudes}{fig:plotDegreeAmplitudes}{Comparison of GRACE solutions (2008-06) with GOCO06s.}
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "plot/plotMisc.h"
#include "plot/plotGraphLayer.h"
#include "plot/plotLegend.h"

/***********************************************/

/** @brief Plot Degree Amplitudes.
* @ingroup programsGroup */
class PlotDegreeAmplitudes
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(PlotDegreeAmplitudes, SINGLEPROCESS, "Plot Degree Amplitudes", Plot, Gravityfield, PotentialCoefficients)

/***********************************************/

void PlotDegreeAmplitudes::run(Config &config)
{
  try
  {
    FileName          fileNamePlot;
    std::string       title;
    std::vector<PlotGraphLayerPtr> layer;
    Double            annotationX=NAN_EXPR, frameX=NAN_EXPR, gridX=NAN_EXPR;
    Double            annotationY=NAN_EXPR, frameY=NAN_EXPR, gridY=NAN_EXPR;
    UInt              minX=INFINITYDEGREE, maxX=INFINITYDEGREE;
    Double            minY=NAN_EXPR, maxY=NAN_EXPR;
    Bool              logX, logY;
    std::string       labelX, labelY, unitY;
    PlotLinePtr       gridLine;
    PlotLegendPtr     legend;
    PlotBasics        plotBasics;

    renameDeprecatedConfig(config, "annotationDegree", "majorTickSpacingDegree", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "frameDegree",      "minorTickSpacingDegree", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "gridDegree",       "gridLineSpacingDegree",  date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "annotationY",      "majorTickSpacingY",      date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "frameY",           "minorTickSpacingY",      date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "gridY",            "gridLineSpacingY",       date2time(2020, 4, 23));

    readConfig(config, "outputfile",             fileNamePlot, Config::MUSTSET,  "",     "*.png, *.jpg, *.eps, ...");
    readConfig(config, "title",                  title,        Config::OPTIONAL, "",     "");
    readConfig(config, "layer",                  layer,        Config::MUSTSET,  "degreeAmplitudes", "");
    readConfig(config, "minDegree",              minX,         Config::OPTIONAL, "0",    "");
    readConfig(config, "maxDegree",              maxX,         Config::OPTIONAL, "",     "");
    readConfig(config, "majorTickSpacingDegree", annotationX,  Config::OPTIONAL, "",     "boundary annotation");
    readConfig(config, "minorTickSpacingDegree", frameX,       Config::OPTIONAL, "",     "frame tick spacing");
    readConfig(config, "gridLineSpacingDegree",  gridX,        Config::OPTIONAL, "",     "gridline spacing");
    readConfig(config, "labelDegree",            labelX,       Config::OPTIONAL, "[degree]", "description of the x-axis");
    readConfig(config, "logarithmicDegree",      logX,         Config::DEFAULT,  "0",    "use logarithmic scale for the x-axis");
    readConfig(config, "minY",                   minY,         Config::OPTIONAL, "1e-6", "");
    readConfig(config, "maxY",                   maxY,         Config::OPTIONAL, "1",    "");
    readConfig(config, "majorTickSpacingY",      annotationY,  Config::OPTIONAL, "",     "boundary annotation");
    readConfig(config, "minorTickSpacingY",      frameY,       Config::OPTIONAL, "",     "frame tick spacing");
    readConfig(config, "gridLineSpacingY",       gridY,        Config::OPTIONAL, "",     "gridline spacing");
    readConfig(config, "unitY",                  unitY,        Config::OPTIONAL, "",     "appended to axis values");
    readConfig(config, "labelY",                 labelY,       Config::OPTIONAL, "geoid height [m]",  "description of the y-axis");
    readConfig(config, "logarithmicY",           logY,         Config::DEFAULT,  "1",    "use logarithmic scale for the y-axis");
    readConfig(config, "gridLine",               gridLine,     Config::OPTIONAL, R"({"solid": {"width":"0.25", "color":"gray"}})", "The style of the grid lines.");
    readConfig(config, "legend",                 legend,       Config::OPTIONAL, "1",    "");
    plotBasics.read(config, "PlotDegreeAmplitudes", fileNamePlot, title, "12", "10");
    if(isCreateSchema(config)) return;

    // ===========================================================

    // calculate min and max
    // ---------------------
    // x-axis
    Double minDegree  =  1e99;
    Double maxDegree  = -1e99;
    for(UInt k = 0; k<layer.size(); k++)
      layer.at(k)->getIntervalX(logX, minDegree, maxDegree);
    if(minX==INFINITYDEGREE) minX = static_cast<UInt>(minDegree);
    if(maxX==INFINITYDEGREE) maxX = static_cast<UInt>(maxDegree);

    // y-axis
    Double minSignal  =  1e99;
    Double maxSignal  = -1e99;
    for(UInt k = 0; k<layer.size(); k++)
      layer.at(k)->getIntervalY(logY, minX, maxX, minSignal, maxSignal);
    if(std::isnan(minY)) minY = minSignal;
    if(std::isnan(maxY)) maxY = maxSignal;

    // create data files
    // -----------------
    logStatus<<"create temporary data files"<<Log::endl;
    for(UInt i=0; i<layer.size(); i++)
      layer.at(i)->writeDataFile(plotBasics.workingDirectory, i, minX, maxX, minY, maxY);

    // create scriptfile
    // -----------------
    logStatus<<"create scriptfile"<<Log::endl;
    {
      OutFile file(plotBasics.fileNameScript());
      file<<plotBasics.scriptHeader();
      file<<PlotBasics::scriptSetVariable("range", minX%"%i/"s+maxX%"%i/"s+minY%"%g/"s+maxY%"%g"s)<<std::endl;
      file<<"gmt psbasemap -Y"<<-plotBasics.height-plotBasics.marginTitle<<"c -R"<<PlotBasics::scriptVariable("range");
      file<<" -JX"<<plotBasics.width<<"c"<<((logX) ? "l" : "")<<"/"<<plotBasics.height<<"c"<<((logY) ? "l" : "");
      file<<" -BWSne";
      file<<" -Bx"<<PlotBasics::axisTicks(logX, minX, maxX, annotationX, frameX, gridX, "",    labelX);
      file<<" -By"<<PlotBasics::axisTicks(logY, minY, maxY, annotationY, frameY, gridY, unitY, labelY);
      if(gridLine)
        file<<" --MAP_GRID_PEN_PRIMARY="<<gridLine->str();
      file<<" -O -K >> groopsPlot.ps"<<std::endl;
      for(UInt i=0; i<layer.size(); i++)
        file<<layer.at(i)->scriptEntry();
      if(legend)
      {
        std::vector<std::string> legendText;
        for(UInt i=0; i<layer.size(); i++)
          legendText.push_back(layer.at(i)->legendEntry());

        legend->writeDataFile(plotBasics.workingDirectory, legendText);
        file<<legend->scriptEntry();
        file<<std::endl;
      }
      file<<plotBasics.scriptTrailer();
    } // end scriptfile

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
