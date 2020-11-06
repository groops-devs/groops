/***********************************************/
/**
* @file plotGraph.cpp
*
* @brief Plot functions (x,y).
*
* @author Enrico Kurtenbach
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2015-10-20
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Generates a two dimensional xy plot using the GMT Generic Mapping Tools (\url{https://www.generic-mapping-tools.org}).
A variety of image file formats are supported (e.g. png, jpg, eps) determined by the extension of \config{outputfile}.

The plotting area is defined by the two axes \configClass{axisX/Y}{plotAxisType}. An alternative \configClass{axisY2}{plotAxisType}
on the right hand side can be added. The content of the graph itself is defined
by one or more \configClass{layer}{plotGraphLayerType}s.

The plot programs create a temporary directory in the path of \config{outputfile}, writes all needed data into it,
generates a batch/shell script with the GMT commands, execute it, and remove the temporary directory.
With setting \config{options:removeFiles}=false the last step is skipped and it is possible to adjust the plot manually
to specific publication needs. Individual GMT settings are adjusted with \config{options:options}="\verb|FORMAT=value|",
see \url{https://docs.generic-mapping-tools.org/latest/gmt.conf.html}.

See also: \program{PlotDegreeAmplitudes}, \program{PlotMap}, \program{PlotMatrix}, \program{PlotSphericalHarmonicsTriangle}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "plot/plotMisc.h"
#include "plot/plotAxis.h"
#include "plot/plotGraphLayer.h"
#include "plot/plotLegend.h"
#include "plot/plotColorbar.h"

/***********************************************/

/** @brief Plot functions (x,y).
* @ingroup programsGroup */
class PlotGraph
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(PlotGraph, SINGLEPROCESS, "Plot functions (x,y)", Plot, TimeSeries)

/***********************************************/

void PlotGraph::run(Config &config)
{
  try
  {
    FileName          fileNamePlot;
    std::string       title;
    std::vector<PlotGraphLayerPtr> layer;
    PlotAxisPtr       axisX, axisY, axisY2;
    PlotColorbarPtr   colorbar;
    PlotLegendPtr     legend;
    PlotBasics        plotBasics;

    readConfig(config, "outputfile", fileNamePlot, Config::MUSTSET,  "",  "*.png, *.jpg, *.eps, ...");
    readConfig(config, "title",      title,        Config::OPTIONAL, "",  "");
    readConfig(config, "layer",      layer,        Config::MUSTSET,  "",  "");
    readConfig(config, "axisX",      axisX,        Config::MUSTSET,  "",  "");
    readConfig(config, "axisY",      axisY,        Config::MUSTSET,  "",  "");
    readConfig(config, "axisY2",     axisY2,       Config::OPTIONAL, "",  "Second y-axis on right hand side");
    readConfig(config, "colorbar",   colorbar,     Config::OPTIONAL, "",  "");
    readConfig(config, "legend",     legend,       Config::OPTIONAL, "1", "");
    plotBasics.read(config, "PlotGraph", fileNamePlot, title, "12", "7", {"TIME_SYSTEM=MJD"s});
    if(isCreateSchema(config)) return;

    // ===========================================================

    // colorbar needed?
    // ----------------
    for(UInt i=0; i<layer.size(); i++)
      if(layer.at(i)->requiresColorBar() && !colorbar)
        throw(Exception("layer need definition of a colorbar"));

    // calculate min and max
    // ---------------------
    // x-axis
    Double minX  =  1e99;
    Double maxX  = -1e99;
    for(UInt k=0; k<layer.size(); k++)
      layer.at(k)->getIntervalX(axisX->isLogarithmic(), minX, maxX);
    axisX->setAutoInterval(minX, maxX);
    logInfo<<"  plot x range  ("<<axisX->getMin()<<" .. "<<axisX->getMax()<<") of ("<<minX<<" .. "<<maxX<<")"<<Log::endl;

    // y-axis
    Double minY  =  1e99;
    Double maxY  = -1e99;
    Double minY2 =  1e99;
    Double maxY2 = -1e99;
    for(UInt k=0; k<layer.size(); k++)
    {
      if(axisY2 && layer.at(k)->drawOnSecondAxis())
        layer.at(k)->getIntervalY(axisY2->isLogarithmic(), axisX->getMin(), axisX->getMax(), minY2, maxY2);
      else
        layer.at(k)->getIntervalY(axisY->isLogarithmic(), axisX->getMin(), axisX->getMax(), minY, maxY);
    }
    axisY->setAutoInterval(minY, maxY);
    if(axisY2)
      axisY2->setAutoInterval(minY2, maxY2);
    logInfo<<"  plot y range  ("<<axisY->getMin()<<" .. "<<axisY->getMax()<<") of ("<<minY<<" .. "<<maxY<<")"<<Log::endl;
    if(axisY2)
      logInfo<<"  plot y2 range ("<<axisY2->getMin()<<" .. "<<axisY2->getMax()<<") of ("<<minY2<<" .. "<<maxY2<<")"<<Log::endl;

    // z-axis
    if(colorbar)
    {
      Double minZ  =  1e99;
      Double maxZ  = -1e99;
      for(UInt k = 0; k<layer.size(); k++)
        if(axisY2 && layer.at(k)->drawOnSecondAxis())
          layer.at(k)->getIntervalZ(colorbar->isLogarithmic(), minX, maxX, minY2, maxY2, minZ, maxZ);
        else
          layer.at(k)->getIntervalZ(colorbar->isLogarithmic(), minX, maxX, minY, maxY, minZ, maxZ);
      colorbar->setAutoInterval(minZ, maxZ);
      logInfo<<"  plot z range  ("<<colorbar->getMin()<<" .. "<<colorbar->getMax()<<") of ("<<minZ<<" .. "<<maxZ<<")"<<Log::endl;
    }

    // create data files
    // -----------------
    logStatus<<"create temporary data files"<<Log::endl;
    for(UInt i=0; i<layer.size(); i++)
      if(axisY2 && layer.at(i)->drawOnSecondAxis())
        layer.at(i)->writeDataFile(plotBasics.workingDirectory, i, minX, maxX, minY2, maxY2);
      else
        layer.at(i)->writeDataFile(plotBasics.workingDirectory, i, minX, maxX, minY, maxY);
    axisX->writeDataFile(plotBasics.workingDirectory, "x");
    axisY->writeDataFile(plotBasics.workingDirectory, "y");
    if(axisY2) axisY2->writeDataFile(plotBasics.workingDirectory, "y2");

    // intervals
    // ---------
    std::string JX, JX2;
    {
      std::stringstream ss;
      ss<<"-R"<<PlotBasics::scriptVariable("range")<<" -JX"<<axisX->direction()*plotBasics.width<<"c"<<axisX->axisModifier()<<"/"<<axisY->direction()*plotBasics.height<<"c"<<axisY->axisModifier();
      JX = ss.str();
    }
    if(axisY2)
    {
      std::stringstream ss;
      ss<<"-R"<<PlotBasics::scriptVariable("rangeY2")<<" -JX"<<axisX->direction()*plotBasics.width<<"c"<<axisX->axisModifier()<<"/"<<axisY2->direction()*plotBasics.height<<"c"<<axisY2->axisModifier();
      JX2 = ss.str();
    }

    // create scriptfile
    // -----------------
    logStatus<<"create scriptfile"<<Log::endl;
    {
      OutFile file(plotBasics.fileNameScript());
      file<<plotBasics.scriptHeader();
      if(colorbar)
        file<<colorbar->scriptColorTable();
      file<<PlotBasics::scriptSetVariable("range", axisX->getMin()%"%.12e/"s+axisX->getMax()%"%.12e/"s+axisY->getMin()%"%.12e/"s+axisY->getMax()%"%.12e"s)<<std::endl;
      if(axisY2)
        file<<PlotBasics::scriptSetVariable("rangeY2", axisX->getMin()%"%.12e/"s+axisX->getMax()%"%.12e/"s+axisY2->getMin()%"%.12e/"s+axisY2->getMax()%"%.12e"s)<<std::endl;
      // x-axis
      file<<"gmt psbasemap -Y"<<-fabs(plotBasics.height)-plotBasics.marginTitle<<"c "<<JX<<" -BSn "<<axisX->scriptEntry("x")<<" -O -K >> groopsPlot.ps"<<std::endl;
      // y-axis
      file<<"gmt psbasemap "<<JX<<" -BW"<<((axisY2) ? " " : "e ")<<axisY->scriptEntry("y")<<" -O -K >> groopsPlot.ps"<<std::endl;
      Bool isAxisY2 = FALSE;
      // second y-axis
      if(axisY2)
      {
        file<<"gmt psbasemap "<<JX2<<" -BE "<<axisY2->scriptEntry("y")<<" -O -K >> groopsPlot.ps"<<std::endl;
        isAxisY2 = TRUE;
      }

      // data sets
      for(UInt i=0; i<layer.size(); i++)
      {
        if(isAxisY2 && !layer.at(i)->drawOnSecondAxis())
        {
          file<<"gmt psxy "<<JX<<" -Bw -T -O -K >> groopsPlot.ps"<<std::endl;
          isAxisY2 = FALSE;
        }
        else if((!isAxisY2) && axisY2 && layer.at(i)->drawOnSecondAxis())
        {
          file<<"gmt psxy "<<JX2<<" -Be -T -O -K >> groopsPlot.ps"<<std::endl;
          isAxisY2 = TRUE;
        }
        file<<layer.at(i)->scriptEntry();
      }
      file<<std::endl;

      // write axis again
      file<<"gmt psbasemap "<<JX<<" -BSn "<<axisX->scriptEntry("x", plotBasics.drawGridOnTop)<<" -O -K >> groopsPlot.ps"<<std::endl;
      file<<"gmt psbasemap "<<JX<<" -BW"<<((axisY2) ? " " : "e ")<<axisY->scriptEntry("y", plotBasics.drawGridOnTop)<<" -O -K >> groopsPlot.ps"<<std::endl;
      if(axisY2)
        file<<"gmt psbasemap "<<JX2<<" -BE "<<axisY2->scriptEntry("y", plotBasics.drawGridOnTop)<<" -O -K >> groopsPlot.ps"<<std::endl;

      if(colorbar)
        file<<colorbar->scriptEntry(plotBasics.width, plotBasics.height, ((axisY2) ? axisY2->getMargin() : 0.0), axisX->getMargin());

      if(legend)
      {
        std::vector<std::string> legendText;
        for(UInt i=0; i<layer.size(); i++)
          legendText.push_back(layer.at(i)->legendEntry());

        legend->writeDataFile(plotBasics.workingDirectory, legendText);
        file<<legend->scriptEntry();
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
