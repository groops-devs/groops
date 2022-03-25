/***********************************************/
/**
* @file plotLegend.cpp
*
* @brief Legend.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2020-05-03
*
*/
/***********************************************/

#define DOCSTRING_PlotLegend

#include "base/import.h"
#include "base/string.h"
#include "config/configRegister.h"
#include "inputOutput/file.h"
#include "plot/plotMisc.h"
#include "plotLegend.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(PlotLegend, "plotLegendType")
GROOPS_READCONFIG_CLASS(PlotLegend, "plotLegendType")

/***********************************************/

PlotLegend::PlotLegend(Config &config, const std::string &name)
{
  try
  {
    height = 0;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "width",       width,     Config::DEFAULT,  "10",   "legend width [cm]");
    readConfig(config, "height",      height,    Config::OPTIONAL, "",     "legend height [cm] (default: estimated)");
    readConfig(config, "positionX",   positionX, Config::DEFAULT,  "1.05", "legend x-position in normalized (0-1) coordinates.");
    readConfig(config, "positionY",   positionY, Config::DEFAULT,  "1.0",  "legend y-position in normalized (0-1) coordinates.");
    readConfig(config, "anchorPoint", anchor,    Config::DEFAULT,  "TL",   "Two character combination of L, C, R (for left, center, or right) and T, M, B for top, middle, or bottom. e.g., TL for top left");
    readConfig(config, "columns",     ncolumns,  Config::DEFAULT,  "1",    "number of columns in legend");
    readConfig(config, "textColor",   textColor, Config::OPTIONAL, "",     "color of the legend text");
    readConfig(config, "fillColor",   fillColor, Config::OPTIONAL, "",     "fill color of the legend box");
    readConfig(config, "edgeLine",    edgeLine,  Config::OPTIONAL, "",     "style of the legend box edge");
    endSequence(config);
    if(isCreateSchema(config)) return;

    _hasEntries = FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotLegend::writeDataFile(const FileName &workingDirectory, const std::vector<std::string> &legendText)
{
  try
  {
    if(legendText.size()==0)
      return;

    _hasEntries = TRUE;
    OutFile legendFile(workingDirectory.append("legend.txt"));

    if(textColor)
      legendFile<<"C "<<textColor->str()<<std::endl;
    legendFile<<"N "<<ncolumns<<std::endl;
    for(UInt i=0; i<legendText.size(); i++)
    {
      std::vector<std::string> lines = String::split(legendText.at(i), '\n');
      for(const std::string &line : lines)
        if(!line.empty())
          legendFile<<line<<std::endl;
    }
  }
  catch(const std::exception& e)
  {
    GROOPS_RETHROW(e)
  }

}

/***********************************************/

std::string PlotLegend::scriptEntry() const
{
  try
  {
    if(!_hasEntries)
      return "";

    std::stringstream ss;
    ss<<"gmt pslegend legend.txt -R0/1/0/1 -JX"<<PlotBasics::scriptVariable("width")<<"/"<<PlotBasics::scriptVariable("height");
    ss<<" -Dn"<<positionX<<"/"<<positionY<<"+w"<<width<<"c/"<<height<<"c+j"<<anchor<<"+l1.5";
    if(fillColor || edgeLine)
    {
      ss<<" -F";
      if(fillColor)
        ss<<"+g"<<fillColor->str();
      if(edgeLine)
        ss<<"+p"<<edgeLine->str();
      else
        ss<<"+p0.0p,"<<fillColor->str();
    }
    ss<<" -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
