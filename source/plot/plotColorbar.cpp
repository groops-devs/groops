/***********************************************/
/**
* @file plotColorbar.cpp
*
* @brief Colorbar.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2020-05-03
*
*/
/***********************************************/

#define DOCSTRING_PlotColorbar

#include "base/import.h"
#include "config/configRegister.h"
#include "plotMisc.h"
#include "plotColorbar.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(PlotColorbar, "plotColorbarType")
GROOPS_READCONFIG_CLASS(PlotColorbar, "plotColorbarType")

/***********************************************/

PlotColorbar::PlotColorbar(Config &config, const std::string &name)
{
  try
  {
    vMin       = NAN_EXPR;
    vMax       = NAN_EXPR;
    annotation = NAN_EXPR;
    Bool reverse;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "min",           vMin,          Config::OPTIONAL, "",    "");
    readConfig(config, "max",           vMax,          Config::OPTIONAL, "",    "");
    readConfig(config, "annotation",    annotation,    Config::OPTIONAL, "",    "boundary annotation");
    readConfig(config, "unit",          unit,          Config::OPTIONAL, "",    "appended to axis values");
    readConfig(config, "label",         label,         Config::OPTIONAL, "",    "description of the axis");
    readConfig(config, "logarithmic",   isLog,         Config::DEFAULT,  "0",   "use logarithmic scale");
    readConfig(config, "triangleLeft",  triangleLeft,  Config::DEFAULT,  "1",   "");
    readConfig(config, "triangleRight", triangleRight, Config::DEFAULT,  "1",   "");
    readConfig(config, "illuminate",    illuminate,    Config::DEFAULT,  "0",   "illuminate");
    readConfig(config, "vertical",      vertical,      Config::DEFAULT,  "0",   "plot vertical color bar on the right");
    readConfig(config, "length",        length,        Config::DEFAULT,  "100", "length of colorbar in percent");
    readConfig(config, "margin",        margin,        Config::DEFAULT,  "0.4", "between colorbar and figure [cm]");
    readConfig(config, "colorTable",    colorTable,    Config::DEFAULT,  "haxby", "name of the color bar");
    readConfig(config, "reverse",       reverse,       Config::DEFAULT,  "0",   "reverse direction");
    readConfig(config, "showColorbar",  isPlot,        Config::DEFAULT,  "1",   "");
    endSequence(config);
    if(isCreateSchema(config)) return;

    if(reverse)
      colorTable += " -I";
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotColorbar::setAutoInterval(Double minAuto, Double maxAuto)
{
  try
  {
    if(std::isnan(vMin)) vMin = minAuto;
    if(std::isnan(vMax)) vMax = maxAuto;

    if(vMin == vMax)
    {
      vMin *= 0.9;
      vMax *= 1.1;
      if(vMin == 0.0)
      {
        vMin = -1.0;
        vMax = +1.0;
      }
    }

    if(vMin > vMax) std::swap(vMin, vMax);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotColorbar::scriptColorTable() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt makecpt -D -Z -C"<<((colorTable.empty()) ? "haxby" : colorTable)<<" -T"<<vMin<<"/"<<vMax<<"/";
    if(isLog)
      ss<<"3 -Qo";
    else
      ss<<((vMax-vMin)/9.)%"%.10g"s;
    ss<<" >groopsPlot.cpt "<<PlotBasics::scriptError2Null()<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotColorbar::scriptEntry(Double width, Double height, Double marginX, Double marginY) const
{
  try
  {
    if(!isPlot)
      return "";
    std::stringstream ss;
    ss<<"gmt psscale -CgroopsPlot.cpt";
    if(!isLog)
    {
      ss<<" -Ba"; if(!std::isnan(annotation)) ss<<annotation;
    }
    else
      ss<<" -Q -Ba1pf3";
    if(!label.empty()) ss<<"+l\""<<label<<"\"";
    if(!unit.empty())  ss<<"+u\""<<(unit=="%" ? "" : " ")<<unit<<"\"";
    if(illuminate)     ss<<" -I";
    if(PlotBasics::gmtVersion() < 531)
    {
      if(triangleLeft)   ss<<" -Eb";
      if(triangleRight)  ss<<" -Ef";
      if(vertical)
        ss<<" -D"<<width+marginX+margin<<"c/"<<height/2<<"c/"<<(height-0.4)*length/100.<<"c/0.4c";
      else
        ss<<" -D"<<width/2<<"c/"<<-margin-marginY<<"c/"<<(width-0.4)*length/100.<<"c/0.4ch";
    }
    else
    {
      if(vertical)
        ss<<" -Dx"<<width+marginX+margin<<"c/"<<height/2<<"c+w"<<(height-0.4)*length/100.<<"c/0.4c+jLM";
      else
        ss<<" -Dx"<<width/2<<"c/"<<-margin-marginY<<"c+w"<<(width-0.4)*length/100.<<"c/0.4c+jCT+h";
      if(triangleLeft || triangleRight) ss<<"+e";
      if(triangleLeft)   ss<<"b";
      if(triangleRight)  ss<<"f";
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
