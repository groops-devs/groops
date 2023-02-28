/***********************************************/
/**
* @file plotMisc.cpp
*
* @brief Miscellaneous classes for plots.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2010-01-13
*
*/
/***********************************************/

#define DOCSTRING_PlotColor
#define DOCSTRING_PlotLine
#define DOCSTRING_PlotSymbol

#include "base/import.h"
#include "base/string.h"
#include "config/configRegister.h"
#include "inputOutput/system.h"
#include "inputOutput/file.h"
#include "files/fileMatrix.h"
#include "files/fileStringTable.h"
#include "plotMisc.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(PlotColor, "plotColorType")
GROOPS_READCONFIG_CLASS(PlotColor, "plotColorType")

/***********************************************/

PlotColor::PlotColor(Config &config, const std::string &name)
{
  try
  {
    const std::vector<std::string> colorStrings = {"black", "red", "blue", "forestgreen", "darkorange", "darkred", "yellow", "green", "gray"};

    std::string choice;
    readConfigChoice(config, name, choice, Config::MUSTSET, "", "color");
    if(readConfigChoiceElement(config, "black",      choice)) colorStr = colorStrings.at(0);
    if(readConfigChoiceElement(config, "red",        choice)) colorStr = colorStrings.at(1);
    if(readConfigChoiceElement(config, "blue",       choice)) colorStr = colorStrings.at(2);
    if(readConfigChoiceElement(config, "green",      choice)) colorStr = colorStrings.at(3);
    if(readConfigChoiceElement(config, "orange",     choice)) colorStr = colorStrings.at(4);
    if(readConfigChoiceElement(config, "darkred",    choice)) colorStr = colorStrings.at(5);
    if(readConfigChoiceElement(config, "yellow",     choice)) colorStr = colorStrings.at(6);
    if(readConfigChoiceElement(config, "lightgreen", choice)) colorStr = colorStrings.at(7);
    if(readConfigChoiceElement(config, "gray",       choice)) colorStr = colorStrings.at(8);
    if(readConfigChoiceElement(config, "rgb", choice))
    {
      UInt r,g,b;
      readConfig(config, "red",   r, Config::MUSTSET, "", "0..255");
      readConfig(config, "green", g, Config::MUSTSET, "", "0..255");
      readConfig(config, "blue",  b, Config::MUSTSET, "", "0..255");
      std::stringstream ss;
      ss<<"#"<<std::setfill('0')<<std::hex<<std::setw(2)<<r<<std::setw(2)<<g<<std::setw(2)<<b;
      colorStr = ss.str();
    }
    if(readConfigChoiceElement(config, "grayscale", choice))
    {
      UInt value;
      readConfig(config, "value",   value , Config::MUSTSET, "", "0..255");
      std::stringstream ss;
      ss<<"#"<<std::setfill('0')<<std::hex<<std::setw(2)<<value<<std::setw(2)<<value<<std::setw(2)<<value;
      colorStr = ss.str();
    }
    if(readConfigChoiceElement(config, "namedColor", choice))
      readConfig(config, "colorName", colorStr, Config::MUSTSET, "", "name after GMT definition");
    if(readConfigChoiceElement(config, "cycler", choice))
    {
      UInt     index;
      FileName fileNameColor;
      readConfig(config, "index",              index,         Config::MUSTSET,  "0", "pick color based on index expression");
      readConfig(config, "inputfileColorList", fileNameColor, Config::OPTIONAL, "{groopsDataDir}/plot/colors.conf", "list of colors as defined by GMT");

      if(!isCreateSchema(config))
      {
        colorStr = colorStrings.at(1  + (index % (colorStrings.size()-1))); // exclude black
        if(!fileNameColor.empty())
        {
          std::vector<std::string> customColors;
          readFileStringList(fileNameColor, customColors);
          if(customColors.size() == 0)
            throw(Exception("File <" + fileNameColor.str() + "> does not contain any colors."));
          colorStr = customColors.at(index % customColors.size());
        }
      }
    }
    endChoice(config);
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotLineSolid = R"(
\subsection{Solid}
Draws a solid line.
)";

class PlotLineSolid : public PlotLine
{
public:
  PlotLineSolid(Config &config)
  {
    style = "solid";
    readConfig(config, "width", width, Config::DEFAULT, "1.5", "line width [p]");
    readConfig(config, "color", color, Config::MUSTSET, "",    "");
  }
};

/***********************************************/

static const char *docstringPlotLineDashed = R"(
\subsection{Dashed}
Draws a dashed line.
)";

class PlotLineDashed : public PlotLine
{
public:
  PlotLineDashed(Config &config)
  {
    style = "dashed";
    readConfig(config, "width", width, Config::DEFAULT, "1.5", "line width [p]");
    readConfig(config, "color", color, Config::MUSTSET, "",    "");
  }
};

/***********************************************/

static const char *docstringPlotLineDotted = R"(
\subsection{Dotted}
Draws a dotted line.
)";

class PlotLineDotted : public PlotLine
{
public:
  PlotLineDotted(Config &config)
  {
    style = "dotted";
    readConfig(config, "width", width, Config::DEFAULT, "1.5", "line width [p]");
    readConfig(config, "color", color, Config::MUSTSET, "", "  ");
  }
};

/***********************************************/

static const char *docstringPlotLineCustom = R"(
\subsection{Custom}
Draws a custom line. The line \config{style} code is described in
\url{https://docs.generic-mapping-tools.org/latest/cookbook/features.html#specifying-pen-attributes}.
)";

class PlotLineCustom : public PlotLine
{
public:
  PlotLineCustom(Config &config)
  {
    readConfig(config, "style", style, Config::MUSTSET, ".-",  "line style code");
    readConfig(config, "width", width, Config::DEFAULT, "1.5", "line width [p]");
    readConfig(config, "color", color, Config::MUSTSET, "",    "");
  }
};

/***********************************************/

GROOPS_REGISTER_CLASS(PlotLine, "plotLineType",
                      PlotLineSolid,
                      PlotLineDashed,
                      PlotLineDotted,
                      PlotLineCustom)

GROOPS_READCONFIG_CLASS(PlotLine, "plotLineType")

/***********************************************/

PlotLinePtr PlotLine::create(Config &config, const std::string &name)
{
  try
  {
    PlotLinePtr line;
    std::string  type;

    readConfigChoice(config, name, type, Config::MUSTSET, "", "line style");
    if(readConfigChoiceElement(config, "solid",   type, "solid line"))
      line = PlotLinePtr(new PlotLineSolid(config));
    if(readConfigChoiceElement(config, "dashed",  type, "dashed line"))
      line = PlotLinePtr(new PlotLineDashed(config));
    if(readConfigChoiceElement(config, "dotted",  type, "dotted line"))
      line = PlotLinePtr(new PlotLineDotted(config));
    if(readConfigChoiceElement(config, "custom",  type, "custom line (e.g. dash-dot)"))
      line = PlotLinePtr(new PlotLineCustom(config));
    endChoice(config);

    return line;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotLine::str() const
{
  try
  {
    std::stringstream ss;
    ss<<width<<"p,"<<color->str()<<","<<style;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(PlotSymbol, "plotSymbolType")
GROOPS_READCONFIG_CLASS(PlotSymbol, "plotSymbolType")

/***********************************************/

PlotSymbol::PlotSymbol(Config &config, const std::string &name)
{
  try
  {
    std::string  symbol;
    PlotColorPtr color;
    Double       size = 3;
    Bool         contour = 0;
    std::string  choice;

    const std::vector<std::pair<std::string, std::string>> symbols =
      {{"circle", "c"},{"star", "a"}, {"cross", "x"}, {"square", "s"}, {"triangle", "t"}, {"diamond",  "d"}, {"dash", "-"}};

    readConfigChoice(config, name, choice, Config::MUSTSET, "", "symbol");
    for(auto &s : symbols)
      if(readConfigChoiceElement(config, s.first, choice))
      {
        readConfig(config, "color",        color,   Config::OPTIONAL, "black", "empty: determined from value");
        readConfig(config, "size",         size,    Config::DEFAULT,  "3",     "size of symbol [point]");
        readConfig(config, "blackContour", contour, Config::DEFAULT,  "0",     "");
        symbol = s.second;
      }
    endChoice(config);
    if(isCreateSchema(config)) return;

    colorFromValue_ = (color == nullptr);

    {
      std::stringstream ss;
      ss<<symbol<<size<<"p ";
      if(contour)
        ss<<"-W"<<size/6.<<"p,black";
      if(color)
        ss<<" -G"<<color->str();
      else
        ss<<" -CgroopsPlot.cpt";
      scriptLine = ss.str();
    }

    {
      std::stringstream ss;
      ss<<symbol<<" "<<std::max(size, 3.)<<"p "<<(color ? color->str() : "white");
      if(contour || !color)
        ss<<" "<<size/6.<<"p,black";
      else
        ss<<" -";
      legendLine = ss.str();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void PlotBasics::read(Config &config, const std::string &nameProgram_, const FileName &fileNamePlot_, const std::string &title_,
                      const std::string &defaultWidth, const std::string &defaultHeight, const std::vector<std::string> &optionsDefault)
{
  try
  {
    std::vector<std::string> options;
    std::string optionsDefaultStr = R"(["FONT_ANNOT_PRIMARY=10p", "MAP_ANNOT_OFFSET_PRIMARY=0.1c", "FONT_LABEL=10p", "MAP_LABEL_OFFSET=0.1c", "FONT_ANNOT_SECONDARY=10p")";
    for(const std::string &str : optionsDefault)
      optionsDefaultStr += ", \""+str+"\"";
    optionsDefaultStr += "]";

    height = NAN_EXPR;

    if(readConfigSequence(config, "options", Config::MUSTSET,  "", "further options..."))
    {
      readConfig(config, "width",          width,          Config::DEFAULT,  defaultWidth,  "in cm");
      readConfig(config, "height",         height,         Config::OPTIONAL, defaultHeight, "in cm");
      readConfig(config, "titleFontSize",  titleSize,      Config::DEFAULT,  "12",  "in pt");
      readConfig(config, "marginTitle",    marginTitle,    Config::DEFAULT,  "0.4", "between title and figure [cm]");
      readConfig(config, "drawGridOnTop",  drawGridOnTop,  Config::DEFAULT,  "0",   "grid lines above all other lines/points");
      readConfig(config, "options",        options,        Config::OPTIONAL, optionsDefaultStr, "");
      readConfig(config, "transparent",    transparent,    Config::DEFAULT,  "1",   "make background transparent");
      readConfig(config, "dpi",            dpi,            Config::DEFAULT,  "300", "use this resolution when rasterizing postscript file");
      readConfig(config, "removeFiles",    removeFiles,    Config::DEFAULT,  "1",   "remove .gmt and script files");
      endSequence(config);
    }
    if(isCreateSchema(config))
      return;

    for(UInt i=0; i<options.size(); i++)
      if(!options.at(i).empty())
        optionsString += options.at(i)+" ";

    nameProgram  = nameProgram_;
    fileNamePlot = fileNamePlot_;
    title        = title_;

    baseDirectory = fileNamePlot.directory();
    if(std::getenv("GROOPS_PLOTDIR") && removeFiles)
      baseDirectory = FileName(std::getenv("GROOPS_PLOTDIR"));

    workingDirectory = baseDirectory.append("groopsPlot_"+fileNamePlot.baseName().str());
    System::createDirectories(workingDirectory);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileName PlotBasics::fileNameScript() const
{
#ifdef _WIN32
  return workingDirectory.append("groopsPlot.bat");
#else
  return workingDirectory.append("groopsPlot.sh");
#endif
}

/***********************************************/

std::string PlotBasics::scriptHeader() const
{
  std::stringstream ss;
#ifdef _WIN32
  ss<<"@echo off"<<std::endl;
  ss<<"REM automatically generated by groops::"<<nameProgram<<std::endl;
#else
  ss<<"#!/bin/bash"<<std::endl;
  ss<<"# automatically generated by groops::"<<nameProgram<<std::endl;
#endif
  ss<<scriptSetVariable("width",  width%"%gc"s)<<std::endl;
  ss<<scriptSetVariable("height", height%"%gc"s)<<std::endl;
  ss<<std::endl;
  ss<<"gmt set PS_MEDIA=a0 PS_PAGE_ORIENTATION=portrait"<<std::endl;
  ss<<"gmt set "<<optionsString<<std::endl;
  ss<<"gmt psxy -Y60c -X40c -JX1c/1c -R0/1/0/1 -K -T > groopsPlot.ps"<<std::endl; // starts upper left

  OutFile titleFile(workingDirectory.append("title.txt"));
  titleFile<<"0.5 0.0 "<<title;
  ss<<"gmt pstext title.txt -F+f"<<titleSize<<"p+jCB -Y-"<<titleSize<<"p -JX"<<width<<"c/"<<titleSize<<"p -R0/1/0/1 -N -O -K >> groopsPlot.ps "<<PlotBasics::scriptError2Null()<<std::endl;
  return ss.str();
}

/***********************************************/

std::string PlotBasics::scriptTrailer() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt psxy -R -J -O -T >> groopsPlot.ps"<<std::endl;

    // convert
    const std::string extension = String::lowerCase(fileNamePlot.fullExtension());
    if(extension != "ps")
    {
      const std::map<std::string, Char> extensionMapping = {{"bmp", 'b'}, {"eps", 'e'}, {"pdf", 'f'}, {"jpg", 'j'}, {"png", 'g'}, {"tif", 't'}, {"ppm", 'm'}};
      if(extensionMapping.find(extension) == extensionMapping.end())
        throw(Exception("Unsupported image format <"+extension+">."));
      Char outputFormat = extensionMapping.at(extension);
      if(outputFormat == 'g' && transparent)
        outputFormat = 'G';
      ss<<"gmt psconvert groopsPlot.ps -Z -A1p -E"<<dpi<<" -Qg1 -Qt4 -T"<<outputFormat<<std::endl;
    }

    // copy plot to destination
#ifdef _WIN32
    std::string cp = "copy";
#else
    std::string cp = "cp";
#endif
    if(baseDirectory.str() == fileNamePlot.directory().str())
      ss<<cp<<" groopsPlot."<<extension<<" "<<FileName("..").append(fileNamePlot.stripDirectory())<<std::endl;
    else
    {
      ss<<"cd \""<<System::currentWorkingDirectory()<<"\" && ";
      ss<<cp<<" \""<<fileNameScript().directory().append("groopsPlot."+extension)<<"\" \""<<fileNamePlot<<"\""<<std::endl;
    }

    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool PlotBasics::runScript() const
{
  try
  {
    std::stringstream ss;
    if(!fileNameScript().directory().empty())
      ss<<"cd \""<<fileNameScript().directory()<<"\" && ";
#ifndef _WIN32
    ss<<"bash ";
#endif
    ss<<fileNameScript().stripDirectory();
    System::exec(ss.str());

    if(removeFiles)
      System::remove(fileNameScript().directory());

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotBasics::axisTicks(Bool isLog, Double /*vmin*/, Double /*vmax*/, Double annotation, Double frame, Double grid, const std::string &unit, const std::string &label)
{
  try
  {
    if(isLog)
    {
      if(std::isnan(annotation)) annotation = 1;
      if(std::isnan(frame))      frame      = 3;
      if(std::isnan(grid))       grid       = 3;
    }

    std::stringstream ss;
    if(annotation != 0.) {ss<<"a"; if(!std::isnan(annotation)) ss<<annotation; if(isLog) ss<<"p";} // major ticks
    if(frame      != 0.) {ss<<"f"; if(!std::isnan(frame))      ss<<frame;}      // minor ticks
    if(grid       != 0.) {ss<<"g"; if(!std::isnan(grid))       ss<<grid;}       // grid line spacing
    if(!unit.empty())  ss<<"+u\""<<(unit=="%" ? "" : " ")<<unit<<"\"";
    if(!label.empty()) ss<<"+l\""<<label<<"\"";
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotBasics::scriptSetVariable(const std::string &var, const std::string &value)
{
#ifdef _WIN32
  return "set "+var+"="+value;
#else
  return var+"="+value;
#endif
}

/***********************************************/

std::string PlotBasics::scriptVariable(const std::string &var)
{
#ifdef _WIN32
  return "%"+var+"%";
#else
  return "$"+var;
#endif
}

/***********************************************/

std::string PlotBasics::scriptError2Null()
{
#ifdef _WIN32
  return "2>nul";
#else
  return "2>/dev/null";
#endif
}

/***********************************************/

UInt PlotBasics::gmtVersion()
{
  try
  {
    std::vector<std::string> output;
    if(!System::exec("gmt --version", output) || (output.size() == 0))
      return MAX_UINT; // GMT does not seem to be installed -> use always latest version

    const UInt version = std::stoul(output.at(0).substr(0,1))*100 + std::stoul(output.at(0).substr(2,1))*10 + std::stoul(output.at(0).substr(4,1));
    return version;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
