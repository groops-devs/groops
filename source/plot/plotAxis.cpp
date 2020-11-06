/***********************************************/
/**
* @file plotAxis.cpp
*
* @brief Axis ticks, limits and labels.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2016-07-23
*
*/
/***********************************************/

#define DOCSTRING_PlotAxis

#include "base/import.h"
#include "config/configRegister.h"
#include "plot/plotMisc.h"
#include "plotAxis.h"
#include "inputOutput/file.h"
#include "inputOutput/system.h"

/***********************************************/

static const char *docstringPlotAxisStandard = R"(
\subsection{Standard}
General axis for arbitrary input data.
)";

/***** CLASS ***********************************/

class PlotAxisStandard : public PlotAxis
{
  Double      annotation, frame, grid;
  std::string unit, label;
  Bool        isLog;

public:
  PlotAxisStandard(Config &config);
  void        setAutoInterval(Double minAuto, Double maxAuto);
  std::string axisModifier() const {return isLog ? "l" : "";}
  std::string scriptEntry(const std::string &axis, Bool withGrid) const;
  Bool        isLogarithmic() const {return isLog;}
};

/***********************************************/

PlotAxisStandard::PlotAxisStandard(Config &config)
{
  try
  {
    vMin       = NAN_EXPR;
    vMax       = NAN_EXPR;
    annotation = NAN_EXPR;
    frame      = NAN_EXPR;
    grid       = NAN_EXPR;

    renameDeprecatedConfig(config, "annotation", "majorTickSpacing", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "frame",      "minorTickSpacing", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "grid",       "gridLineSpacing",  date2time(2020, 4, 23));

    readConfig(config, "min",              vMin,            Config::OPTIONAL, "",  "The minimum value of the axis. If no value is given, the minimum scale value is set automatically.");
    readConfig(config, "max",              vMax,            Config::OPTIONAL, "",  "The maximum value of the axis. If no value is given, the maximum scale value is set automatically.");
    readConfig(config, "majorTickSpacing", annotation,      Config::OPTIONAL, "",  "The boundary annotation.");
    readConfig(config, "minorTickSpacing", frame,           Config::OPTIONAL, "",  "The spacing of the frame tick intervals.");
    readConfig(config, "gridLineSpacing",  grid,            Config::OPTIONAL, "",  "The spacing of the grid line intervals");
    readConfig(config, "gridLine",         gridLine,        Config::OPTIONAL, R"({"solid": {"width":"0.25", "color":"gray"}})", "The style of the grid lines.");
    readConfig(config, "unit",             unit,            Config::OPTIONAL, "",  "Naming unit to append to the axis values.");
    readConfig(config, "label",            label,           Config::OPTIONAL, "",  "The description of the axis.");
    readConfig(config, "logarithmic",      isLog,           Config::DEFAULT,  "0", "If set to 'yes', a logarithmic scale is used for the axis.");
    readConfig(config, "color",            color,           Config::MUSTSET,  "",  "Setting the color of the axis bars and labels.");
    readConfig(config, "changeDirection",  changeDirection, Config::DEFAULT,  "0", "If set to 'yes', the directions right/up are changed to left/down.");
    if(isCreateSchema(config)) return;

    margin = 0.5;
    if(!label.empty())
      margin += 0.3;

    if(!gridLine)
      grid = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotAxisStandard::setAutoInterval(Double minAuto, Double maxAuto)
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

std::string PlotAxisStandard::scriptEntry(const std::string &axis, Bool withGrid) const
{
  try
  {
    std::stringstream ss;
    ss<<" --MAP_DEFAULT_PEN=+"<<color->str()<<" --FONT_ANNOT_PRIMARY="<<color->str()<<" --FONT_LABEL="<<color->str();
    if(gridLine) ss<<" --MAP_GRID_PEN_PRIMARY="<<gridLine->str();
    ss<<" -B"<<axis<<PlotBasics::axisTicks(isLog, vMin, vMax, annotation, frame, (withGrid) ? grid : 0, unit, label);
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotAxisTime = R"(
\subsection{Time}
The input data are interpreted as MJD (modified Julian date).
The unit of the tick spacings should be appenend to the number and can be any of
\begin{itemize}
\item Y (year, plot with 4 digits)
\item y (year, plot with 2 digits)
\item O (month, plot using \verb|FORMAT_DATE_MAP|)
\item o (month, plot with 2 digits)
\item U (ISO week, plot using \verb|FORMAT_DATE_MAP|)
\item u (ISO week, plot using 2 digits)
\item r (Gregorian week, 7-day stride from start of week \verb|TIME_WEEK_START|)
\item K (ISO weekday, plot name of day)
\item D (date, plot using \verb|FORMAT_DATE_MAP|)
\item d (day, plot day of month 0-31 or year 1-366, via \verb|FORMAT_DATE_MAP|)
\item R (day, same as d, aligned with \verb|TIME_WEEK_START|)
\item H (hour, plot using \verb|FORMAT_CLOCK_MAP|)
\item h (hour, plot with 2 digits)
\item M (minute, plot using \verb|FORMAT_CLOCK_MAP|)
\item m (minute, plot with 2 digits)
\item S (second, plot using \verb|FORMAT_CLOCK_MAP|)
\item s (second, plot with 2 digits).
\end{itemize}

A secondary time axis can be added to specify larger intervals (e.g dates of hourly data).

Examples: Settings for Fig.~\ref{plotAxisType:plotAxisTime1}: \config{majorTickSpacing}=\verb|6H|, secondary: \config{majorTickSpacing}=\verb|1D|.
\fig{!hb}{1.0}{plotAxisTime1}{plotAxisType:plotAxisTime1}{Time axis for daily data.}

Settings for Fig.~\ref{plotAxisType:plotAxisTime2}: \config{majorTickSpacing}=\verb|2d|, secondary: \config{majorTickSpacing}=\verb|1O|, \config{options}=\verb|FORMAT_DATE_MAP="o yyyy"|.
\fig{!hb}{1.0}{plotAxisTime2}{plotAxisType:plotAxisTime2}{Time axis for monthly data.}

Settings for Fig.~\ref{plotAxisType:plotAxisTime3}: \config{majorTickSpacing}=\verb|1o|, secondary: \config{majorTickSpacing}=\verb|1Y|, \config{options}=\verb|FORMAT_DATE_MAP="mm"|.
\fig{!hb}{1.0}{plotAxisTime3}{plotAxisType:plotAxisTime3}{Time axis for yearly data.}

)";

/***** CLASS ***********************************/

class PlotAxisTime : public PlotAxis
{
  std::string annotation,  grid,  frame;
  std::string annotation2, grid2, frame2;
  std::string optionsString;

public:
  PlotAxisTime(Config &config);
  void        setAutoInterval(Double minAuto, Double maxAuto);
  std::string axisModifier() const {return "t";}
  std::string scriptEntry(const std::string &axis, Bool withGrid) const;
};

/***********************************************/

PlotAxisTime::PlotAxisTime(Config &config)
{
  try
  {
    Time minTime, maxTime;
    std::vector<std::string> options;

    renameDeprecatedConfig(config, "annotation", "majorTickSpacing", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "frame",      "minorTickSpacing", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "grid",       "gridLineSpacing",  date2time(2020, 4, 23));

    readConfig(config, "min",                minTime,       Config::OPTIONAL, "",   "The minimum value of the time axis. If no value is given, the minimum scale value is set automatically.");
    readConfig(config, "max",                maxTime,       Config::OPTIONAL, "",   "The maximum value of the time axis. If no value is given, the maximum scale value is set automatically.");
    readConfig(config, "majorTickSpacing",   annotation,    Config::OPTIONAL, "2o", "Y: year, o: month");
    readConfig(config, "minorTickSpacing",   frame,         Config::OPTIONAL, "1o", "D: date, d: day");
    readConfig(config, "gridLineSpacing",    grid,          Config::OPTIONAL, "",   "H: clock, h: hour, m: minute, s: second");
    if(readConfigSequence(config, "secondary", Config::OPTIONAL, "", "secondary time axis"))
    {
      renameDeprecatedConfig(config, "annotation", "majorTickSpacing", date2time(2020, 4, 23));
      renameDeprecatedConfig(config, "frame",      "minorTickSpacing", date2time(2020, 4, 23));
      renameDeprecatedConfig(config, "grid",       "gridLineSpacing",  date2time(2020, 4, 23));

      readConfig(config, "majorTickSpacing", annotation2, Config::OPTIONAL, "1Y", "Y: year, o: month");
      readConfig(config, "minorTickSpacing", frame2,      Config::OPTIONAL, "",   "D: date, d: day");
      readConfig(config, "gridLineSpacing",  grid2,       Config::OPTIONAL, "1Y", "H: clock, h: hour, m: minute, s: second");
      endSequence(config);
    }
    readConfig(config, "color",           color,           Config::MUSTSET,  "",  "color of axis bars and labels");
    readConfig(config, "gridLine",        gridLine,        Config::OPTIONAL, R"({"solid": {"width":"0.25", "color":"gray"}})", "The style of the grid lines.");
    readConfig(config, "changeDirection", changeDirection, Config::DEFAULT,  "0", "right->left / up->down");
    readConfig(config, "options",         options,         Config::OPTIONAL, R"(["FORMAT_DATE_MAP=yyyy-mm-dd", "FORMAT_CLOCK_MAP=hh:mm"])", "adjust date format");
    if(isCreateSchema(config)) return;

    vMin = NAN_EXPR;
    vMax = NAN_EXPR;
    if(maxTime != Time()) vMax = maxTime.mjd();
    if(minTime != Time()) vMin = minTime.mjd();

    for(UInt i=0; i<options.size(); i++)
      if(!options.at(i).empty())
        optionsString += " --"+options.at(i);

    margin = 0.5;
    if(!(annotation2.empty() && frame2.empty() && grid2.empty()))
      margin +=0.3;

    if(!gridLine)
    {
      grid  = "";
      grid2 = "";
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotAxisTime::setAutoInterval(Double minAuto, Double maxAuto)
{
  try
  {
    if(std::isnan(vMin)) vMin = minAuto;
    if(std::isnan(vMax)) vMax = maxAuto;
    if(vMin > vMax) std::swap(vMin, vMax);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotAxisTime::scriptEntry(const std::string &axis, Bool withGrid) const
{
  try
  {
    std::string intervals;
    {
      std::stringstream ss;
      if(!annotation.empty())       ss<<"a"<<annotation;
      if(!frame.empty())            ss<<"f"<<frame;
      if(!grid.empty() && withGrid) ss<<"g"<<grid;
      intervals = ss.str();
    }
    std::string intervals2;
    {
      std::stringstream ss;
      if(!annotation2.empty())       ss<<"a"<<annotation2;
      if(!frame2.empty())            ss<<"f"<<frame2;
      if(!grid2.empty() && withGrid) ss<<"g"<<grid2;
      intervals2 = ss.str();
    }

    std::stringstream ss;
    ss<<" --MAP_DEFAULT_PEN=+"<<color->str()<<" --FONT_ANNOT_PRIMARY="<<color->str()<<" --FONT_LABEL="<<color->str();
    if(gridLine) ss<<" --MAP_GRID_PEN_PRIMARY="<<gridLine->str();
    ss<<optionsString;
    if(intervals2.empty())
      ss<<" -B"<<axis<<intervals;
    else
      ss<<" -Bp"<<axis<<intervals<<" -Bs"<<axis<<intervals2;

    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotAxisLabeled = R"(
\subsection{Labeled}
Axis with string labels. The coordinate system is based on the label indices (e.g. 0, 1, 2).
)";

/***** CLASS ***********************************/

class PlotAxisLabeled : public PlotAxis
{
  std::vector<std::string> labels;
  Bool     orthoLabels;
  FileName dataFileName;

public:
  PlotAxisLabeled(Config &config);
  void        setAutoInterval(Double /*minAuto*/, Double /*maxAuto*/) {}
  void        writeDataFile(const FileName &workingDirectory, const std::string &axis);
  std::string axisModifier() const {return "";}
  std::string scriptEntry(const std::string &axis, Bool withGrid) const;
};

/***********************************************/

PlotAxisLabeled::PlotAxisLabeled(Config &config)
{
  try
  {
    ExpressionVariablePtr exprMin, exprMax;

    readConfig(config, "labels",           labels,          Config::MUSTSET,  "",             "tick labels (ticks are placed at their index. e.g. 0, 1, ..., 5)");
    readConfig(config, "min",              exprMin,         Config::DEFAULT,  "0",            "minimum value of the axis");
    readConfig(config, "max",              exprMax,         Config::DEFAULT,  "labelCount-1", "maximum values of the axis");
    readConfig(config, "orthogonalLabels", orthoLabels,     Config::DEFAULT,  "0",            "labels are oriented orthogonal to axis");
    readConfig(config, "gridLine",         gridLine,        Config::OPTIONAL, R"({"solid": {"width":"0.25", "color":"gray"}})", "The style of the grid lines.");
    readConfig(config, "color",            color,           Config::MUSTSET,  "",             "set the color of the axis and labels");
    readConfig(config, "changeDirection",  changeDirection, Config::DEFAULT,  "0",            "If set to 'yes', the directions right/up are changed to left/down.");
    if(isCreateSchema(config)) return;

    margin = 0.5;

    VariableList fileNameVariableList;
    addVariable("labelCount", fileNameVariableList);
    fileNameVariableList["labelCount"]->setValue(static_cast<Double>(labels.size()));

    vMin = exprMin->evaluate(fileNameVariableList) - 0.5;
    vMax = exprMax->evaluate(fileNameVariableList) + 0.5;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotAxisLabeled::writeDataFile(const FileName &workingDirectory, const std::string &axis)
{
  try
  {
    dataFileName = "labels."+axis+".txt";
    OutFile file(workingDirectory.append(dataFileName));
    for(UInt k=0; k<labels.size(); k++)
      if((k >= vMin) && (k <= vMax))
      {
        file<<k<<" a "<<labels.at(k)<<std::endl;
        if(gridLine) file<<k<<" g"<<std::endl;
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/***********************************************/

std::string PlotAxisLabeled::scriptEntry(const std::string &axis, Bool /*withGrid*/) const
{
  try
  {
    std::stringstream ss;
    ss<<" --MAP_DEFAULT_PEN=+"<<color->str()<<" --FONT_ANNOT_PRIMARY="<<color->str()<<" --FONT_LABEL="<<color->str()<<" -B"<<axis<<"c"<<dataFileName;
    if(gridLine) ss<<" --MAP_GRID_PEN_PRIMARY="<<gridLine->str();
    if(orthoLabels)
      ss<<" --MAP_ANNOT_ORTHO=ns";
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GROOPS_REGISTER_CLASS(PlotAxis, "plotAxisType",
                      PlotAxisStandard,
                      PlotAxisTime,
                      PlotAxisLabeled)

GROOPS_READCONFIG_CLASS(PlotAxis, "plotAxisType")

/***********************************************/

PlotAxisPtr PlotAxis::create(Config &config, const std::string &name)
{
  try
  {
    PlotAxisPtr plotAxis;
    std::string  type;

    readConfigChoice(config, name, type, Config::MUSTSET, "", "axis limits, labels and annotation");
    if(readConfigChoiceElement(config, "standard", type, "generic x/y axis"))
      plotAxis = PlotAxisPtr(new PlotAxisStandard(config));
    if(readConfigChoiceElement(config, "time",     type, "intepret x-values as MJD"))
      plotAxis = PlotAxisPtr(new PlotAxisTime(config));
    if(readConfigChoiceElement(config, "labeled",  type, "major ticks with string labels"))
      plotAxis = PlotAxisPtr(new PlotAxisLabeled(config));
    endChoice(config);

    return plotAxis;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
