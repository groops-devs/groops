/***********************************************/
/**
* @file plotGraphLayer.cpp
*
* @brief Lines, points and polygons in 2d plots.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2016-07-23
*
*/
/***********************************************/

#define DOCSTRING_PlotGraphLayer

#include "base/import.h"
#include "parser/stringParser.h"
#include "parser/dataVariables.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "inputOutput/file.h"
#include "inputOutput/system.h"
#include "files/fileMatrix.h"
#include "classes/gravityfield/gravityfield.h"
#include "plot/plotMisc.h"
#include "plotGraphLayer.h"

/***********************************************/

static const char *docstringPlotGraphLayerLinesAndPoints = R"(
\subsection{LinesAndPoints}\label{plotGraphLayerType:linesAndPoints}
Draws a \configClass{line}{plotLineType} and/or points (\configClass{symbol}{plotSymbolType})
of xy data. The standard \reference{dataVariables}{general.parser:dataVariables}
are available to select the data columns of \configFile{inputfileMatrix}{matrix}.
If no \configClass{color}{plotColorType} of the \configClass{symbol}{plotSymbolType}
is given a \configClass{colorbar}{plotColorbarType}
is required and the color is determined by \config{valueZ}.
Additionally a vertical error bar can be plotted at each data point with
size \config{valueErrorBar}.

See \program{Gravityfield2AreaMeanTimeSeries} for an example plot.
)";

class PlotGraphLayerLinesAndPoints : public PlotGraphLayer
{
protected:
  std::pair<std::string, VariableList> description;
  PlotLinePtr   line;
  PlotSymbolPtr symbol;
  Bool          hasZValues, hasErrors;

public:
  PlotGraphLayerLinesAndPoints(Config &config);
  Bool requiresColorBar()   const override {return hasZValues;}
  std::string scriptEntry() const override;
  std::string legendEntry() const override;
};

/***********************************************/

PlotGraphLayerLinesAndPoints::PlotGraphLayerLinesAndPoints(Config &config)
{
  try
  {
    FileName fileName;
    ExpressionVariablePtr exprX, exprY, exprZ, exprError;

    readConfig(config, "inputfileMatrix",  fileName,     Config::MUSTSET,  "",      "each line contains x,y");
    readConfig(config, "valueX",           exprX,        Config::OPTIONAL, "data0", "expression for x-values (input columns are named data0, data1, ...)");
    readConfig(config, "valueY",           exprY,        Config::MUSTSET,  "data1", "expression for y-values (input columns are named data0, data1, ...)");
    readConfig(config, "valueZ",           exprZ,        Config::OPTIONAL, "",      "expression for the colorbar");
    readConfig(config, "valueErrorBar",    exprError,    Config::OPTIONAL, "",      "expression for error bars (input columns are named data0, data1, ...)");
    readConfig(config, "description",      description,  Config::OPTIONAL, "",      "text of the legend");
    readConfig(config, "line",             line,         Config::OPTIONAL, "solid", "");
    readConfig(config, "symbol",           symbol,       Config::OPTIONAL, "",      "");
    readConfig(config, "plotOnSecondAxis", onSecondAxis, Config::DEFAULT,  "0",     "draw dataset on a second Y-axis (if available).");
    if(isCreateSchema(config)) return;

    hasZValues = symbol && symbol->requiresColorBar() && exprZ;
    hasErrors  = (exprError != nullptr);

    // tests
    if(!line && !symbol)
      throw(Exception("At least one of line and symbol must be set."));
    if(symbol && symbol->requiresColorBar() && !exprZ)
      throw(Exception("valueZ is needed to determine color of line/symbol"));
    if(!hasZValues)
      exprZ = nullptr;

    // check if file exists
    // --------------------
    if(!System::exists(fileName))
    {
      if(description.first.empty())
        description.first = fileName.str();
      description.first += " (file not found)";
      logWarning<<"file <"<<fileName<<"> not found!"<<Log::endl;
      return;
    }

    // read data
    // ---------
    Matrix A;
    readFileMatrix(fileName, A);

    if(!A.size())
    {
      if(description.first.empty())
        description.first = fileName.str();
      description.first += " (empty file)";
      logWarning<<"file <"<<fileName<<"> is empty!"<<Log::endl;
      return;
    }

    // create data variables
    // ---------------------
    VariableList varList;
    addDataVariables(A, varList);
    try {description.second += varList; description.first = StringParser::parse(description.first, description.second);} catch(std::exception &) {}
    for(ExpressionVariablePtr expr : {exprX, exprY, exprZ, exprError})
      if(expr) expr->simplify(varList);

    // evaluate expressions
    // --------------------
    data = Matrix(A.rows(), 2+hasZValues+hasErrors);
    for(UInt i=0; i<A.rows(); i++)
    {
      UInt idx = 0;
      evaluateDataVariables(A, i, varList);
      if(!exprX)
        data(i, idx++) = static_cast<Double>(i); // default: index
      for(auto expression : {exprX, exprY, exprZ, exprError})
        if(expression)
          data(i, idx++) = expression->evaluate(varList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerLinesAndPoints::scriptEntry() const
{
  try
  {
    if(!data.size())
      return std::string();

    std::stringstream ss;
    if(hasErrors)
      ss<<"gmt psxy "<<dataFileName<<" -bi"<<data.columns()<<"d -J -R -i0,1,"<<(hasZValues ? 3 : 2)<<" -Sc1p -Ey/1p -O -K >> groopsPlot.ps"<<std::endl;
    if(line)
      ss<<"gmt psxy "<<dataFileName<<" -bi"<<data.columns()<<"d -J -R -W"<<line->str()<<" -O -K >> groopsPlot.ps"<<std::endl;
    if(symbol)
      ss<<"gmt psxy "<<dataFileName<<" -bi"<<data.columns()<<"d -J -R -S"<<symbol->str()<<" -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerLinesAndPoints::legendEntry() const
{
  try
  {
    std::stringstream ss;
    if(description.first.empty() || (!line && !symbol))
      return ss.str();
    if(line)
      ss<<"S 0.3c - 0.5c - "<<line->str()<<" 0.7c " ;
    if(symbol)
    {
      if(line)
        ss<<std::endl<<"G -1l"<<std::endl;
      ss<<"S 0.3c "<<symbol->legendStr()<<"\t0.7c\t";
    }
    ss<<description.first<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotGraphLayerErrorEnvelope = R"(
\subsection{ErrorEnvelope}
Draws a symmetrical envelope around \config{valueY} as function of \config{valueX}
using deviations \config{valueErrors}.
The standard \reference{dataVariables}{general.parser:dataVariables}
are available to select the data columns of \configFile{inputfileMatrix}{matrix}.
The data line itself is not plotted but must be added as extra
\configClass{layer:linesAndPoints}{plotGraphLayerType:linesAndPoints}.
)";

class PlotGraphLayerErrorEnvelope : public PlotGraphLayer
{
protected:
  std::pair<std::string, VariableList> description;
  PlotColorPtr fillColor;
  PlotLinePtr  edgeLine;

public:
  PlotGraphLayerErrorEnvelope(Config &config);
  std::string scriptEntry() const override;
  std::string legendEntry() const override;
};

/***********************************************/

PlotGraphLayerErrorEnvelope::PlotGraphLayerErrorEnvelope(Config &config)
{
  try
  {
    FileName fileName;
    ExpressionVariablePtr exprX, exprY, exprErrors;

    readConfig(config, "inputfileMatrix",  fileName,     Config::MUSTSET,  "",      "each line contains x,y");
    readConfig(config, "valueX",           exprX,        Config::OPTIONAL, "data0", "expression for x-values (input columns are named data0, data1, ...)");
    readConfig(config, "valueY",           exprY,        Config::MUSTSET,  "data1", "expression for y-values (input columns are named data0, data1, ...)");
    readConfig(config, "valueErrors",      exprErrors,   Config::MUSTSET,  "data2", "expression for error values");
    readConfig(config, "description",      description,  Config::OPTIONAL, "",      "text of the legend");
    readConfig(config, "fillColor",        fillColor,    Config::OPTIONAL, "gray",  "fill color of the envelope");
    readConfig(config, "edgeLine",         edgeLine,     Config::OPTIONAL, "",      "edge line style of the envelope");
    readConfig(config, "plotOnSecondAxis", onSecondAxis, Config::DEFAULT,  "0",     "draw dataset on a second Y-axis (if available).");
    if(isCreateSchema(config)) return;

    if(!fillColor && !edgeLine)
      throw(Exception("At least one of fillColor and edgeLine must be set."));

    // check if file exists
    // --------------------
    if(!System::exists(fileName))
    {
      if(description.first.empty())
        description.first = fileName.str();
      description.first += " (file not found)";
      logWarning<<"file <"<<fileName<<"> not found!"<<Log::endl;
      return;
    }

    // read data
    // ---------
    Matrix A;
    readFileMatrix(fileName, A);

    if(!A.size())
    {
      if(description.first.empty())
        description.first = fileName.str();
      description.first += " (empty file)";
      logWarning<<"file <"<<fileName<<"> is empty!"<<Log::endl;
      return;
    }

    // create data variables
    // ---------------------
    VariableList varList;
    addDataVariables(A, varList);
    try {description.second += varList; description.first = StringParser::parse(description.first, description.second);} catch(std::exception &) {}

    std::vector<ExpressionVariablePtr> expressions = {exprX, exprY, exprErrors};
    for(ExpressionVariablePtr expr : expressions)
      if(expr) expr->simplify(varList);

    // evaluate expressions
    // --------------------
    data = Matrix(A.rows(), 3);
    for(UInt i=0; i<A.rows(); i++)
    {
      evaluateDataVariables(A, i, varList);
      data(i, 0) = static_cast<Double>(i); // default: index
      for(UInt k=0; k<expressions.size(); k++)
        if(expressions.at(k))
          data(i, k) = expressions.at(k)->evaluate(varList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerErrorEnvelope::scriptEntry() const
{
  try
  {
    if(!data.size())
      return std::string();

    std::stringstream ss;
    ss<<"gmt psxy "<<dataFileName<<" -bi"<<data.columns()<<"d -J -R -L+d";
    if(edgeLine)
      ss<<"+p"<<edgeLine->str();
    if(fillColor)
      ss<<" -G"<<fillColor->str();
    ss<<" -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerErrorEnvelope::legendEntry() const
{
  try
  {
    if(description.first.empty())
      return std::string();

    std::stringstream ss;
    ss<<"S 0.3c s 0.5c "<<(fillColor ? fillColor->str() : "-"s)<<" ";
    if(edgeLine)
      ss<<edgeLine->str();
    else
      ss<<"0p,"<<fillColor->str();
    ss<<"\t0.7c\t"<<description.first<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotGraphLayerBars = R"(
\subsection{Bars}
Creates a bar plot with vertical or \config{horizontal} bars out of the given
x- and y-values. The standard \reference{dataVariables}{general.parser:dataVariables}
are available to select the data columns of \configFile{inputfileMatrix}{matrix}.
The bars ranges from \config{valueBase} (can be also an expression) to the \config{valueY}.
If no \configClass{color}{plotColorType} is given a \configClass{colorbar}{plotColorbarType}
is required and the color is determined by \config{valueZ}.

See \program{Instrument2Histogram} for an example plot.
)";

class PlotGraphLayerBars : public PlotGraphLayer
{
private:
  Double       barWidth;
  Bool         horizontal;
  PlotColorPtr color;
  PlotLinePtr  edgeLine;
  std::pair<std::string, VariableList> description;

public:
  PlotGraphLayerBars(Config &config);
  std::string scriptEntry() const override;
  std::string legendEntry() const override;
};

/***********************************************/

PlotGraphLayerBars::PlotGraphLayerBars(Config &config)
{
  try
  {
    FileName fileName;
    ExpressionVariablePtr exprX, exprY, exprZ, exprBase, exprWidth;

    horizontal = FALSE;

    readConfig(config, "inputfileMatrix",  fileName,     Config::MUSTSET,  "",      "each line contains x,y");
    readConfig(config, "valueX",           exprX,        Config::OPTIONAL, "data0", "expression for x-values (input columns are named data0, data1, ...)");
    readConfig(config, "valueY",           exprY,        Config::MUSTSET,  "data1", "expression for y-values (input columns are named data0, data1, ...)");
    readConfig(config, "valueZ",           exprZ,        Config::OPTIONAL, "",      "expression for the colorbar");
    readConfig(config, "valueBase",        exprBase,     Config::OPTIONAL, "",      "base value of bars (default: minimum y-value)");
    readConfig(config, "width",            exprWidth,    Config::OPTIONAL, "",      "width of bars (default: minimum x-gap)");
    readConfig(config, "horizontal",       horizontal,   Config::OPTIONAL, "",      "draw horizontal bars instead of vertical");
    readConfig(config, "description",      description,  Config::OPTIONAL, "",      "text of the legend");
    readConfig(config, "color",            color,        Config::OPTIONAL, "black", "");
    readConfig(config, "edgeLine",         edgeLine,     Config::OPTIONAL, "",      "line");
    readConfig(config, "plotOnSecondAxis", onSecondAxis, Config::DEFAULT,  "0",     "draw dataset on a second Y-axis (if available).");
    if(isCreateSchema(config)) return;

    // check if file exists
    // --------------------
    if(!System::exists(fileName))
    {
      if(description.first.empty())
        description.first = fileName.str();
      description.first += " (file not found)";
      logWarning<<"file <"<<fileName<<"> not found!"<<Log::endl;
      return;
    }

    // read data
    // ---------
    Matrix A;
    readFileMatrix(fileName, A);

    if(!A.size())
    {
      if(description.first.empty())
        description.first = fileName.str();
      description.first += " (empty file)";
      logWarning<<"file <"<<fileName<<"> is empty!"<<Log::endl;
      return;
    }

    // create data variables
    // ---------------------
    VariableList varList;
    addDataVariables(A, varList);
    try {description.second += varList; description.first = StringParser::parse(description.first, description.second);} catch(std::exception &) {}
    for(ExpressionVariablePtr expr : {exprX, exprY, exprZ, exprBase, exprBase, exprWidth})
      if(expr) expr->simplify(varList);

    // evaluate expressions
    // --------------------
    data = Matrix(A.rows(), exprZ ? 4 : 3);
    for(UInt i=0; i<A.rows(); i++)
    {
      evaluateDataVariables(A, i, varList);
      data(i, 0) = static_cast<Double>(i); // default: index
      if(exprX)    data(i, 0) = exprX->evaluate(varList);
      if(exprY)    data(i, 1) = exprY->evaluate(varList);
      if(exprZ)    data(i, 2) = exprZ->evaluate(varList);
      if(exprBase) data(i, exprZ ? 3 : 2) = exprBase->evaluate(varList);
    }

    // compute bar width
    // -----------------
    if(!exprWidth)
    {
      barWidth = std::numeric_limits<Double>::infinity();
      for(UInt i=1; i<data.rows(); i++)
        barWidth = std::min(barWidth, std::abs(data(i, horizontal ? 1 : 0) - (data(i-1, horizontal ? 1 : 0))));
    }
    else
      barWidth = exprWidth->evaluate(varList);

    // compute base value
    // ------------------
    if(!exprBase)
      copy(Vector(data.rows(), min(data.column(horizontal ? 0 : 1))), data.column(2 +  (exprZ ? 1 : 0)));

    // test z-values
    // -------------
    if(!color && !exprZ)
      throw(Exception("valueZ is needed to determine color of line/symbol"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerBars::scriptEntry() const
{
  try
  {
    if(!data.size())
      return std::string();

    std::stringstream ss;
    ss<<"gmt psxy "<<dataFileName<<" -bi"<<data.columns()<<"d -J -R -S"<<(horizontal ? "B" : "b")<<barWidth<<"ub";
    if(edgeLine)
      ss<<" -W"<<edgeLine->str();
    if(color)
      ss<<" -G"<<color->str();
    else
      ss<<" -CgroopsPlot.cpt";
    ss<<" -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerBars::legendEntry() const
{
  try
  {
    if(description.first.empty())
      return std::string();

    if(!color)
    {
      logWarning<<"No legend entry can be generated for bar layers when using <fromValue> as color."<<Log::endl;
      return std::string();
    }

    std::stringstream ss;
    ss<<"S 0.3c r 5p "<<color->str()<<" ";
    if(edgeLine) ss<<edgeLine->str();
    else ss<<"-";
    ss<<"\t0.7c\t"<<description.first<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotGraphLayerGridded = R"(
\subsection{Gridded}
Creates a regular grid of yxz values. The standard \reference{dataVariables}{general.parser:dataVariables}
are available to select the data columns of \configFile{inputfileMatrix}{matrix}.
Empty grid cells are not plotted. Cells with more than one value will be set to the mean value.
The grid spacing is determined by the median spacing of the input data or set by \config{incrementX/Y}.

See \program{Orbit2ArgumentOfLatitude} for an example plot.
)";

/***** CLASS ***********************************/

class PlotGraphLayerGridded : public PlotGraphLayer
{
  Double incX, incY;

public:
  PlotGraphLayerGridded(Config &config);
  Bool requiresColorBar() const override {return TRUE;}
  std::string scriptEntry() const override;
  Double bufferX() const override {return 0.5 * incX;}
  Double bufferY() const override {return 0.5 * incY;}
};

/***********************************************/

PlotGraphLayerGridded::PlotGraphLayerGridded(Config &config)
{
  try
  {
    FileName    fileName;
    std::string choice;
    ExpressionVariablePtr exprX, exprY, exprZ;
    incX = incY = NAN_EXPR;

    readConfig(config, "inputfileMatrix",  fileName,     Config::MUSTSET,   "",      "each line contains x,y,z");
    readConfig(config, "valueX",           exprX,        Config::OPTIONAL,  "data0", "expression for x-values (input columns are named data0, data1, ...)");
    readConfig(config, "valueY",           exprY,        Config::MUSTSET,   "data1", "expression for y-values (input columns are named data0, data1, ...)");
    readConfig(config, "valueZ",           exprZ,        Config::MUSTSET,   "data2", "expression for the colorbar");
    readConfig(config, "incrementX",       incX,         Config::OPTIONAL,  "",      "the grid spacing");
    readConfig(config, "incrementY",       incY,         Config::OPTIONAL,  "",      "the grid spacing");
    readConfig(config, "plotOnSecondAxis", onSecondAxis, Config::DEFAULT,   "0",     "draw dataset on a second Y-axis (if available).");
    if(isCreateSchema(config)) return;

    // check if file exists
    // --------------------
    if(!System::exists(fileName))
    {
      logWarning<<"file <"<<fileName<<"> not found!"<<Log::endl;
      return;
    }

    // read data
    // ---------
    Matrix A;
    readFileMatrix(fileName, A);

    if(!A.size())
    {
      logWarning<<"file <"<<fileName<<"> is empty!"<<Log::endl;
      return;
    }

    // create data variables
    // ---------------------
    VariableList varList;
    addDataVariables(A, varList);

    std::vector<ExpressionVariablePtr> expressions = {exprX, exprY, exprZ};
    for(ExpressionVariablePtr expr : expressions)
      if(expr) expr->simplify(varList);

    // evaluate expressions
    // --------------------
    data = Matrix(A.rows(), (exprZ ? 3 : 2));
    for(UInt i=0; i<A.rows(); i++)
    {
      evaluateDataVariables(A, i, varList);
      data(i, 0) = static_cast<Double>(i); // default: index
      for(UInt k=0; k<expressions.size(); k++)
        if(expressions.at(k))
          data(i, k) = expressions.at(k)->evaluate(varList);
    }

    // determine sampling (median)
    // ---------------------------
    auto computeIncrement = [](std::vector<Double> x)
    {
      std::sort(x.begin(), x.end());
      auto it = std::unique(x.begin(), x.end());
      std::vector<Double> dx;
      std::adjacent_difference(x.begin(), it, std::back_inserter(dx));
      std::nth_element(dx.begin(), dx.begin()+dx.size()/2, dx.end());
      return dx[dx.size()/2];
    };

    if(std::isnan(incX))
      incX = computeIncrement(Vector(data.column(0)));
    if(std::isnan(incY))
      incY = computeIncrement(Vector(data.column(1)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerGridded::scriptEntry() const
{
  try
  {
    if(!data.size())
      return std::string();

    const Double minX = min(data.column(0)) - bufferX();
    const Double maxX = max(data.column(0)) + bufferX();
    const Double minY = min(data.column(1)) - bufferY();
    const Double maxY = max(data.column(1)) + bufferY();

    std::stringstream ss;
    ss<<"gmt xyz2grd "<<dataFileName<<" -bi3d -G"<<dataFileName<<".grd -Vn";
    ss<<" -r -I"<<incX<<"=/"<<incY<<"= -R"<<minX<<"/"<<maxX<<"/"<<minY<<"/"<<maxY<<std::endl;
    ss<<"gmt grdimage "<<dataFileName<<".grd -Q -J -R"<<PlotBasics::scriptVariable("range")<<" -CgroopsPlot.cpt";
    ss<<" -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotGraphLayerRectangle = R"(
\subsection{Rectangle}
Plots a rectangle to highlight an area.
)";

/***** CLASS ***********************************/

class PlotGraphLayerRectangle : public PlotGraphLayer
{
  PlotLinePtr  edgeLine;
  PlotColorPtr fillColor;
  std::string  description;

public:
  PlotGraphLayerRectangle(Config &config);
  void        writeDataFile(const FileName &workingDirectory, UInt idxLayer, Double minX, Double maxX, Double minY, Double maxY) override;
  std::string scriptEntry() const override;
  std::string legendEntry() const override;
};

/***********************************************/

PlotGraphLayerRectangle::PlotGraphLayerRectangle(Config &config)
{
  try
  {
    Double x1=NAN_EXPR, y1=NAN_EXPR, x2=NAN_EXPR, y2=NAN_EXPR;

    readConfig(config, "minX",             x1,           Config::OPTIONAL, "",  "empty: left");
    readConfig(config, "maxX",             x2,           Config::OPTIONAL, "",  "empty: right");
    readConfig(config, "minY",             y1,           Config::OPTIONAL, "",  "empty: bottom");
    readConfig(config, "maxY",             y2,           Config::OPTIONAL, "",  "empty: top");
    readConfig(config, "description",      description,  Config::OPTIONAL, "",  "text of the legend");
    readConfig(config, "edgeLine",         edgeLine,     Config::OPTIONAL, "",  "");
    readConfig(config, "fillColor",        fillColor,    Config::OPTIONAL, "",  "");
    readConfig(config, "plotOnSecondAxis", onSecondAxis, Config::DEFAULT,  "0", "draw dataset on a second Y-axis (if available).");
    if(isCreateSchema(config)) return;

    if(!fillColor && !edgeLine)
      throw(Exception("At least one of fillColor and edgeLine must be set."));

    data = Matrix(4, 2);
    copy(Vector({x1, x1, x2, x2}), data.column(0));
    copy(Vector({y1, y2, y2, y1}), data.column(1));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotGraphLayerRectangle::writeDataFile(const FileName &workingDirectory, UInt idxLayer, Double minX, Double maxX, Double minY, Double maxY)
{
  try
  {
    if(std::isnan(data(0, 0))) data(0, 0) = data(1, 0) = minX;
    if(std::isnan(data(2, 0))) data(2, 0) = data(3, 0) = maxX;
    if(std::isnan(data(0, 1))) data(0, 1) = data(3, 1) = minY;
    if(std::isnan(data(1, 1))) data(1, 1) = data(2, 1) = maxY;
    PlotGraphLayer::writeDataFile(workingDirectory, idxLayer, minX, maxX, minY, maxY);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerRectangle::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt psclip "<<dataFileName<<" -bi2d -J -R -O -K >> groopsPlot.ps"<<std::endl;
    ss<<"gmt psxy "<<dataFileName<<" -bi2d -A -L -J -R";
    if(fillColor) ss<<" -G"<<fillColor->str();
    if(edgeLine) ss<<" -W"<<edgeLine->str();
    ss<<" -O -K >> groopsPlot.ps"<<std::endl;
    ss<<"gmt psclip -C -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerRectangle::legendEntry() const
{
  try
  {
    std::stringstream ss;
    if(!description.empty())
    {
      ss<<"S 0.3c s 0.5c "<<(fillColor ? fillColor->str() : "-")<<" ";
      if(edgeLine)
        ss<<edgeLine->str();
      else
        ss<<"0p,"<<fillColor->str();
      ss<<"\t0.7c\t"<<description<<std::endl;
    }
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotGraphLayerText = R"(
\subsection{Text}
Writes a \config{text} at \config{originX} and \config{originY} position in the graph.
With \config{clip} the text is cutted at the boundaries of the plotting area.
)";

class PlotGraphLayerText : public PlotGraphLayer
{
  Double       xOffset, yOffset;
  Double       fontSize;
  PlotColorPtr fontColor;
  std::string  text;
  std::string  alignment;
  Bool         clip;

public:
  PlotGraphLayerText(Config &config);
  void writeDataFile(const FileName &, UInt, Double, Double, Double, Double) override;
  std::string scriptEntry() const override;
};

/***********************************************/

PlotGraphLayerText::PlotGraphLayerText(Config &config)
{
  try
  {
    data = Matrix(1, 2);
    readConfig(config, "text",             text,         Config::MUSTSET, "",   "");
    readConfig(config, "originX",          data(0, 0),   Config::MUSTSET, "",   "");
    readConfig(config, "originY",          data(0, 1),   Config::MUSTSET, "",   "");
    readConfig(config, "offsetX",          xOffset,      Config::DEFAULT, "0",  "[cm] x-offset from origin");
    readConfig(config, "offsetY",          yOffset,      Config::DEFAULT, "0",  "[cm] y-offset from origin");
    readConfig(config, "alignment",        alignment,    Config::DEFAULT, "BL", "L, C, R (left, center, right) and T, M, B (top, middle, bottom)");
    readConfig(config, "fontSize",         fontSize,     Config::DEFAULT, "10", "[pt]");
    readConfig(config, "fontColor",        fontColor,    Config::MUSTSET, "",   "");
    readConfig(config, "clip",             clip,         Config::DEFAULT, "1",  "clip at boundaries");
    readConfig(config, "plotOnSecondAxis", onSecondAxis, Config::DEFAULT, "0",  "draw dataset on a second Y-axis (if available).");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotGraphLayerText::writeDataFile(const FileName &workingDirectory, UInt idxLayer, Double /*minX*/, Double /*maxX*/, Double /*minY*/, Double /*maxY*/)
{
  try
  {
    dataFileName = "text."+idxLayer%"%i.txt"s;
    OutFile file(workingDirectory.append(dataFileName));
    file<<data(0, 0)<<" "<<data(0, 1)<<" "<<text;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerText::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt pstext "<<dataFileName<<" -F+f"<<fontSize<<"p,,"<<fontColor->str()<<"+j"<<alignment<<" -D"<<xOffset<<"/"<<yOffset<<" -J -R "<<(clip ? "" : "-N")<<" -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotGraphLayerDegreeAmplitudes = R"(
\subsection{DegreeAmplitudes}
Plot degree amplitudes of potential coefficients computed by \program{Gravityfield2DegreeAmplitudes}
or \program{PotentialCoefficients2DegreeAmplitudes}.
The standard \reference{dataVariables}{general.parser:dataVariables} are available to select
the data columns of \configFile{inputfileMatrix}{matrix}. It plots a solid line for the
\config{valueSignal} and a dotted line for the \config{valueError} per default.
)";

class PlotGraphLayerDegreeAmplitudes : public PlotGraphLayer
{
  std::pair<std::string, VariableList> description;
  PlotLinePtr lineSignal, lineErrors;

public:
  PlotGraphLayerDegreeAmplitudes(Config &config);
  std::string scriptEntry() const override;
  std::string legendEntry() const override;
};

/***********************************************/

PlotGraphLayerDegreeAmplitudes::PlotGraphLayerDegreeAmplitudes(Config &config)
{
  try
  {
    FileName fileName;
    ExpressionVariablePtr exprDegree, exprSignal, exprErrors;

    readConfig(config, "inputfileMatrix",  fileName,     Config::MUSTSET,  "",      "degree amplitudes");
    readConfig(config, "valueDegree",      exprDegree,   Config::OPTIONAL, "data0", "expression for x-values (degrees) (input columns are named data0, data1, ...)");
    readConfig(config, "valueSignal",      exprSignal,   Config::OPTIONAL, "data1", "expression for y-values (signal) (input columns are named data0, data1, ...)");
    readConfig(config, "valueErrors",      exprErrors,   Config::OPTIONAL, "data2", "expression for y-values (formal errors)");
    readConfig(config, "description",      description,  Config::OPTIONAL, "",      "text of the legend");
    readConfig(config, "lineSignal",       lineSignal,   Config::OPTIONAL, "solid", "");
    readConfig(config, "lineErrors",       lineErrors,   Config::OPTIONAL, R"({"custom": {"style":"5_2:0"}})", "");
    readConfig(config, "plotOnSecondAxis", onSecondAxis, Config::DEFAULT,  "0",     "draw dataset on a second Y-axis (if available).");
    if(isCreateSchema(config)) return;

    // check if file exists
    // --------------------
    if(!System::exists(fileName))
    {
      if(description.first.empty())
        description.first = fileName.str();
      description.first += " (file not found)";
      logWarning<<"file <"<<fileName<<"> not found!"<<Log::endl;
      return;
    }

    // read data
    // ---------
    Matrix A;
    readFileMatrix(fileName, A);

    if(!A.size())
    {
      if(description.first.empty())
        description.first = fileName.str();
      description.first += " (empty file)";
      logWarning<<"file <"<<fileName<<"> is empty!"<<Log::endl;
      return;
    }

    // create data variables
    // ---------------------
    VariableList varList;
    addDataVariables(A, varList);
    try {description.second += varList; description.first = StringParser::parse(description.first, description.second);} catch(std::exception &) {}

    std::vector<ExpressionVariablePtr> expressions = {exprDegree, exprSignal, exprErrors};
    for(ExpressionVariablePtr expr : expressions)
      if(expr) expr->simplify(varList);

    // evaluate expressions
    // --------------------
    data = Matrix(A.rows(), 3, NAN_EXPR);
    for(UInt i=0; i<A.rows(); i++)
    {
      evaluateDataVariables(A, i, varList);
      data(i, 0) = static_cast<Double>(i); // default: index
      for(UInt k=0; k<expressions.size(); k++)
        if(expressions.at(k))
          data(i, k) = expressions.at(k)->evaluate(varList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerDegreeAmplitudes::scriptEntry() const
{
  try
  {
    if(!data.size())
      return std::string();

    std::stringstream ss;
    if(lineSignal)
      ss<<"gmt psxy "<<dataFileName<<" -bi"<<data.columns()<<"d -J -R -i0,1 -W"<<lineSignal->str()<<" -O -K >> groopsPlot.ps"<<std::endl;
    if(lineErrors)
      ss<<"gmt psxy "<<dataFileName<<" -bi"<<data.columns()<<"d -J -R -i0,2 -W"<<lineErrors->str()<<" -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotGraphLayerDegreeAmplitudes::legendEntry() const
{
  try
  {
    if(description.empty() || (!lineErrors && !lineSignal))
      return std::string();

    std::stringstream ss;
    ss<<"S 0.3c - 0.5c - "<<(lineSignal ? lineSignal->str() : lineErrors->str())<<" 0.7c "<<description<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GROOPS_REGISTER_CLASS(PlotGraphLayer, "plotGraphLayerType",
                      PlotGraphLayerLinesAndPoints,
                      PlotGraphLayerErrorEnvelope,
                      PlotGraphLayerBars,
                      PlotGraphLayerGridded,
                      PlotGraphLayerRectangle,
                      PlotGraphLayerText,
                      PlotGraphLayerDegreeAmplitudes)

GROOPS_READCONFIG_CLASS(PlotGraphLayer, "plotGraphLayerType")

/***********************************************/

PlotGraphLayerPtr PlotGraphLayer::create(Config &config, const std::string &name)
{
  try
  {
    PlotGraphLayerPtr plotGraphLayer;
    std::string  type;

    readConfigChoice(config, name, type, Config::MUSTSET, "", "lines, points and polygons");
    if(readConfigChoiceElement(config, "linesAndPoints",   type, "line/points"))
      plotGraphLayer = PlotGraphLayerPtr(new PlotGraphLayerLinesAndPoints(config));
    if(readConfigChoiceElement(config, "errorEnvelope",    type, "error envelope for line plots"))
      plotGraphLayer = PlotGraphLayerPtr(new PlotGraphLayerErrorEnvelope(config));
    if(readConfigChoiceElement(config, "bars",             type, "bar graph"))
      plotGraphLayer = PlotGraphLayerPtr(new PlotGraphLayerBars(config));
    if(readConfigChoiceElement(config, "gridded",          type, "mesh data"))
      plotGraphLayer = PlotGraphLayerPtr(new PlotGraphLayerGridded(config));
    if(readConfigChoiceElement(config, "rectangle",        type, "draw rectangle"))
      plotGraphLayer = PlotGraphLayerPtr(new PlotGraphLayerRectangle(config));
    if(readConfigChoiceElement(config, "text",             type, "text"))
      plotGraphLayer = PlotGraphLayerPtr(new PlotGraphLayerText(config));
    if(readConfigChoiceElement(config, "degreeAmplitudes", type, "degree amplitudes of a gravity field"))
      plotGraphLayer = PlotGraphLayerPtr(new PlotGraphLayerDegreeAmplitudes(config));
    endChoice(config);

    return plotGraphLayer;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotGraphLayer::getIntervalX(Bool isLogarithmic, Double &minX, Double &maxX) const
{
  try
  {
    for(UInt i=0; i<data.rows(); i++)
      if(!std::isnan(data(i, 0)) && (!isLogarithmic || (data(i, 0) > 0)))
      {
        minX = std::min(minX, data(i, 0)-bufferX());
        maxX = std::max(maxX, data(i, 0)+bufferX());
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotGraphLayer::getIntervalY(Bool isLogarithmic, Double minX, Double maxX, Double &minY, Double &maxY) const
{
  try
  {
    for(UInt i=0; i<data.rows(); i++)
      if((minX <= data(i, 0)) && (data(i, 0) <= maxX))
        if(!std::isnan(data(i, 1)) && (!isLogarithmic || (data(i, 1) > 0)))
        {
          minY = std::min(minY, data(i, 1)-bufferY());
          maxY = std::max(maxY, data(i, 1)+bufferY());
        }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotGraphLayer::getIntervalZ(Bool isLogarithmic, Double minX, Double maxX, Double minY, Double maxY, Double &minZ, Double &maxZ) const
{
  try
  {
    if(!requiresColorBar())
      return;

    UInt   count =  0;
    Double avg   =  0.;
    for(UInt i=0; i<data.rows(); i++)
      if((minX <= data(i, 0)) && (data(i, 0) <= maxX) && (minY <= data(i, 1)) && (data(i, 1) <= maxY))
        if(!std::isnan(data(i, 2)) && (!isLogarithmic || (data(i, 2) > 0)))
        {
          minZ  = std::min(minZ, data(i, 2));
          maxZ  = std::max(maxZ, data(i, 2));
          avg  += std::fabs(data(i, 2));
          count++;
        }
    avg /= count;

    if(!isLogarithmic)
    {
      minZ = (minZ > 0) ? 0 : -3*avg;
      maxZ = +3*avg;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotGraphLayer::writeDataFile(const FileName &workingDirectory, UInt idxLayer, Double /*minX*/, Double /*maxX*/, Double /*minY*/, Double /*maxY*/)
{
  try
  {
    dataFileName = "data."+idxLayer%"%i.dat"s;
    OutFile file(workingDirectory.append(dataFileName), std::ios::out | std::ios::binary);
    for(UInt i=0; i<data.rows(); i++)
      for(UInt k=0; k<data.columns(); k++)
        file.write(reinterpret_cast<char *>(&data(i, k)), sizeof(Double));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
