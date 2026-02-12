/***********************************************/
/**
* @file plotMapLayer.cpp
*
* @brief plot layers of maps.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2015-10-23
*
*/
/***********************************************/

#define DOCSTRING_PlotMapLayer

#include "base/import.h"
#include "parser/dataVariables.h"
#include "config/configRegister.h"
#include "inputOutput/system.h"
#include "files/filePolygon.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "misc/miscGriddedData.h"
#include "plot/plotMisc.h"
#include "plotMapLayer.h"

/***********************************************/

// Latex documentation
static const char *docstringPlotMapLayerGrid = R"(
\subsection{GriddedData}
Creates a regular grid of xyz values. The standard \reference{dataVariables}{general.parser:dataVariables}
are available to select the data column of \configFile{inputfileGriddedData}{griddedData}.
Empty grid cells are not plotted. Cells with more than one value will be set to the mean value.
The grid spacing can be determined automatically for regular rectangular grids otherwise
it must be set with \config{increment}. To get a better display together with some projections
the grid should be internally \config{resample}d to higher resolution.
It is assumed that the points of \configFile{inputfileGriddedData}{griddedData} represents centers of grid cells.
This assumption can be changed with \config{gridlineRegistered} (e.g. if the data starts at the north pole).
)";

class PlotMapLayerGrid : public PlotMapLayer
{
  Angle  incrementLat, incrementLon;
  Bool   illuminate;
  Bool   isGridline;
  Double intermediateDpi, threshold;
  Char   interpolationMethod;
  Bool   resample;

public:
  PlotMapLayerGrid(Config &config);
  Bool        requiresColorBar() const override {return TRUE;}
  std::string scriptEntry() const override;
};

/***********************************************/

PlotMapLayerGrid::PlotMapLayerGrid(Config &config)
{
  try
  {
    FileName              fileNameGrid;
    ExpressionVariablePtr exprValue;
    std::string           choice;

    resample = FALSE;

    readConfig(config, "inputfileGriddedData", fileNameGrid, Config::MUSTSET,  "",      "");
    readConfig(config, "value",                exprValue,    Config::MUSTSET,  "data0", "expression to compute values (input columns are named data0, data1, ...)");
    readConfig(config, "increment",            incrementLon, Config::OPTIONAL, "",      "the grid spacing [degrees]");
    readConfig(config, "illuminate",           illuminate,   Config::DEFAULT,  "0",     "illuminate grid");
    if(readConfigSequence(config, "resample", Config::OPTIONAL, "", ""))
    {
      resample = TRUE;
      readConfig(config, "intermediateDpi", intermediateDpi, Config::DEFAULT, "100", "oversample grid for a smoother visual effect");
      if(readConfigChoice(config, "interpolationMethod",  choice, Config::MUSTSET, "", "interpolation method for oversampling"))
      {
        if(readConfigChoiceElement(config, "bspline",  choice, "B-Spline interpolation"))             interpolationMethod = 'b';
        if(readConfigChoiceElement(config, "bicubic",  choice, "bicubic interpolation"))              interpolationMethod = 'c';
        if(readConfigChoiceElement(config, "bilinear", choice, "bilinear interpolation"))             interpolationMethod = 'l';
        if(readConfigChoiceElement(config, "nearest",  choice, "nearest neighbour interpolation"))    interpolationMethod = 'n';
        endChoice(config);
      }
      readConfig(config, "threshold", threshold, Config::DEFAULT, "0.5", "A threshold of 1.0 requires all (4 or 16) nodes involved in interpolation to be non-NaN. 0.5 will interpolate about half way from a non-NaN value; 0.1 will go about 90% of the way.");
      endSequence(config);
    }
    readConfig(config, "gridlineRegistered", isGridline, Config::DEFAULT,  "0", "treat input as point values instead of cell means");
    if(isCreateSchema(config)) return;

    GriddedData grid;
    readFileGriddedData(fileNameGrid, grid);
    points = grid.points;
    areas  = grid.areas;

    if(grid.values.size()==0)
      throw(Exception("<"+fileNameGrid.str()+"> has no values."));

    // try to define grid spacing
    // --------------------------
    incrementLat = incrementLon;
    if(incrementLon <= 0)
    {
      std::vector<Angle>  lambda, phi;
      std::vector<Double> radius;
      if(!grid.isRectangle(lambda, phi, radius))
         throw(Exception("'increment' must be set for non rectangular grids"));

      Vector dLambda(lambda.size()-1);
      for(UInt k=0; k<dLambda.size(); k++)
        dLambda(k) = std::fabs(std::remainder(lambda.at(k+1)-lambda.at(k), 2*PI));
      incrementLon = mean(dLambda);

      Angle  lon;
      Double h;
      for(UInt i=0; i<phi.size(); i++)
        grid.ellipsoid(polar(Angle(0.), phi.at(i), radius.at(i)), lon, phi.at(i), h);  // geocentric -> ellipsoidal

      Vector dPhi(phi.size()-1);
      for(UInt i=0; i<dPhi.size(); i++)
        dPhi(i) = std::fabs(phi.at(i+1)-phi.at(i));
      incrementLat = mean(dPhi);
    }

    // evaluate expression
    // -------------------
    VariableList varList;
    addDataVariables(grid, varList);
    exprValue->simplify(varList);

    data = Matrix(points.size(), 1);
    for(UInt i=0; i<points.size(); i++)
    {
      evaluateDataVariables(grid, i, varList);
      data(i, 0) = exprValue->evaluate(varList);
    }

    if(!isGridline)
    {
      bufferLon = 0.5 * incrementLon;
      bufferLat = 0.5 * incrementLat;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotMapLayerGrid::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt xyz2grd -bi3d "<<dataFileName<<" -G"<<dataFileName<<".grd -Vn -I"<<incrementLon*RAD2DEG*3600.<<"s/"<<incrementLat*RAD2DEG*3600.<<"s -R"<<PlotBasics::scriptVariable("region")<<(isGridline ? "" : " -r")<<std::endl;
    if(illuminate)
      ss<<"gmt grdgradient "<<dataFileName<<".grd -Nt0.8 -A45/315 -G"<<dataFileName<<".intense.grd"<<std::endl;
    ss<<"gmt grdimage "<<dataFileName<<".grd";
    if(resample)
      ss<<" -E"<<intermediateDpi<<" -n"<<interpolationMethod<<"+bg+t"<<threshold;
    ss<<" -Q -J -R -B -CgroopsPlot.cpt";
    if(illuminate)
      ss<<" -I"<<dataFileName<<".intense.grd";
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

// Latex documentation
static const char *docstringPlotMapLayerPoints = R"(
\subsection{Points}\label{plotMapLayerType:points}
Draws points (\configClass{symbol}{plotSymbolType}) and/or \configClass{line}{plotLineType}s
between the points. If no \configClass{color}{plotColorType} of the \configClass{symbol}{plotSymbolType}
is given a \configClass{colorbar}{plotColorbarType} is required and the color is determined
by the \config{value} expression. The standard \reference{dataVariables}{general.parser:dataVariables}
are available to select the data column of \configFile{inputfileGriddedData}{griddedData}.
)";

class PlotMapLayerPoints : public PlotMapLayer
{
  PlotLinePtr   line;
  PlotSymbolPtr symbol;
  Bool greatCircle;

public:
  PlotMapLayerPoints(Config &config);
  Bool        requiresColorBar() const override {return data.columns();}
  std::string scriptEntry() const override;
};

/***********************************************/

PlotMapLayerPoints::PlotMapLayerPoints(Config &config)
{
  try
  {
    FileName              fileNameGrid;
    ExpressionVariablePtr exprValue;

    readConfig(config, "inputfileGriddedData", fileNameGrid, Config::MUSTSET,  "",      "");
    readConfig(config, "value",                exprValue,    Config::OPTIONAL, "data0", "expression to compute color (input columns are named data0, data1, ...)");
    readConfig(config, "symbol",               symbol,       Config::OPTIONAL, "1",     "");
    readConfig(config, "line",                 line,         Config::OPTIONAL, "",      "style of connecting lines");
    readConfig(config, "drawLineAsGreatCircle",greatCircle,  Config::DEFAULT,  "1",     "draw connecting lines as great circles (otherwise, a straight line is drawn instead)");
    if(isCreateSchema(config)) return;

    // tests
    if(!line && !symbol)
      throw(Exception("At least one of line and symbol must be set."));
    if(symbol && symbol->requiresColorBar() && !exprValue)
      throw(Exception("value is needed to determine color of line/symbol"));
    if(!symbol || !symbol->requiresColorBar())
      exprValue = nullptr;

    GriddedData grid;
    readFileGriddedData(fileNameGrid, grid);
    points = grid.points;

    // evaluate expression
    // -------------------
    if(exprValue)
    {
      VariableList varList;
      addDataVariables(grid, varList);
      exprValue->simplify(varList);

      data = Matrix(points.size(), 1);
      for(UInt i=0; i<points.size(); i++)
      {
        evaluateDataVariables(grid, i, varList);
        data(i, 0) = exprValue->evaluate(varList);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotMapLayerPoints::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    if(line)
    {
      ss<<"gmt psxy -bi"<<2+data.columns()<<"d "<<dataFileName<<" -J -R -W"<<line->str();
      if(!greatCircle) ss<<" -A";
      ss<<" -O -K >> groopsPlot.ps"<<std::endl;
    }
    if(symbol)
      ss<<"gmt psxy -bi"<<2+data.columns()<<"d "<<dataFileName<<" -J -R -S"<<symbol->str()<<" -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static const char *docstringPlotMapLayerArrows = R"(
\subsection{Arrows}
Draws an arrow for each point in \configFile{inputfileGriddedData}{griddedData}.
The arrows are defined by the expressions \config{valueNorth/East}.
The standard \reference{dataVariables}{general.parser:dataVariables}
are available to select the correspondent data columns of \configFile{inputfileGriddedData}{griddedData}.
The \config{scale} factor converts the input units to cm in the plot.
If no \configClass{color}{plotColorType} is given a \configClass{colorbar}{plotColorbarType} is required
and the color is determined by the \config{value} expression.
With \config{scaleArrow} a reference arrow as legend can be plotted inside or outside the map.
)";

class PlotMapLayerArrows : public PlotMapLayer
{
  Double       scale;
  Double       penSize;
  Double       headSize;
  PlotColorPtr penColor;
  Double       scaleArrowOriginX;
  Double       scaleArrowOriginY;
  Double       scaleArrowLength;
  std::string  scaleArrowUnit;
  std::string  scaleArrowLabel;
  Bool         drawScaleArrow;
  Bool         hasZValues;

public:
  PlotMapLayerArrows(Config &config);
  Bool        requiresColorBar() const override {return hasZValues;}
  std::string scriptEntry() const override;
  std::string legendEntry(const FileName &workingDirectory, UInt idxLayer) const override;
};

/***********************************************/

PlotMapLayerArrows::PlotMapLayerArrows(Config &config)
{
  try
  {
    drawScaleArrow = FALSE;

    FileName fileNameGrid;
    ExpressionVariablePtr exprValueNorth, exprValueEast, exprValue;

    readConfig(config, "inputfileGriddedData", fileNameGrid,   Config::MUSTSET,  "",      "grid file with north and east values for arrows");
    readConfig(config, "valueNorth",           exprValueNorth, Config::MUSTSET,  "data0", "expression to compute north values (input columns are named data0, data1, ...)");
    readConfig(config, "valueEast",            exprValueEast,  Config::MUSTSET,  "data1", "expression to compute east values (input columns are named data0, data1, ...)");
    readConfig(config, "value",                exprValue,      Config::OPTIONAL, "",      "expression to compute arrow color (input columns are named data0, data1, ...)");
    readConfig(config, "scale",                scale,          Config::DEFAULT,  "50",    "[cm per input unit] length scale factor");
    readConfig(config, "penSize",              penSize,        Config::MUSTSET,  "1",     "[pt] width of arrow shaft");
    readConfig(config, "headSize",             headSize,       Config::MUSTSET,  "5",     "[pt] size of arrow head, 0: no head, negative: reverse head");
    readConfig(config, "color",                penColor,       Config::OPTIONAL, "",      "empty: from value");
    if(readConfigSequence(config, "scaleArrow", Config::OPTIONAL, "1", "draw an arrow for scale reference"))
    {
      readConfig(config, "originX", scaleArrowOriginX, Config::MUSTSET,  "0.055", "[0-1] 0: left, 1: right");
      readConfig(config, "originY", scaleArrowOriginY, Config::MUSTSET,  "0.065", "[0-1] 0: bottom, 1: top");
      readConfig(config, "length",  scaleArrowLength,  Config::MUSTSET,  "",      "in same unit as valueNorth and valueEast");
      readConfig(config, "unit",    scaleArrowUnit,    Config::MUSTSET,  "",      "displayed unit text (e.g. 1 cm)");
      readConfig(config, "label",   scaleArrowLabel,   Config::OPTIONAL, "",      "description of the arrows");
      endSequence(config);
      drawScaleArrow = TRUE;
    }
    if(isCreateSchema(config)) return;

    // tests
    hasZValues = !penColor && exprValue;
    if(!penColor && !exprValue)
      throw(Exception("value is needed to determine color of arrows"));
    if(!hasZValues)
      exprValue = nullptr;

    // read data
    GriddedData grid;
    readFileGriddedData(fileNameGrid, grid);
    points = grid.points;

    // create data variables
    // ---------------------
    VariableList varList;
    addDataVariables(grid, varList);

    std::vector<ExpressionVariablePtr> expressions = {exprValue, exprValueNorth, exprValueEast};
    for(ExpressionVariablePtr expr : expressions)
      if(expr) expr->simplify(varList);

    // evaluate expressions
    // --------------------
    data = Matrix(points.size(), exprValue ? 3 : 2);
    for(UInt i=0; i<points.size(); i++)
    {
      evaluateDataVariables(grid, i, varList);
      UInt idx = 0;
      for(auto expression : expressions)
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

std::string PlotMapLayerArrows::scriptEntry() const
{
  try
  {
    std::stringstream ss;

    // arrows from grid file
    ss<<"gmt psxy "<<dataFileName<<" -bi"<<2+data.columns()<<"d -J -R -A -SV"<<headSize<<"p+ea+z"<<scale<<"c";
    if(hasZValues)
      ss<<" -W"<<penSize<<"p+cl -CgroopsPlot.cpt";
    else
      ss<<" -W"<<penSize<<"p,"<<penColor->str()<<" -G"<<penColor->str();
    ss<<" -O -K >> groopsPlot.ps"<<std::endl;

    if(drawScaleArrow)
    {
      // scale arrow
      ss<<"echo "<<scaleArrowOriginX<<" "<<scaleArrowOriginY<<" 90 "<<scale*scaleArrowLength<<"c";
      ss<<" | gmt psxy -JX"<<PlotBasics::scriptVariable("width")<<"/"<<PlotBasics::scriptVariable("height")<<" -R0/1/0/1 -A -SV"<<headSize<<"p+ea+jc";
      if(penColor)
        ss<<" -W"<<penSize<<"p,"<<penColor->str()<<" -G"<<penColor->str();
      else
        ss<<" -W"<<penSize<<"p,black -Gblack";
      ss<<" -O -K >> groopsPlot.ps"<<std::endl;

      // scale arrow labels
      ss<<"echo "<<scaleArrowOriginX<<" "<<scaleArrowOriginY<<" CT "<<scaleArrowUnit<<" | gmt pstext -F+j -Dj4p -J -R -N -O -K >> groopsPlot.ps"<<std::endl;
      if(!scaleArrowLabel.empty())
        ss<<"echo "<<scaleArrowOriginX<<" "<<scaleArrowOriginY<<" CB "<<scaleArrowLabel<<" | gmt pstext -F+j -Dj4p -J -R -N -O -K >> groopsPlot.ps"<<std::endl;
      ss<<"gmt psbasemap -J"<<PlotBasics::scriptVariable("projection")<<" -R"<<PlotBasics::scriptVariable("region")<<" -B -O -K >> groopsPlot.ps"<<std::endl;
    }
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotMapLayerArrows::legendEntry(const FileName &workingDirectory, UInt idxLayer) const
{
  try
  {
    if(!drawScaleArrow)
      return "";

    FileName arrowDataFileName("arrow"+idxLayer%"%i.txt"s);
    FileName arrowLabelFileName("arrowLabel"+idxLayer%"%i.txt"s);

    std::stringstream ss;

    // scale arrow
    OutFile arrowFile(workingDirectory.append(arrowDataFileName));
    arrowFile<<scaleArrowOriginX<<" "<<scaleArrowOriginY<<" 90 "<<scale*scaleArrowLength<<"c";
    ss<<"gmt psxy "<<arrowDataFileName<<" -JX"<<PlotBasics::scriptVariable("width")<<"/"<<PlotBasics::scriptVariable("height")<<" -R0/1/0/1 -A -SV"<<headSize<<"p+ea+jc";

    if(penColor)
      ss<<" -W"<<penSize<<"p,"<<penColor->str()<<" -G"<<penColor->str();
    else
      ss<<" -W"<<penSize<<"p,black -Gblack";
    ss<<" -O -K >> groopsPlot.ps"<<std::endl;

    // scale arrow labels
    OutFile arrowLabelFile(workingDirectory.append(arrowLabelFileName));
    arrowLabelFile<<scaleArrowOriginX<<" "<<scaleArrowOriginY<<" CT "<<scaleArrowUnit;
    if(!scaleArrowLabel.empty())
      arrowLabelFile<<std::endl<<scaleArrowOriginX<<" "<<scaleArrowOriginY<<" CB "<<scaleArrowLabel;

    ss<<"gmt pstext "<<arrowLabelFileName<<" -F+j -Dj4p -J -R -N -O -K >> groopsPlot.ps"<<std::endl;

    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

// Latex documentation
static const char *docstringPlotMapLayerPolygon = R"(
\subsection{Polygon}
Draws a \configFile{inputfilePolygon}{polygon}.
If \configClass{fillColor}{plotColorType} is not set and a \config{value}
is given the fill color is taken from a \configClass{colorbar}{plotColorbarType}.
)";

class PlotMapLayerPolygon : public PlotMapLayer
{
  std::vector<Polygon> polygons;
  PlotLinePtr          line;
  PlotColorPtr         fillColor;
  Double               value;
  Bool                 greatCircle;

public:
  PlotMapLayerPolygon(Config &config);
  void boundary(const Ellipsoid &ellipsoid, Angle &minL, Angle &maxL, Angle &minB, Angle &maxB) const override;
  void writeDataFile(const Ellipsoid &ellipsoid, const FileName &workingDirectory, UInt idxLayer) override;
  Bool requiresColorBar() const override {return !fillColor && !std::isnan(value);}
  std::string scriptEntry() const override;
};

/***********************************************/

PlotMapLayerPolygon::PlotMapLayerPolygon(Config &config)
{
  try
  {
    FileName fileName;
    value = NAN_EXPR;

    readConfig(config, "inputfilePolygon", fileName,  Config::MUSTSET,  "{groopsDataDir}/border/", "");
    readConfig(config, "line",             line,      Config::OPTIONAL, "solid", "style of border lines");
    readConfig(config, "fillColor",        fillColor, Config::OPTIONAL, "", "polygon fill color (no fill color: determine from value if given, else: no fill)");
    readConfig(config, "value",            value,     Config::OPTIONAL, "", "value to compute fill color from a colorbar (ignored if a fillColor is given)");
    readConfig(config, "drawLineAsGreatCircle", greatCircle, Config::DEFAULT, "1", "draw connecting lines as great circles (otherwise, a straight line is drawn instead)");
    if(isCreateSchema(config)) return;

    // tests
    if(!line && !(fillColor || !std::isnan(value)) )
      throw(Exception("At least one of line, fillColor or value must be set."));

    readFilePolygon(fileName, polygons);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotMapLayerPolygon::boundary(const Ellipsoid &/*ellipsoid*/, Angle &minL, Angle &maxL, Angle &minB, Angle &maxB) const
{
  try
  {
    for(auto &p : polygons)
      for(UInt i=0; i<p.L.size(); i++)
      {
        minL = std::min(minL, Angle(p.L(i)));
        maxL = std::max(maxL, Angle(p.L(i)));
        minB = std::min(minB, Angle(p.B(i)));
        maxB = std::max(maxB, Angle(p.B(i)));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotMapLayerPolygon::writeDataFile(const Ellipsoid &/*ellipsoid*/, const FileName &workingDirectory, UInt idxLayer)
{
  try
  {
    dataFileName = "polygon."+idxLayer%"%i.dat"s;
    OutFile file(workingDirectory.append(dataFileName));
    if(!std::isnan(value)) file<<">-Z"<<value<<std::endl;
    for(auto &p : polygons)
    {
      for(UInt i=0; i<p.L.size(); i++)
      {
        file<<p.L(i)*RAD2DEG<<" "<<p.B(i)*RAD2DEG<<std::endl;
      }
      file<<">"<<std::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotMapLayerPolygon::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt psxy "<<dataFileName<<" -L -J -R";
    if(fillColor) ss<<" -G"<<fillColor->str();
    else if(!std::isnan(value)) ss<<" -CgroopsPlot.cpt";
    if(line)
    {
      ss<<" -W"<<line->str();
      if(!greatCircle) ss<<" -A";
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
/***********************************************/

// Latex documentation
static const char *docstringPlotMapLayerCoast = R"(
\subsection{Coast}
Plots coastlines. GMT provides them in different \config{resolution}s.
Features with an area smaller than \config{minArea} in $km^2$ will not be plotted.
)";

class PlotMapLayerCoast : public PlotMapLayer
{
  std::string  resolution;
  PlotLinePtr  line;
  PlotColorPtr landColor;
  PlotColorPtr oceanColor;
  UInt         minArea;

public:
  PlotMapLayerCoast(Config &config);
  std::string scriptEntry() const override;
};

/***********************************************/

PlotMapLayerCoast::PlotMapLayerCoast(Config &config)
{
  try
  {
    std::string choice;
    if(readConfigChoice(config, "resolution", choice, Config::MUSTSET, "medium", ""))
    {
      if(readConfigChoiceElement(config, "crude",  choice)) resolution = "c";
      if(readConfigChoiceElement(config, "low",    choice)) resolution = "l";
      if(readConfigChoiceElement(config, "medium", choice)) resolution = "i";
      if(readConfigChoiceElement(config, "high",   choice)) resolution = "h";
      if(readConfigChoiceElement(config, "full",   choice)) resolution = "f";
      endChoice(config);
    }
    readConfig(config, "line",       line,       Config::OPTIONAL, R"({"solid": {"width": "1"}})", "line style for coastlines");
    readConfig(config, "landColor",  landColor,  Config::OPTIONAL, "",      "fill land area");
    readConfig(config, "oceanColor", oceanColor, Config::OPTIONAL, "",      "fill ocean area");
    readConfig(config, "minArea",    minArea,    Config::DEFAULT,  "5000", "[km^2] features with a smaller area than this are dropped");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotMapLayerCoast::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt pscoast -J -R -D"<<resolution<<"+ -A"<<minArea;
    if(line)
      ss<<" -W"<<line->str();
    if(landColor)
      ss<<" -G"<<landColor->str();
    if(oceanColor)
      ss<<" -S"<<oceanColor->str();
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

// Latex documentation
static const char *docstringPlotMapLayerRivers = R"(
\subsection{Rivers}
Plots rivers and lakes. GMT provides different classes
(\url{https://docs.generic-mapping-tools.org/latest/coast.html}).
)";

class PlotMapLayerRivers : public PlotMapLayer
{
  std::vector<std::string> riverclass;
  PlotLinePtr line;

public:
  PlotMapLayerRivers(Config &config);
  std::string scriptEntry() const override;
};

/***********************************************/

PlotMapLayerRivers::PlotMapLayerRivers(Config &config)
{
  try
  {
    std::string choice;
    if(readConfigChoice(config, "class", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "riversCanalsLakes",     choice, "")) riverclass.push_back("a");
      if(readConfigChoiceElement(config, "riversCanals",          choice, "")) riverclass.push_back("A");
      if(readConfigChoiceElement(config, "permanentRiversLakes",  choice, "")) riverclass.push_back("r");
      if(readConfigChoiceElement(config, "permanentRivers",       choice, "")) riverclass.push_back("R");
      if(readConfigChoiceElement(config, "intermittentRivers",    choice, "")) riverclass.push_back("i");
      if(readConfigChoiceElement(config, "canals",                choice, "")) riverclass.push_back("c");
      if(readConfigChoiceElement(config, "singleClass",           choice, ""))
      {
        std::vector<UInt> types;
        readConfig(config, "class", types, Config::MUSTSET, "", "0-10. See GMT documentation");

        for(UInt k=0; k<types.size(); k++)
          riverclass.push_back(types.at(k)%"%i"s);
      }
      endChoice(config);
    }
    readConfig(config, "line", line, Config::MUSTSET, R"({"solid": {"width": "1"}})", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotMapLayerRivers::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt pscoast -J -R";
    for(UInt k=0; k<riverclass.size(); k++)
      ss<<" -I"<<riverclass.at(k)<<"/"<<line->str();
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

// Latex documentation
static const char *docstringPlotMapLayerPolitical = R"(
\subsection{PoliticalBoundary}
Plots national boundaries. GMT provides them in different \config{resolution}s.
)";

class PlotMapLayerPolitical : public PlotMapLayer
{
  std::string resolution;
  PlotLinePtr line;

public:
  PlotMapLayerPolitical(Config &config);
  std::string scriptEntry() const override;
};

/***********************************************/

PlotMapLayerPolitical::PlotMapLayerPolitical(Config &config)
{
  try
  {
    std::string choice;
    if(readConfigChoice(config, "resolution", choice, Config::MUSTSET, "medium", ""))
    {
      if(readConfigChoiceElement(config, "crude",  choice, "")) resolution = "c";
      if(readConfigChoiceElement(config, "low",    choice, "")) resolution = "l";
      if(readConfigChoiceElement(config, "medium", choice, "")) resolution = "i";
      if(readConfigChoiceElement(config, "high",   choice, "")) resolution = "h";
      if(readConfigChoiceElement(config, "full",   choice, "")) resolution = "f";
      endChoice(config);
    }
    readConfig(config, "line", line, Config::MUSTSET, R"({"solid": {"width": "0.25"}})", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotMapLayerPolitical::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt pscoast -J -R -N1/"<<line->str()<<" -D"<<resolution<<"+ -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

// Latex documentation
static const char *docstringPlotMapLayerBlueMarble = R"(
\subsection{BlueMarble}
An image of the Earth's surface as seen from outer space -
the image is known as \emph{blue marble}. The directory of \config{inputfileChannels}
contains several files in different resolutions representing the Earth's surface each
month throughout a year.

\fig{!hb}{0.8}{blueMarble}{fig:blueMarbleMap}{The blue marble.}
)";

class PlotMapLayerBlueMarble : public PlotMapLayer
{
  FileName fileNameImage;
  Double   brightness;

  Bool     illuminate;
  FileName fileNameTopography;
  Angle    azimuth, elevation;
  Double   amplitude;
  Double   ambient, diffuse, specular, shine;

public:
  PlotMapLayerBlueMarble(Config &config);
  void        writeDataFile(const Ellipsoid &ellipsoid, const FileName &workingDirectory, UInt idxLayer) override;
  std::string scriptEntry() const override;
};

/***********************************************/

PlotMapLayerBlueMarble::PlotMapLayerBlueMarble(Config &config)
{
  try
  {
    readConfig(config, "inputfileImage",    fileNameImage,    Config::MUSTSET, "{groopsDataDir}/plot/bluemarble/1800x900/bluemarble.05.jpg", "Blue Marble image file");
    readConfig(config, "brightness",        brightness,       Config::DEFAULT, "0.0", "brightness of bitmap [-1, 1]");
    if(readConfigSequence(config, "illuminate", Config::OPTIONAL, "", "add hillshade based on topography"))
    {
      readConfig(config, "inputfileTopography", fileNameTopography, Config::MUSTSET, "{groopsDataDir}/plot/bluemarble/1800x900/bluemarble.topography.grd", "GMT grid file containing topography.");
      readConfig(config, "azimuth",             azimuth,            Config::DEFAULT, "315", "direction of lighting source [deg]");
      readConfig(config, "elevation",           elevation,          Config::DEFAULT, "45",  "direction of lighting source [deg]");
      readConfig(config, "ambient",             ambient,            Config::DEFAULT, "0.1", "ambient lighting");
      readConfig(config, "diffuse",             diffuse,            Config::DEFAULT, "0.6", "diffuse lighting");
      readConfig(config, "specular",            specular,           Config::DEFAULT, "0.0", "specular reflection");
      readConfig(config, "shine",               shine,              Config::DEFAULT, "0.0", "surface shine");
      readConfig(config, "amplitude",           amplitude,          Config::DEFAULT, "0.1", "scale gradient by factor");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    illuminate = !fileNameTopography.empty();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotMapLayerBlueMarble::writeDataFile(const Ellipsoid &/*ellipsoid*/, const FileName &/*workingDirectory*/, UInt idxLayer)
{
  dataFileName = "blueMarble."+idxLayer%"%i.dat"s;
}

/***********************************************/

std::string PlotMapLayerBlueMarble::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    if(illuminate)
    {
      ss<<"gmt grdgradient "<<fileNameTopography<<" -G"<<dataFileName;
      ss<<" -Em"<<azimuth*RAD2DEG<<"/"<<elevation*RAD2DEG<<"+a"<<ambient<<"+d"<<diffuse<<"+p"<<specular<<"+s"<<shine<<std::endl;
      ss<<"gmt grdmath "<<dataFileName<<" "<<amplitude<<" MUL = "<<dataFileName<<std::endl;
    }

    ss<<"gmt grdimage "<<fileNameImage;
    if(!illuminate)
      ss<<" -I"<<brightness;
    else
      ss<<" -I"<<dataFileName;
    ss<<" -Q -J -R -O -K >> groopsPlot.ps"<<std::endl;

    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

// Latex documentation
static const char *docstringPlotMapLayerText = R"(
\subsection{Text}
Writes a \config{text} at \config{originLongitude} and \config{originLatitude} position in the map.
With \config{clip} the text is cutted at the boundaries of the plotting area.
)";

class PlotMapLayerText : public PlotMapLayer
{
  Angle        lat, lon;
  Double       xOffset, yOffset;
  Double       fontSize;
  PlotColorPtr fontColor;
  std::string  text;
  std::string  alignment;
  Bool         clip;

public:
  PlotMapLayerText(Config &config);
  void        writeDataFile(const Ellipsoid &ellipsoid, const FileName &workingDirectory, UInt idxLayer) override;
  std::string scriptEntry() const override;
};

/***********************************************/

PlotMapLayerText::PlotMapLayerText(Config &config)
{
  try
  {
    readConfig(config, "text",            text,      Config::MUSTSET,  "",   "");
    readConfig(config, "originLongitude", lon,       Config::MUSTSET,  "",   "[deg]");
    readConfig(config, "originLatitude",  lat,       Config::MUSTSET,  "",   "[deg]");
    readConfig(config, "offsetX",         xOffset,   Config::DEFAULT,  "0",  "[cm] x-offset from origin");
    readConfig(config, "offsetY",         yOffset,   Config::DEFAULT,  "0",  "[cm] y-offset from origin");
    readConfig(config, "alignment",       alignment, Config::DEFAULT,  "LB", "L, C, R (left, center, right) and T, M, B (top, middle, bottom)");
    readConfig(config, "fontSize",        fontSize,  Config::DEFAULT,  "10", "");
    readConfig(config, "fontColor",       fontColor, Config::MUSTSET,  "",   "");
    readConfig(config, "clip",            clip,      Config::DEFAULT,  "1",  "clip at boundaries");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotMapLayerText::writeDataFile(const Ellipsoid &/*ellipsoid*/, const FileName &workingDirectory, UInt idxLayer)
{
  dataFileName = "text."+idxLayer%"%i.txt"s;
  OutFile file(workingDirectory.append(dataFileName));

  file<<lon*RAD2DEG<<" "<<lat*RAD2DEG<<" "<<text;
}

/***********************************************/

std::string PlotMapLayerText::scriptEntry() const
{
  try
  {
    std::stringstream ss;
    ss<<"gmt pstext "<<dataFileName<<" -F+f"<<fontSize<<"p,,"<<fontColor->str()<<"+j"<<alignment<<" -D"<<xOffset<<"/"<<yOffset<<" -J -R "<<(clip?"":"-N")<<" -O -K >> groopsPlot.ps"<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GROOPS_REGISTER_CLASS(PlotMapLayer, "plotMapLayerType",
                      PlotMapLayerGrid,
                      PlotMapLayerPoints,
                      PlotMapLayerArrows,
                      PlotMapLayerPolygon,
                      PlotMapLayerCoast,
                      PlotMapLayerRivers,
                      PlotMapLayerPolitical,
                      PlotMapLayerBlueMarble,
                      PlotMapLayerText)

GROOPS_READCONFIG_CLASS(PlotMapLayer, "plotMapLayerType")

/***********************************************/

PlotMapLayerPtr PlotMapLayer::create(Config &config, const std::string &name)
{
  try
  {
    PlotMapLayerPtr plotMapLayer;
    std::string     type;

    readConfigChoice(config, name, type, Config::MUSTSET, "", "plot layers");
    if(readConfigChoiceElement(config, "griddedData",       type, "data with regular sampling"))
      plotMapLayer = PlotMapLayerPtr(new PlotMapLayerGrid(config));
    if(readConfigChoiceElement(config, "points",            type, "colored points"))
      plotMapLayer = PlotMapLayerPtr(new PlotMapLayerPoints(config));
    if(readConfigChoiceElement(config, "arrows",            type, "colored arrows"))
      plotMapLayer = PlotMapLayerPtr(new PlotMapLayerArrows(config));
    if(readConfigChoiceElement(config, "polygon",           type, "boundaries"))
      plotMapLayer = PlotMapLayerPtr(new PlotMapLayerPolygon(config));
    if(readConfigChoiceElement(config, "coast",             type, "coast lines, ..."))
      plotMapLayer = PlotMapLayerPtr(new PlotMapLayerCoast(config));
    if(readConfigChoiceElement(config, "rivers",            type, "major rivers"))
      plotMapLayer = PlotMapLayerPtr(new PlotMapLayerRivers(config));
    if(readConfigChoiceElement(config, "politicalBoundary", type, "political boundaries"))
      plotMapLayer = PlotMapLayerPtr(new PlotMapLayerPolitical(config));
    if(readConfigChoiceElement(config, "blueMarble",        type, "blue marble"))
      plotMapLayer = PlotMapLayerPtr(new PlotMapLayerBlueMarble(config));
    if(readConfigChoiceElement(config, "text",              type, "text"))
      plotMapLayer = PlotMapLayerPtr(new PlotMapLayerText(config));
    endChoice(config);

    return plotMapLayer;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotMapLayer::boundary(const Ellipsoid &ellipsoid, Angle &minL, Angle &maxL, Angle &minB, Angle &maxB) const
{
  try
  {
    for(UInt k=0; k<points.size(); k++)
    {
      Angle  lon, lat;
      Double h;
      ellipsoid(points.at(k), lon, lat, h);
      minL = std::min(minL, lon-bufferLon);
      maxL = std::max(maxL, lon+bufferLon);
      minB = std::min(minB, lat-bufferLat);
      maxB = std::max(maxB, lat+bufferLat);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotMapLayer::getIntervalZ(Bool isLogarithmic, Double &minZ, Double &maxZ) const
{
  try
  {
    if(!requiresColorBar())
      return;

    UInt   count =  0;
    Double avg   =  0.;
    for(UInt i=0; i<data.rows(); i++)
      if(!std::isnan(data(i, 0)) && (!isLogarithmic || (data(i, 0) > 0)))
      {
        minZ  = std::min(minZ, data(i, 0));
        maxZ  = std::max(maxZ, data(i, 0));
        avg  += std::fabs(data(i, 0));
        count++;
      }
    avg /= count;

    if(!isLogarithmic)
    {
      minZ = (minZ >= 0) ? 0 : -3*avg;
      maxZ = +3*avg;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlotMapLayer::writeDataFile(const Ellipsoid &ellipsoid, const FileName &workingDirectory, UInt idxLayer)
{
  try
  {
    if(!points.size())
      return;

    dataFileName = "data."+idxLayer%"%i.dat"s;
    OutFile file(workingDirectory.append(dataFileName), std::ios::out | std::ios::binary);
    for(UInt i=0; i<points.size(); i++)
    {
      Angle  lon, lat;
      Double h;
      ellipsoid(points.at(i), lon, lat, h);
      std::vector<Double> line = {lon*RAD2DEG, lat*RAD2DEG};
      file.write(reinterpret_cast<char*>(line.data()), line.size()*sizeof(Double));
      for(UInt k=0; k<data.columns(); k++)
        file.write(reinterpret_cast<char*>(&data(i, k)), sizeof(Double));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotMapLayer::scriptStatisticsInfo(UInt fontSize, Double width, const FileName &workingDirectory, UInt idxLayer) const
{
  try
  {
    if(!requiresColorBar())
      return std::string();

    Double rms, avg, vmin, vmax, mean;
    MiscGriddedData::statistics(Vector(data.column(0)), areas, rms, avg, vmin, vmax, mean);

    FileName fileNameStatistics("statistics"+idxLayer%"%i.txt"s);

    OutFile statisticsFile(workingDirectory.append(fileNameStatistics));
    statisticsFile<<"0.5 0.0 "<<fontSize<<"p CB min="<<vmin<<", max="<<vmax<<", mean="<<mean<<", rms="<<rms;

    std::stringstream ss;
    ss<<"gmt pstext "<<fileNameStatistics<<" -F+f+j -Y-"<<fontSize<<"p -JX"<<width<<"c/"<<fontSize<<"p -R0/1/0/1 -N -O -K >> groopsPlot.ps "<<PlotBasics::scriptError2Null()<<std::endl;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string PlotMapLayer::legendEntry(const FileName &/*workingDirectory*/, UInt /*idxLayer*/) const
{
  try
  {
    return "";
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
