/***********************************************/
/**
* @file plotMatrix.cpp
*
* @brief Plot coefficients of a Matrix.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2008-08-29
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Plot the coefficients of a \configFile{inputfileMatrix}{matrix}
using the GMT Generic Mapping Tools (\url{https://www.generic-mapping-tools.org}).
A variety of image file formats are supported (e.g. png, jpg, eps) determined by the extension of \config{outputfile}.

The plot programs create a temporary directory in the path of \config{outputfile}, writes all needed data into it,
generates a batch/shell script with the GMT commands, execute it, and remove the temporary directory.
With setting \config{options:removeFiles}=false the last step is skipped and it is possible to adjust the plot manually
to specific publication needs. Individual GMT settings are adjusted with \config{options:options}="\verb|FORMAT=value|",
see \url{https://docs.generic-mapping-tools.org/latest/gmt.conf.html}.

\fig{!hb}{0.6}{plotMatrix}{fig:plotMatrix}{Upper left part of the DDK filter matrix.}
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileMatrix.h"
#include "plot/plotColorbar.h"
#include "plot/plotMisc.h"

/***********************************************/

/** @brief Plot coefficients of a Matrix.
* @ingroup programsGroup */
class PlotMatrix
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(PlotMatrix, SINGLEPROCESS, "Plot coefficients of a Matrix.", Plot, Matrix)

/***********************************************/

void PlotMatrix::run(Config &config)
{
  try
  {
    FileName        fileNamePlot;
    FileName        fileNameMatrix;
    std::string     title;
    UInt            minY, maxY=NULLINDEX, minX, maxX=NULLINDEX;
    Double          annotationX=NAN_EXPR, frameX=NAN_EXPR, gridX=NAN_EXPR;
    Double          annotationY=NAN_EXPR, frameY=NAN_EXPR, gridY=NAN_EXPR;
    PlotLinePtr     gridLine;
    PlotColorbarPtr colorbar;
    PlotBasics      plotBasics;

    renameDeprecatedConfig(config, "annotationX", "majorTickSpacingX", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "frameX",      "minorTickSpacingX", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "gridX",       "gridLineSpacingX",  date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "annotationY", "majorTickSpacingY", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "frameY",      "minorTickSpacingY", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "gridY",       "gridLineSpacingY",  date2time(2020, 4, 23));

    readConfig(config, "outputfile",        fileNamePlot,   Config::MUSTSET,  "", "*.png, *.jpg, *.eps, ...");
    readConfig(config, "title",             title,          Config::OPTIONAL, "",  "");
    readConfig(config, "inputfileMatrix",   fileNameMatrix, Config::MUSTSET,  "",  "");
    readConfig(config, "minColumn",         minX,           Config::DEFAULT,  "0", "minimum column index to plot");
    readConfig(config, "maxColumn",         maxX,           Config::OPTIONAL, "",  "maximum column index to plot");
    readConfig(config, "majorTickSpacingX", annotationX,    Config::OPTIONAL, "",  "boundary annotation");
    readConfig(config, "minorTickSpacingX", frameX,         Config::OPTIONAL, "",  "frame tick spacing");
    readConfig(config, "gridLineSpacingX",  gridX,          Config::OPTIONAL, "",  "gridline spacing");
    readConfig(config, "minRow",            minY,           Config::DEFAULT,  "0", "minimum row index to plot");
    readConfig(config, "maxRow",            maxY,           Config::OPTIONAL, "",  "maximum row index to plot");
    readConfig(config, "majorTickSpacingY", annotationY,    Config::OPTIONAL, "",  "boundary annotation");
    readConfig(config, "minorTickSpacingY", frameY,         Config::OPTIONAL, "",  "frame tick spacing");
    readConfig(config, "gridLineSpacingY",  gridY,          Config::OPTIONAL, "",  "gridline spacing");
    readConfig(config, "gridLine",          gridLine,       Config::OPTIONAL, R"({"solid": {"width":"0.25", "color":"gray"}})", "The style of the grid lines.");
    readConfig(config, "colorbar",          colorbar,       Config::MUSTSET,  "",  "");
    plotBasics.read(config, "PlotMatrix", fileNamePlot, title, "12", "");
    if(isCreateSchema(config)) return;

    // ===========================================================

    // read matrix
    // -----------
    logStatus<<"read matrix"<<Log::endl;
    Matrix M;
    readFileMatrix(fileNameMatrix, M);
    if(M.getType()==Matrix::SYMMETRIC)
      fillSymmetric(M);

    // calculate min and max
    // ---------------------
    if(maxX==NULLINDEX) maxX = M.columns()-1;
    if(maxY==NULLINDEX) maxY = M.rows()-1;
    Double minZ =  1e99;
    Double maxZ = -1e99;
    for(UInt z=0; z<M.rows(); z++)
      for(UInt s=0; s<M.columns(); s++)
      {
        if(colorbar->isLogarithmic())
          M(z,s) = std::fabs(M(z,s));
        if(M(z,s) || !colorbar->isLogarithmic())
        {
          minZ = std::min(minZ, M(z,s));
          maxZ = std::max(maxZ, M(z,s));
        }
      }
    colorbar->setAutoInterval(minZ, maxZ);

    logInfo<<"  matrix dimension: ("<<M.rows()<<" x "<<M.columns()<<")"<<Log::endl;
    logInfo<<"  plot z range  ("<<colorbar->getMin()<<" .. "<<colorbar->getMax()<<") of ("<<minZ<<" .. "<<maxZ<<")"<<Log::endl;

    // determine annotation x
    const Double scaleX = std::pow(10.0, 1-std::ceil(std::log10(maxX-minX))); // only one significant digit
    const Double scaleY = std::pow(10.0, 1-std::ceil(std::log10(maxY-minY))); // only one significant digit
    if(std::isnan(gridX)) gridX = std::max(1., std::round((maxX-minX)*scaleX)/(scaleX*10));
    if(std::isnan(gridY)) gridY = std::max(1., std::round((maxY-minY)*scaleY)/(scaleY*10));

    if(std::isnan(plotBasics.height))
      plotBasics.height = plotBasics.width*(maxY-minY)/(maxX-minX);

    // write data files
    // ----------------
    logStatus<<"create temporary data files"<<Log::endl;
    {
      OutFile file(plotBasics.workingDirectory.append("data.dat"), std::ios::out | std::ios::binary);
      for(UInt z=minY; z<=std::min(maxY, M.rows()-1); z++)
        for(UInt s=minX; s<=std::min(maxX, M.columns()-1); s++)
          if(M(z,s))
          {
            std::vector<Double> line = {Double(s), Double(z), M(z,s)};
            file.write(reinterpret_cast<char*>(line.data()), line.size()*sizeof(Double));
          }
    }

    // create scriptfile
    // -----------------
    logStatus<<"create scriptfile"<<Log::endl;
    {
      OutFile file(plotBasics.fileNameScript());
      file<<plotBasics.scriptHeader();
      file<<colorbar->scriptColorTable();
      file<<"gmt psbasemap -Y"<<-plotBasics.height-plotBasics.marginTitle<<"c -R"<<minX-0.5<<"/"<<maxX+0.5<<"/"<<minY-0.5<<"/"<<maxY+0.5<<" -JX"<<plotBasics.width<<"c/"<<-plotBasics.height<<"c";
      file<<" -BWSne -Bx"<<PlotBasics::axisTicks(FALSE, minX, maxX, annotationX, frameX, 0)<<" -By"<<PlotBasics::axisTicks(FALSE, minY, maxY, annotationY, frameY, 0);
      file<<" -O -K >> groopsPlot.ps"<<std::endl;
      file<<"gmt xyz2grd data.dat -bi3d -Gdata.grd -r -V -I1 -R"<<std::endl;
      file<<"gmt grdimage  data.grd --COLOR_NAN=255/255/255 -Q -R -J -CgroopsPlot.cpt";
      if(gridLine)
        file<<" -Bxg"<<gridX<<"-0.5 -Byg"<<gridY<<"-0.5 --MAP_GRID_PEN_PRIMARY="<<gridLine->str();
      file<<" -O -K >> groopsPlot.ps"<<std::endl;
      file<<colorbar->scriptEntry(plotBasics.width, plotBasics.height, 0.0/*marginX*/, 0.8/*marginY*/);
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
