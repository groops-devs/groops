/***********************************************/
/**
* @file plotSphericalHarmonicsTriangle.cpp
*
* @brief Plot triangle of potential coefficients.
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
Plot the potential coefficients of a spherical harmonic expansion
using the GMT Generic Mapping Tools (\url{https://www.generic-mapping-tools.org}).
A variety of image file formats are supported (e.g. png, jpg, eps) determined by the extension of \config{outputfile}.

This program plots the formal errors (sigmas).
If \configClass{gravityfield}{gravityfieldType} provides no sigmas
e.g. with \config{setSigmasToZero} in \configClass{gravityfield:potentialCoefficients}{gravityfieldType:potentialCoefficients}
the coefficients itself are plotted instead.

The plot programs create a temporary directory in the path of \config{outputfile}, writes all needed data into it,
generates a batch/shell script with the GMT commands, execute it, and remove the temporary directory.
With setting \config{options:removeFiles}=false the last step is skipped and it is possible to adjust the plot manually
to specific publication needs. Individual GMT settings are adjusted with \config{options:options}="\verb|FORMAT=value|",
see \url{https://docs.generic-mapping-tools.org/latest/gmt.conf.html}.

\fig{!hb}{0.8}{plotSphericalHarmonicsTriangle}{fig:plotSphericalHarmonicsTriangle}{Formal errors of GOCO06s.}
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "classes/gravityfield/gravityfield.h"
#include "plot/plotColorbar.h"
#include "plot/plotMisc.h"

/***********************************************/

/** @brief Plot triangle of potential coefficients.
* @ingroup programsGroup */
class PlotSphericalHarmonicsTriangle
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(PlotSphericalHarmonicsTriangle, SINGLEPROCESS, "Plot triangle of potential coefficients", Plot, PotentialCoefficients)

/***********************************************/

void PlotSphericalHarmonicsTriangle::run(Config &config)
{
  try
  {
    FileName        fileNamePlot;
    std::string     title;
    GravityfieldPtr gravityfield;
    Double          annotation=NAN_EXPR, frame=NAN_EXPR, grid=NAN_EXPR;
    Time            time;
    UInt            minDegree=0, maxDegree=INFINITYDEGREE;
    PlotLinePtr     gridLine;
    PlotColorbarPtr colorbar;
    PlotBasics      plotBasics;

    renameDeprecatedConfig(config, "annotation", "majorTickSpacing", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "frame",      "minorTickSpacing", date2time(2020, 4, 23));
    renameDeprecatedConfig(config, "grid",       "gridLineSpacing",  date2time(2020, 4, 23));

    readConfig(config, "outputfile",       fileNamePlot,  Config::MUSTSET,  "",  "*.png, *.jpg, *.eps, ...");
    readConfig(config, "title",            title,         Config::OPTIONAL, "",  "");
    readConfig(config, "gravityfield",     gravityfield,  Config::MUSTSET,  "",  "use sigmas, if not given use signal (cnm,snm)");
    readConfig(config, "time",             time,          Config::OPTIONAL, "",  "at this time the gravity field will be evaluated");
    readConfig(config, "minDegree",        minDegree,     Config::OPTIONAL, "",  "");
    readConfig(config, "maxDegree",        maxDegree,     Config::OPTIONAL, "",  "");
    readConfig(config, "majorTickSpacing", annotation,    Config::OPTIONAL, "",  "boundary annotation");
    readConfig(config, "minorTickSpacing", frame,         Config::OPTIONAL, "",  "frame tick spacing");
    readConfig(config, "gridLineSpacing",  grid,          Config::DEFAULT,  "0", "gridline spacing");
    readConfig(config, "gridLine",         gridLine,      Config::OPTIONAL, R"({"solid": {"width":"0.25", "color":"gray"}})", "The style of the grid lines.");
    readConfig(config, "colorbar",         colorbar,      Config::MUSTSET,  R"({"logarithmic":"1"})", "");
    plotBasics.read(config, "PlotSphericalHarmonicsTriangle", fileNamePlot, title, "12", "");
    if(isCreateSchema(config)) return;

    // ===========================================================

    logStatus<<"use accuracies, if not given use signal"<<Log::endl;
    SphericalHarmonics harm = gravityfield->sphericalHarmonics(time, maxDegree, minDegree);
    maxDegree = harm.maxDegree();
    Matrix cnm = harm.sigma2cnm();
    Matrix snm = harm.sigma2snm();
    if(cnm.size() == 0)
      cnm = snm = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    for(UInt n=minDegree; n<=maxDegree; n++)
      for(UInt m=0; m<=n; m++)
      {
        cnm(n,m) = cnm(n,m) ? std::sqrt(cnm(n,m)) : harm.cnm()(n,m);
        snm(n,m) = snm(n,m) ? std::sqrt(snm(n,m)) : harm.snm()(n,m);
      }

    // calculate min and max
    // ---------------------
    Double minZ =  1e99;
    Double maxZ = -1e99;
    for(UInt n=minDegree; n<=maxDegree; n++)
      for(UInt m=0; m<=n; m++)
      {
        if(colorbar->isLogarithmic())
        {
          cnm(n,m) = std::fabs(cnm(n,m));
          snm(n,m) = std::fabs(snm(n,m));
        }
        if(cnm(n,m))
        {
          minZ = std::min(minZ, cnm(n,m));
          maxZ = std::max(maxZ, cnm(n,m));
        }
        if(snm(n,m))
        {
          minZ = std::min(minZ, snm(n,m));
          maxZ = std::max(maxZ, snm(n,m));
        }
      }
    colorbar->setAutoInterval(minZ, maxZ);

    if(std::isnan(plotBasics.height))
      plotBasics.height = plotBasics.width*(maxDegree+1-minDegree)/(2*maxDegree+1);
    const Double shiftX = plotBasics.width*maxDegree/(2*maxDegree+1);

    // write data files
    // ----------------
    logStatus<<"create temporary data files"<<Log::endl;
    {
      OutFile file1(plotBasics.workingDirectory.append("data.cnm.dat"), std::ios::out | std::ios::binary);
      OutFile file2(plotBasics.workingDirectory.append("data.snm.dat"), std::ios::out | std::ios::binary);
      for(UInt n=minDegree; n<=maxDegree; n++)
        for(UInt m=0; m<=n; m++)
        {
          if(cnm(n,m))
          {
            const std::vector<Double> line = {Double(n), Double(m), cnm(n,m)};
            file1.write(reinterpret_cast<const char*>(line.data()), line.size()*sizeof(Double));
          }
          if(snm(n,m))
          {
            const std::vector<Double> line = {Double(n), Double(m), snm(n,m)};
            file2.write(reinterpret_cast<const char*>(line.data()), line.size()*sizeof(Double));
          }
        }
    }

    // create scriptfile
    // -----------------
    logStatus<<"create scriptfile"<<Log::endl;
    {
      OutFile file(plotBasics.fileNameScript());
      file<<plotBasics.scriptHeader();
      if(colorbar)
        file<<colorbar->scriptColorTable();
      file<<std::endl;
      file<<"gmt psbasemap -Y"<<-plotBasics.height-plotBasics.marginTitle<<"c -R0.5/"<<maxDegree+0.5<<"/-0.5/"<<maxDegree+0.5<<" -JX"<<-shiftX<<"c/"<<-plotBasics.height<<"c -BWSn";
      file<<" -Bx"<<PlotBasics::axisTicks(FALSE, 0.5, maxDegree+0.5, annotation, frame, 0, "", "s@-nm@-");
      file<<" -By"<<PlotBasics::axisTicks(FALSE, 0.5, maxDegree+0.5, annotation, frame, 0, "", "degree");
      file<<" -O -K >> groopsPlot.ps"<<std::endl;
      file<<"gmt xyz2grd   data.snm.dat -Gdata.snm.grd -: -r -bi3d -I1 -R"<<std::endl;
      file<<"gmt grdimage  data.snm.grd --COLOR_NAN=255/255/255 -Q -R -J -CgroopsPlot.cpt -B -O -K >> groopsPlot.ps"<<std::endl;
      if(grid > 0)
        file<<"gmt psbasemap -R -J -Bwsn -Bg"<<grid<<"-0.5 --MAP_GRID_PEN_PRIMARY="<<gridLine->str()<<" -O -K >> groopsPlot.ps"<<std::endl;
      file<<std::endl;
      file<<"gmt psbasemap -X"<<shiftX<<"c -Y0c -R-0.5/"<<maxDegree+0.5<<"/-0.5/"<<maxDegree+0.5<<" -JX"<<plotBasics.width*(maxDegree+1)/(2*maxDegree+1)<<"c/-"<<plotBasics.height<<"c -BSEn";
      file<<" -Bx"<<PlotBasics::axisTicks(FALSE, 0.5, maxDegree+0.5, annotation, frame, 0, "", "c@-nm@-");
      file<<" -By"<<PlotBasics::axisTicks(FALSE, 0.5, maxDegree+0.5, annotation, frame, 0, "", "");
      file<<" -O -K >> groopsPlot.ps"<<std::endl;
      file<<"gmt xyz2grd   data.cnm.dat -Gdata.cnm.grd -: -r -bi3d  -I1 -R"<<std::endl;
      file<<"gmt grdimage  data.cnm.grd --COLOR_NAN=255/255/255 -Q -R -J -CgroopsPlot.cpt -B -O -K >> groopsPlot.ps"<<std::endl;
      if(grid > 0)
        file<<"gmt psbasemap -R -J -Besn -Bg"<<grid<<"-0.5 --MAP_GRID_PEN_PRIMARY="<<gridLine->str()<<" -O -K >> groopsPlot.ps"<<std::endl;
      file<<"gmt psxy -X"<<-shiftX<<"c -R -J -T -O -K >> groopsPlot.ps"<<std::endl;
      file<<std::endl;
      if(colorbar)
        file<<colorbar->scriptEntry(plotBasics.width, plotBasics.height, 0.8/*marginX*/, 0.8/*marginY*/);
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
