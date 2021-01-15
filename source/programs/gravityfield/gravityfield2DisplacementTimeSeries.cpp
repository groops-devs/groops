/***********************************************/
/**
* @file gravityfield2DisplacementTimeSeries.cpp
*
* @brief Generates a time series of station displacements.
*
* @author Torsten Mayer-Guerr
* @date 2010-11-01
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes a time series of displacements of a list of stations (\configClass{grid}{gridType})
due to the effect of time variable loading masses. The displacement~$\M u$ of a station is calculated according to
\begin{equation}\label{eq:displacement}
\M u(\M r) = \frac{1}{\gamma}\sum_{n=0}^\infty \left[\frac{h_n}{1+k_n}V_n(\M r)\,\M e_{up}
+ R\frac{l_n}{1+k_n}\left(
 \frac{\partial V_n(\M r)}{\partial \M e_{north}}\M e_{north}
+\frac{\partial V_n(\M r)}{\partial \M e_{east}} \M e_{east}\right)\right],
\end{equation}
where $\gamma$ is the normal gravity, the load Love and Shida numbers $h_n,l_n$ are given by
\configFile{inputfileDeformationLoadLoveNumber}{matrix} and the load Love numbers $k_n$ are given by
\configFile{inputfilePotentialLoadLoveNumber}{matrix}. The $V_n$ are the spherical harmonics expansion of
the full time variable gravitational potential (potential of the loading mass + deformation potential):
\begin{equation}
V(\M r) = \sum_{n=0}^\infty V_n(\M r).
\end{equation}
Deformations due to Earth tide and due to polar tides are computed using the IERS conventions.
Eq.~\eqref{eq:displacement} is not used in these cases.

The \config{outputfileTimeSeries} is an \file{instrument file}{instrument}, MISCVALUES.
The data columns contain the deformation of each station in $x,y,z$ in a global terrestrial
reference frame or alternatively in a local elliposidal frame (north, east, up)
if \config{localReferenceFrame} is set.
)";

/***********************************************/

#include "programs/program.h"
#include "base/planets.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "classes/grid/grid.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/tides/tides.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"

/***** CLASS ***********************************/

/** @brief Generates a time series of station displacements.
* @ingroup programsGroup */
class Gravityfield2DisplacementTimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2DisplacementTimeSeries, SINGLEPROCESS, "generates a time series of station displacements", Gravityfield, TimeSeries)
GROOPS_RENAMED_PROGRAM(DisplacementTimeSeries, Gravityfield2DisplacementTimeSeries, date2time(2020, 2, 9))

/***********************************************/

void Gravityfield2DisplacementTimeSeries::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName         outputName;
    FileName         deformationName, potentialName;
    GridPtr          grid;
    TimeSeriesPtr    timeSeries;
    GravityfieldPtr  gravityfield;
    TidesPtr         tides;
    EarthRotationPtr earthRotation;
    EphemeridesPtr   ephemerides;
    Bool             removeMean, useLocalFrame;

    readConfig(config, "outputfileTimeSeries",               outputName,      Config::MUSTSET,  "", "x,y,z [m] per station");
    readConfig(config, "grid",                               grid,            Config::MUSTSET,  "",  "station list");
    readConfig(config, "timeSeries",                         timeSeries,      Config::MUSTSET,  "", "");
    readConfig(config, "gravityfield",                       gravityfield,    Config::DEFAULT,  "", "");
    readConfig(config, "tides",                              tides,           Config::DEFAULT,  "", "");
    readConfig(config, "earthRotation",                      earthRotation,   Config::MUSTSET,  "", "");
    readConfig(config, "ephemerides",                        ephemerides,     Config::OPTIONAL, "jpl", "");
    readConfig(config, "inputfileDeformationLoadLoveNumber", deformationName, Config::MUSTSET,  "{groopsDataDir}/loading/deformationLoveNumbers_CM_Gegout97.txt", "");
    readConfig(config, "inputfilePotentialLoadLoveNumber",   potentialName,   Config::OPTIONAL, "{groopsDataDir}/loading/loadLoveNumbers_Gegout97.txt", "if full potential is given and not only loading potential");
    readConfig(config, "removeMean",                         removeMean,      Config::DEFAULT,  "0", "remove the temporal mean of each coordinate");
    readConfig(config, "localReferenceFrame",                useLocalFrame,   Config::DEFAULT,  "0", "local left handed reference frame (north, east, up)");
    if(isCreateSchema(config)) return;

    std::vector<Time>     times  = timeSeries->times();
    std::vector<Vector3d> points = grid->points();

    // deformation load love numbers
    // -----------------------------
    Matrix love;
    readFileMatrix(deformationName, love);
    Vector hn = love.column(0);
    Vector ln = love.column(1);

    // models contain the total mass (loading mass & deformation mass effect)
    if(!potentialName.empty())
    {
      Vector kn;
      readFileMatrix(potentialName, kn);
      for(UInt n=2; n<std::min(kn.rows(), hn.rows()); n++)
        hn(n) /= (1.+kn(n));
      for(UInt n=2; n<std::min(kn.rows(), ln.rows()); n++)
        ln(n) /= (1.+kn(n));
    }

    // normal gravity
    // --------------
    Vector gravity(points.size());
    for(UInt i=0; i<points.size(); i++)
      gravity(i) = Planets::normalGravity(points.at(i));

    // Earth rotation
    // --------------
    logStatus<<"Compute Earth rotation"<<Log::endl;
    std::vector<Rotary3d> rotEarth(times.size());
    Single::forEach(times.size(), [&](UInt i)
    {
      rotEarth.at(i) = earthRotation->rotaryMatrix(times.at(i));
    });

    // displacements
    // -------------
    logStatus<<"compute station displacements"<<Log::endl;
    std::vector< std::vector<Vector3d> > disp(points.size());
    for(UInt k=0; k<points.size(); k++)
      disp.at(k).resize(times.size());
    gravityfield->deformation(times, points, gravity, hn, ln, disp);
    tides->deformation(times, points, rotEarth, earthRotation, ephemerides, gravity, hn, ln, disp);

    // local frame
    // -----------
    if(useLocalFrame)
    {
      logStatus<<"rotate into local left handed reference frame (north, east, up)"<<Log::endl;
      for(UInt k=0; k<points.size(); k++)
      {
        const Transform3d lnof2trf = localNorthEastUp(points.at(k));
        for(UInt i=0; i<times.size(); i++)
          disp.at(k).at(i) = lnof2trf.inverseTransform(disp.at(k).at(i));
      }
    }

    // sort into matrix
    // ----------------
    Matrix A(times.size(), 1+3*points.size());
    for(UInt i=0; i<times.size(); i++)
    {
      A(i,0) = times.at(i).mjd();
      for(UInt k=0; k<points.size(); k++)
      {
        A(i,3*k+1) =  disp.at(k).at(i).x();
        A(i,3*k+2) =  disp.at(k).at(i).y();
        A(i,3*k+3) =  disp.at(k).at(i).z();
      }
    }

    // remove temporal mean
    // --------------------
    if(removeMean)
    {
      logStatus<<"remove mean"<<Log::endl;
      for(UInt k=1; k<A.columns(); k++)
        A.column(k) -= mean(A.column(k));
    }

    // write file
    // ----------
    if(!outputName.empty())
    {
      logStatus<<"write time series to file <"<<outputName<<">"<<Log::endl;
      InstrumentFile::write(outputName, Arc(times, A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
