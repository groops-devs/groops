/***********************************************/
/**
* @file doodsonHarmonics2GriddedAmplitudeAndPhase.cpp
*
* @brief Amplitude and phase of a harmonic tidal constituent on a grid.
*
* @author Torsten Mayer-Guerr
* @date 2008-08-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads a \configFile{inputfileDoodsonHarmonics}{doodsonHarmonic} and evaluates a single tidal
constituent selected by \config{dooddson} (Doodson number or Darwin´s name, e.g. 255.555 or M2).
This program computes the amplitude and phase from the cos and sin coefficients on
a given \configClass{grid}{gridType}. The type of functional (e.g gravity anomalies or geoid heights)
can be choosen with \configClass{kernel}{kernelType}.
The values will be saved together with points expressed as ellipsoidal coordinates (longitude, latitude, height)
based on a reference ellipsoid with parameters \config{R} and \config{inverseFlattening}.
To visualize the results use \program{PlotMap}.

\fig{!hb}{1.}{doodsonHarmonics2GriddedAmplitudeAndPhase}{fig:doodsonHarmonics2GriddedAmplitudeAndPhase}{M2 amplitude and phase of FES2014b.}
)";

/***********************************************/

#include "programs/program.h"
#include "base/doodson.h"
#include "files/fileGriddedData.h"
#include "files/fileDoodsonHarmonic.h"
#include "classes/grid/grid.h"
#include "classes/kernel/kernel.h"
#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilter.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Amplitude and phase of a harmonic tidal constituent on a grid.
* @ingroup programsGroup */
class DoodsonHarmonics2GriddedAmplitudeAndPhase
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(DoodsonHarmonics2GriddedAmplitudeAndPhase, PARALLEL, "amplitude and phase of a harmonic tidal constituent on a grid.", DoodsonHarmonics, Grid)

/***********************************************/

void DoodsonHarmonics2GriddedAmplitudeAndPhase::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName  fileNameGrid;
    FileName  inputName;
    Doodson   doodson;
    Double    factor;
    SphericalHarmonicsFilterPtr filter;
    GridPtr   grid;
    KernelPtr kernel;
    Double    a, f;
    UInt      minDegree, maxDegree = INFINITYDEGREE;

    readConfig(config, "outputfileGrid",            fileNameGrid, Config::MUSTSET,   "",    "ampl, phase [-pi,pi], cos, sin");
    readConfig(config, "inputfileDoodsonHarmonics", inputName,    Config::MUSTSET,   "{groopsDataDir}/tides/", "");
    readConfig(config, "doodson",                   doodson,      Config::MUSTSET,   "",    "tidal constituent");
    readConfig(config, "filter",                    filter,       Config::DEFAULT,   "",    "");
    readConfig(config, "grid",                      grid,         Config::MUSTSET,   "",    "");
    readConfig(config, "kernel",                    kernel,       Config::MUSTSET,   "",    "");
    readConfig(config, "minDegree",                 minDegree,    Config::DEFAULT,   "0",   "");
    readConfig(config, "maxDegree",                 maxDegree,    Config::OPTIONAL, "",    "");
    readConfig(config, "factor",                    factor,       Config::DEFAULT,   "1.0", "the values on grid are multplied by this factor");
    readConfig(config, "R",                         a,            Config::DEFAULT,   STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",         f,            Config::DEFAULT,   STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // read ocean tide file
    // --------------------
    logStatus<<"read doodson harmonics file <"<<inputName<<">"<<Log::endl;
    SphericalHarmonics harmCos, harmSin;
    {
      DoodsonHarmonic d;
      readFileDoodsonHarmonic(inputName, d);

      Bool found = FALSE;
      for(UInt i=0; i<d.doodson.size(); i++)
        if(doodson == d.doodson.at(i))
        {
          harmCos = filter->filter(SphericalHarmonics(d.GM, d.R, d.cnmCos.at(i), d.snmCos.at(i)).get(maxDegree, minDegree));
          harmSin = filter->filter(SphericalHarmonics(d.GM, d.R, d.cnmSin.at(i), d.snmSin.at(i)).get(maxDegree, minDegree));
          found = TRUE;
          break;
        }
      if(!found)
        throw(Exception("tide "+doodson.name()+" not included in <"+inputName.str()+">"));
    }

    // create grid values
    // ------------------
    logStatus<<"create values on grid (cos, sin)"<<Log::endl;
    const std::vector<Vector3d> points = grid->points();
    std::vector<std::vector<Double>> values(4); // ampl, phase, cos, sin
    values.at(0).resize(points.size());
    values.at(1).resize(points.size());
    values.at(2) = MiscGriddedData::synthesisSphericalHarmonics(harmCos, points, kernel, comm);
    values.at(3) = MiscGriddedData::synthesisSphericalHarmonics(harmSin, points, kernel, comm);

    if(Parallel::isMaster(comm))
    {
      for(UInt i=0; i<points.size(); i++)
      {
        values.at(2).at(i) *= factor; // cos
        values.at(3).at(i) *= factor; // sin
        values.at(0).at(i)  = std::sqrt(pow(values.at(2).at(i),2) + pow(values.at(3).at(i),2)); // ampl
        values.at(1).at(i)  = 0; // phase
        if(values.at(0).at(i)>0)
          values.at(1).at(i) = std::atan2(values.at(3).at(i), values.at(2).at(i));
      }

      logStatus<<"write grid <"<<fileNameGrid<<">"<<Log::endl;
      GriddedData griddedData(Ellipsoid(a,f), points, grid->areas(), values);
      writeFileGriddedData(fileNameGrid, griddedData);
      MiscGriddedData::printStatistics(griddedData);
    } // if(Parallel::isMaster(comm))
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
