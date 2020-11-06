/***********************************************/
/**
* @file gravityfieldCovariancesPropagation2GriddedData.cpp
*
* @brief Covariances of values of a gravity field on a grid.
*
* @author Torsten Mayer-Guerr
* @date 2008-11-06
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the covariance between a source point given
by longitude/latitude (\config{L}, \config{B}) and the points of a \configClass{grid}{gridType}
in terms of the functional given by \configClass{kernel}{kernelType} from the variance-covariance
matrix of a \configClass{gravityfield}{gravityfieldType} evaluated at \config{time}.

If \config{computeCorrelation} is set, the program returns the correlation according to
\begin{equation}
r_{ij} = \frac{\sigma_{ij}}{\sigma_i \sigma_j}
\end{equation}
in the range of [-1, 1] instead ofthe covariance.

See also \program{Gravityfield2GridCovarianceMatrix}, \program{GravityfieldVariancesPropagation2GriddedData}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Covariances of values of a gravity field on a grid.
* @ingroup programsGroup */
class GravityfieldCovariancesPropagation2GriddedData
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GravityfieldCovariancesPropagation2GriddedData, PARALLEL, "covariances of values of a gravity field on a grid", Gravityfield, Grid)

/***********************************************/

void GravityfieldCovariancesPropagation2GriddedData::run(Config &config)
{
  try
  {
    FileName        fileNameGrid;
    GridPtr         grid;
    KernelPtr       kernel;
    GravityfieldPtr gravityfield;
    Time            time;
    Angle           L, B;
    Double          height;
    Double          a, f;
    Bool            calcCorrelation;

    readConfig(config, "outputfileGriddedData", fileNameGrid,    Config::MUSTSET,  "", "gridded data file containing the covariance betwenn source point and grid points");
    readConfig(config, "grid",                  grid,            Config::MUSTSET,  "", "");
    readConfig(config, "kernel",                kernel,          Config::MUSTSET, "", "functional");
    readConfig(config, "gravityfield",          gravityfield,    Config::MUSTSET,  "", "");
    readConfig(config, "time",                  time,            Config::OPTIONAL, "",  "at this time the gravity field will be evaluated");
    readConfig(config, "L",                     L,               Config::DEFAULT,  "0", "longitude of variance point");
    readConfig(config, "B",                     B,               Config::DEFAULT,  "0", "latitude of variance point");
    readConfig(config, "height",                height,          Config::DEFAULT,  "0", "ellipsoidal height of source point");
    readConfig(config, "computeCorrelation",    calcCorrelation, Config::DEFAULT,  "0", "compute correlations instead of covariances");
    readConfig(config, "R",                     a,               Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",     f,               Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // Create grid
    // -----------
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();
    Vector3d              point0 = Ellipsoid(a,f)(L,B,height);

    // Compute covariances
    // -------------------
    logStatus<<"calculate covariances on grid"<<Log::endl;
    std::vector<Double> field(points.size());
    Parallel::forEach(field, [&](UInt i) {return gravityfield->covariance(time, point0, points.at(i), *kernel);});

    std::vector<Double> sigma(points.size());
    if(calcCorrelation)
    {
      logStatus<<"calculate standard deviations on grid"<<Log::endl;
      Parallel::forEach(sigma, [&](UInt i) {return sqrt(gravityfield->variance(time, points.at(i), *kernel));});
    }

    if(Parallel::isMaster())
    {
      if(calcCorrelation)
      {
        Double sigma0 = std::sqrt(gravityfield->variance(time, point0, *kernel));
        for(UInt i=0; i<field.size(); i++)
          field.at(i) /= (sigma0*sigma.at(i));
      }

      logStatus<<"save values to file <"<<fileNameGrid<<">"<<Log::endl;
      GriddedData griddedData(Ellipsoid(a,f), points, areas, {field});
      writeFileGriddedData(fileNameGrid, griddedData);
      MiscGriddedData::printStatistics(griddedData);
    } // if(Parallel::isMaster())
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
