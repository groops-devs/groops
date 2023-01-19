/***********************************************/
/**
* @file programTemplate.cpp
*
* @brief Short description.
*
* 1. rename class ProgramTemplate
* 2. give a short description of the program
* 3. add this file to the sources list (Makefile.source or sourcesCXX.txt)
*
* @author Author
* @author Second author
* @date yyyy-mm-dd
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Put documenation here...

References to other programs with \program{Gravityfield2GriddedData}.

Describe elements with \config{factor}.

Description of config elements with links with \configClass{configElement}{gravityfieldType}
contains the name of config element and name of the class.

Description of file formats: The \configFile{outputfileEOP}{matrix} is a matrix \ldots
or alternatively: The \config{outputfileEOP} is a \file{matrix}{matrix}  \ldots
(file types: admittance, arcList, doodsonEarthOrientationParameter, doodsonHarmonic, earthOrientationParameter, earthTide,
gnssAntennaDefinition, gnssIonosphereMaps, gnssReceiverDefinition, gnssSignalBias, platform, gradioAccelerometerCalibration,
griddedData, instrument, matrix, meanPolarMotion, normalEquation, oceanPoleTide, parameterName, polygon, potentialCoefficients,
satelliteModel, stringList, stringTable, tideGeneratingPotential, timeSplinesCovariance, timeSplinesGravityField, variationalEquation)

Description with links to subsections:
\configClass{gravityfield:potentialCoefficients}{gravityfieldType:potentialCoefficients}
Maybe a label must be added to the subclass docstring.

References to other parts of the documentation:
\reference{dataVariables}{general.parser:dataVariables}.

%Figures: \fig{position}{width}{fileName without .png}{label}{caption}
Example: \fig{!hb}{0.5}{regionalGeoidTopography}{fig:regionalGeoidTopographyTemplate}{Topography and geoid heights}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Short description.
* Detailed description.
* @ingroup programsGroup */
class ProgramTemplate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

// Give some tags of the program for description.
// The first tag should agree with the first part of the name/directory.
// For a full list of tags see program.h (GROOPS_TAGS),
// Conversion, Provisional, Deprecated
// Covariance, DoodsonHarmonics, Gnss, Grace, Gravityfield, Grid, Instrument, KalmanFilter, Matrix, Misc,
// Noise, NormalEquation, Orbit, Plot, PotentialCoefficients, Preprocessing, Residuals, Simulation, SpatialTimeSeries,
// Statistics, System, TimeSeries, TimeSplines, VariationalEquation
GROOPS_REGISTER_PROGRAM(ProgramTemplate, PARALLEL/*or SINGLEPROCESS*/, "short description.", Provisional/*, further tags, ...*/)

/***********************************************/

void ProgramTemplate::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOut;
    UInt     count;

    readConfig(config, "outputfile", fileNameOut, Config::MUSTSET,  "", "description");
    readConfig(config, "count",      count,       Config::DEFAULT,  "100" /* default value */, "description");
    if(isCreateSchema(config)) return;

    // a short function (see c++ lambda functions)
    auto compute = [](UInt i) {return i*i;};

    // non-parallel loop
    // -----------------
    logStatus<<"starting non parallel loop"<<Log::endl;
    std::vector<Double> result(count);
    Single::forEach(count, [&](UInt i)
    {
      result.at(i) = compute(i);
    });

    // same with paralleization
    // ------------------------
    logStatus<<"starting parallel loop"<<Log::endl;
    Parallel::forEach(result, [&](UInt i) {return compute(i);}, comm);
    // or insertion of code directly as alternative:
    Parallel::forEach(result, [](UInt i) {return i*i;}, comm);
    // for functions without return:
    Parallel::forEach(count, [&](UInt i) {compute(i);}, comm);

    // Write results
    // -------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"writing output file <"<<fileNameOut<<">"<<Log::endl;
      writeFileMatrix(fileNameOut, Vector(result));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
