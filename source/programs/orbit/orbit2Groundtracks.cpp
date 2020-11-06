/***********************************************/
/**
* @file orbit2Groundtracks.cpp
*
* @brief Groundtracks of a satellite as gridded data.
*
* @author Torsten Mayer-Guerr
* @date 2003-05-21
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program write \file{satellites positions}{instrument} as \file{gridded data}{griddedData}
(\config{outputfileTrackGriddedData}) in a terrestrial reference frame. The points are expressed as ellipsoidal coordinates
(longitude, latitude, height) based on a reference ellipsoid with parameters \config{R} and
\config{inverseFlattening}. The orbit data are given in the celestial frame so \configClass{earthRotation}{earthRotationType}
is needed to transform the data into the terrestrial frame.
The data of \configFile{inputfileInstrument}{instrument} are appended as values to each point.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileGriddedData.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Groundtracks of a satellite as gridded data.
* @ingroup programsGroup */
class Orbit2Groundtracks
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Orbit2Groundtracks, SINGLEPROCESS, "groundtracks of a satellite as gridded data.", Orbit, Instrument, Grid)
GROOPS_RENAMED_PROGRAM(ArcGroundtracks, Orbit2Groundtracks, date2time(2020, 05, 25))

/***********************************************/

void Orbit2Groundtracks::run(Config &config)
{
  try
  {
    FileName fileNameOut, fileNameOrbit;
    std::vector<FileName> fileNamesInstrument;
    EarthRotationPtr      earthRotation;
    Double                a, f;

    readConfig(config, "outputfileGriddedData", fileNameOut,         Config::MUSTSET,  "", "positions as gridded data");
    readConfig(config, "inputfileOrbit",        fileNameOrbit,       Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrument",   fileNamesInstrument, Config::OPTIONAL, "", "values at grid points");
    readConfig(config, "earthRotation",         earthRotation,       Config::MUSTSET,  "", "transformation from CRF to TRF");
    readConfig(config, "R",                     a,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",     f,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    // positions
    // ---------
    logStatus<<"read orbit file <"<<fileNameOrbit<<">"<<Log::endl;
    OrbitArc orbit = InstrumentFile::read(fileNameOrbit);
    std::vector<Vector3d> points;
    for(UInt i=0; i<orbit.size(); i++)
      points.push_back(earthRotation->rotaryMatrix(orbit.at(i).time).rotate(orbit.at(i).position));

    // values
    // ------
    std::vector<std::vector<Double>> values;
    for(auto &fileName : fileNamesInstrument)
    {
      logStatus<<"read instrument file <"<<fileName<<">"<<Log::endl;
      Arc arc = InstrumentFile::read(fileName);
      Arc::checkSynchronized({arc, orbit});

      const Matrix A = arc.matrix();
      for(UInt k=1; k<A.columns(); k++)
        values.push_back(Vector(A.column(k)));
    }

    // write gridded data
    // ------------------
    logStatus<<"write gridded data <"<<fileNameOut<<">"<<Log::endl;
    writeFileGriddedData(fileNameOut, GriddedData(Ellipsoid(a,f), points, {}, values));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
