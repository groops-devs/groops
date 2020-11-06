/***********************************************/
/**
* @file orbit2ThermosphericState.cpp
*
* @brief Thermospheric state along orbit.
*
* @author Torsten Mayer-Guerr
* @date 2020-03-07
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the thermosperic state (density, temperature, wind)
based on emprical models along an \file{orbit}{instrument}
and writes it as \file{instrument file}{instrument} (MISCVALUES).
The wind is given in an celestial reference frame (CRF).
The data of \configFile{inputfileInstrument}{instrument} are appended as values to each point.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/thermosphere/thermosphere.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Thermospheric state along orbit.
* @ingroup programsGroup */
class Orbit2ThermosphericState
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Orbit2ThermosphericState, PARALLEL, "Thermospheric state along orbit.", Orbit, Instrument)
GROOPS_RENAMED_PROGRAM(InstrumentOrbit2ThermosphericState, Orbit2ThermosphericState, date2time(2020, 05, 25))

/***********************************************/

void Orbit2ThermosphericState::run(Config &config)
{
  try
  {
    FileName              fileNameOut, fileNameOrbit;
    std::vector<FileName> fileNamesInstrument;
    ThermospherePtr       thermosphere;
    EarthRotationPtr      earthRotation;

    readConfig(config, "outputfileThermosphericState", fileNameOut,         Config::MUSTSET,  "", "density, temperature, wind(x,y,z)(CRF)");
    readConfig(config, "inputfileOrbit",               fileNameOrbit,       Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrument",          fileNamesInstrument, Config::OPTIONAL, "", "data are appended to output file");
    readConfig(config, "thermosphere",                 thermosphere,        Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",                earthRotation,       Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =======================

    logStatus<<"computing thermospheric state"<<Log::endl;
    InstrumentFile orbitFile(fileNameOrbit);
    UInt dataCount = 6; // time, density, temperature, wind(x,yz)
    std::vector<InstrumentFilePtr> instrumentFile;
    for(auto &fileName : fileNamesInstrument)
    {
      instrumentFile.push_back(InstrumentFile::newFile(fileName));
      InstrumentFile::checkArcCount({orbitFile, *instrumentFile.back()});
      dataCount += instrumentFile.back()->dataCount(TRUE/*mustDefined*/);
    }

    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&] (UInt arcNo)
    {
      const OrbitArc orbit = orbitFile.readArc(arcNo);
      Matrix A(orbit.size(), dataCount);
      for(UInt i=0; i<orbit.size(); i++)
      {
        const Rotary3d rotEarth = earthRotation->rotaryMatrix(orbit.at(i).time);
        Double   density, temperature;
        Vector3d wind;
        thermosphere->state(orbit.at(i).time, rotEarth.rotate(orbit.at(i).position), density, temperature, wind);
        wind = rotEarth.inverseRotate(wind) + crossProduct(earthRotation->rotaryAxis(orbit.at(i).time), orbit.at(i).position);

        A(i, 1) = density;
        A(i, 2) = temperature;
        A(i, 3) = wind.x();
        A(i, 4) = wind.y();
        A(i, 5) = wind.z();
      }

      UInt idx = 6;
      for(auto &file: instrumentFile)
      {
        Arc arc = file->readArc(arcNo);
        Arc::checkSynchronized({orbit, arc});
        Matrix B = arc.matrix();
        copy(B.column(1, B.columns()-1), A.column(idx, B.columns()-1));
        idx += B.columns()-1;
      }

      return Arc(orbit.times(), A);
    });

    // write results
    // -------------
    if(Parallel::isMaster())
    {
      logStatus<<"write thermospheric state to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
