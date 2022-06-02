/***********************************************/
/**
* @file orbit2MagneticField.cpp
*
* @brief Magentic field vector along orbit.
*
* @author Torsten Mayer-Guerr
* @date 2022-05-20
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the magentic field vector($x, y, z$ $[Tesla = kg/A/s^2]$ in CRF))
along an \file{orbit}{instrument} and writes it as \file{instrument file}{instrument} (MISCVALUES).
The data of \configFile{inputfileInstrument}{instrument} are appended as data columns to each epoch.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/magnetosphere/magnetosphere.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Magentic field vector along orbit.
* @ingroup programsGroup */
class Orbit2MagneticField
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Orbit2MagneticField, PARALLEL, "Thermospheric state along orbit.", Orbit, Instrument)

/***********************************************/

void Orbit2MagneticField::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName              fileNameOut, fileNameOrbit;
    std::vector<FileName> fileNamesInstrument;
    MagnetospherePtr      magnetosphere;
    EarthRotationPtr      earthRotation;

    readConfig(config, "outputfileMagneticField", fileNameOut,         Config::MUSTSET,  "", "instrument file (x,y,z in CRF [Tesla = kg/A/s^2]), ...)");
    readConfig(config, "inputfileOrbit",          fileNameOrbit,       Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrument",     fileNamesInstrument, Config::OPTIONAL, "", "data are appended to output file");
    readConfig(config, "magnetosphere",           magnetosphere,       Config::MUSTSET,  "",  "");
    readConfig(config, "earthRotation",           earthRotation,       Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =======================

    logStatus<<"computing thermospheric state"<<Log::endl;
    InstrumentFile orbitFile(fileNameOrbit);
    UInt dataCount = 3; // (x,y, z)
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
        Vector3d field = rotEarth.inverseRotate(magnetosphere->magenticFieldVector(orbit.at(i).time, rotEarth.rotate(orbit.at(i).position)));

        A(i, 1) = field.x();
        A(i, 2) = field.y();
        A(i, 3) = field.z();
      }

      UInt idx = 3;
      for(auto &file: instrumentFile)
      {
        Arc arc = file->readArc(arcNo);
        Arc::checkSynchronized({orbit, arc});
        Matrix B = arc.matrix();
        copy(B.column(1, B.columns()-1), A.column(idx, B.columns()-1));
        idx += B.columns()-1;
      }

      return Arc(orbit.times(), A);
    }, comm);

    // write results
    // -------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write magnetic field vector to file <"<<fileNameOut<<">"<<Log::endl;
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
