/***********************************************/
/**
* @file observationMisc.cpp
*
* @brief Right hand sides and calibration parameter.
*
* @author Torsten Mayer-Guerr
* @date 2008-07-28
*
*/
/***********************************************/

#define DOCSTRING_PodRightSide
#define DOCSTRING_SstRightSide
#define DOCSTRING_SggRightSide

#include "base/import.h"
#include "config/configRegister.h"
#include "files/fileInstrument.h"
#include "classes/tides/tides.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/miscAccelerations/miscAccelerations.h"
#include "observationMisc.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(PodRightSide, "podRightSideType")
GROOPS_READCONFIG_CLASS(PodRightSide, "podRightSideType")

PodRightSide::PodRightSide(Config &config, const std::string &name)
{
  try
  {
    FileName orbitName, accName;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "inputfileOrbit",         orbitName, Config::MUSTSET,  "", "kinematic positions of satellite as observations");
    readConfig(config, "inputfileAccelerometer", accName,   Config::OPTIONAL, "", "non-gravitational forces in satellite reference frame");
    readConfig(config, "forces",                 forces,    Config::MUSTSET,  "", "");
    endSequence(config);
    if(isCreateSchema(config)) return;

    orbitFile         = InstrumentFile::newFile(orbitName);
    accelerometerFile = InstrumentFile::newFile(accName);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(SstRightSide, "sstRightSideType")
GROOPS_READCONFIG_CLASS(SstRightSide, "sstRightSideType")

SstRightSide::SstRightSide(Config &config, const std::string &name)
{
  try
  {
    std::vector<FileName> sstName;
    FileName              orbit1Name, acc1Name;
    FileName              orbit2Name, acc2Name;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "inputfileSatelliteTracking", sstName,      Config::OPTIONAL, "", "ranging observations and corrections");
    readConfig(config, "inputfileOrbit1",            orbit1Name,   Config::OPTIONAL, "", "kinematic positions of satellite A as observations");
    readConfig(config, "inputfileOrbit2",            orbit2Name,   Config::OPTIONAL, "", "kinematic positions of satellite B as observations");
    readConfig(config, "inputfileAccelerometer1",    acc1Name,     Config::OPTIONAL, "", "non-gravitational forces in satellite reference frame A");
    readConfig(config, "inputfileAccelerometer2",    acc2Name,     Config::OPTIONAL, "", "non-gravitational forces in satellite reference frame B");
    readConfig(config, "forces",                     forces,       Config::MUSTSET,  "", "");
    endSequence(config);
    if(isCreateSchema(config)) return;

    sstFile.resize(sstName.size());
    for(UInt i=0; i<sstName.size(); i++)
      sstFile.at(i) = InstrumentFile::newFile(sstName.at(i));
    orbit1File = InstrumentFile::newFile(orbit1Name);
    orbit2File = InstrumentFile::newFile(orbit2Name);
    accelerometer1File  = InstrumentFile::newFile(acc1Name);
    accelerometer2File  = InstrumentFile::newFile(acc2Name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(SggRightSide, "sggRightSideType")
GROOPS_READCONFIG_CLASS(SggRightSide, "sggRightSideType")

SggRightSide::SggRightSide(Config &config, const std::string &name)
{
  try
  {
    FileName              fileNameGradiometer;
    std::vector<FileName> fileNameReference;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "inputfileGradiometer",          fileNameGradiometer, Config::MUSTSET,  "", "observed gravity gradients");
    readConfig(config, "inputfileReferenceGradiometer", fileNameReference,   Config::OPTIONAL, "", "precomputed gradients at orbit positions");
    readConfig(config, "referencefield",                referencefield,      Config::DEFAULT,  "", "");
    readConfig(config, "tides",                         tides,               Config::DEFAULT,  "", "");
    endSequence(config);
    if(isCreateSchema(config)) return;

    gradiometerFile = InstrumentFile::newFile(fileNameGradiometer);
    referenceFile.resize(fileNameReference.size());
    for(UInt i=0; i<fileNameReference.size(); i++)
      referenceFile.at(i) = InstrumentFile::newFile(fileNameReference.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
