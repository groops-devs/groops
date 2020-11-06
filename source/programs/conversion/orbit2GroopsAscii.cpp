/***********************************************/
/**
* @file orbit2GroopsAscii.cpp
*
* @brief write orbit positions to TVGOGO ASCII format.
*
* @author Norbert Zehentner
* @date 2013-08-06
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert groops orbits and corresponding covariance information to ASCII format.
The format is used to publish TUG orbits. It contains a two line header
with a short description of the orbit defined in \config{firstLine}.
The orbit is rotated to the Earth fixed frame (TRF) with \configClass{earthRotation}{earthRotationType} and given as one line per epoch.
The epoch lines containe time [MJD GPS time], position x, y and z [m], and the epoch covariance xx, yy, zz, xy, xz and yz [$m^2$].

See also \program{GroopsAscii2Orbit}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief write orbit positions into TVGOGO ASCII file.
* @ingroup programsConversionGroup */
class Orbit2GroopsAscii
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Orbit2GroopsAscii, SINGLEPROCESS, "write orbit positions to TUG ASCII format", Conversion, Orbit, Covariance, Instrument)
GROOPS_RENAMED_PROGRAM(Orbit2Ascii, Orbit2GroopsAscii, date2time(2020, 6, 14))

/***********************************************/

void Orbit2GroopsAscii::run(Config &config)
{
  try
  {
    FileName         outName, inName, inNameCov;
    std::string      textLine1, textDatum;
    EarthRotationPtr earthRotation;

    readConfig(config, "outputfile",          outName,        Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit",      inName,         Config::MUSTSET,  "", "");
    readConfig(config, "inputfileCovariance", inNameCov,      Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",       earthRotation,  Config::MUSTSET,  "", "");
    readConfig(config, "firstLine",           textLine1,      Config::OPTIONAL, "", "Text for first line");
    if(isCreateSchema(config)) return;

    logStatus<<"read orbit file <"<<inName<<">"<<Log::endl;
    OrbitArc        orbit = InstrumentFile::read(inName);
    Covariance3dArc cov   = InstrumentFile::read(inNameCov);
    Arc::checkSynchronized({orbit, cov});

    logStatus<<"write file <"<<outName<<">"<<Log::endl;
    OutFile file(outName);
    file<<"# "<<textLine1<<textDatum<<std::endl;
    file<<"# MJD GPS time [days]        x [m]           y [m]           z [m]        xx [m^2]        yy [m^2]        zz [m^2]        xy [m^2]        xz [m^2]        yz [m^2]"<<std::endl;

    logTimerStart;
    for(UInt i=0; i<orbit.size(); i++)
    {
      logTimerLoop(i, orbit.size());
      const Rotary3d rot      = earthRotation->rotaryMatrix(orbit.at(i).time);
      const Vector3d position = rot.rotate(orbit.at(i).position);
      Tensor3d cv;
      if(cov.size())
        cv = rot.rotate(cov.at(i).covariance);

      file<<orbit.at(i).time%"%21.15f"s;
      file<<position.x()%" %15.4f"s<<position.y()%" %15.4f"s<<position.z()%" %15.4f"s;
      file<<cv.xx()%" %15.8e"s<<cv.yy()%" %15.8e"s<<cv.zz()%" %15.8e"s<<cv.xy()%" %15.8e"s<<cv.xz()%" %15.8e"s<<cv.yz()%" %15.8e"s;
      file<<std::endl;
    }
    logTimerLoopEnd(orbit.size());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
