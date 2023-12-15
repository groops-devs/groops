/***********************************************/
/**
* @file orbit2EarthFixedOrbit.cpp
*
* @brief Rotate an orbit into a rotation earth fixed frame.
*
* @author Torsten Mayer-Guerr
* @date 2023-12-07
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Normally the orbits in GROOPS are given in the celestial reference frame (CRF) with the
origin in the center of mass (CoM). This program rotates the orbit with
\configClass{earthRotation}{earthRotationType} from CRF to the TRF.

To additionally tranform into the center of solid Earth (CE) frame (or center of Figure (CF)),
a correction can be applied by providing degree one coefficients of a
\configClass{gravityfield}{gravityfieldType} (e.g. ocean tides).

If \config{celestial2terrestrial} is set to no, the inverse transformation is applied.

See also \program{InstrumentRotate}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Rotate an orbit into a rotation earth fixed frame.
* @ingroup programsGroup */
class Orbit2EarthFixedOrbit
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Orbit2EarthFixedOrbit, PARALLEL, "rotate an orbit into a rotation earth fixed frame", Orbit)

/***********************************************/

void Orbit2EarthFixedOrbit::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName         fileNameOut, fileNameIn;
    EarthRotationPtr earthRotation;
    GravityfieldPtr  gravityfield;
    Bool             crf2trf;

    readConfig(config, "outputfileOrbit",       fileNameOut,   Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit",        fileNameIn,    Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",         earthRotation, Config::MUSTSET,  "", "transformation from CRF to TRF");
    readConfig(config, "gravityfield",          gravityfield,  Config::DEFAULT,  R"([{"tides": {"tides": {"doodsonHarmonicTide": {"minDegree":1, "maxDegree":1}}}}])", "degree 1 fluid mantle for CM2CE correction");
    readConfig(config, "celestial2terrestrial", crf2trf,       Config::DEFAULT,  "1", "yes: crf->trf, no: trf->crf");
    if(isCreateSchema(config)) return;

    logStatus<<"read orbit file <"<<fileNameIn<<">"<<Log::endl;
    InstrumentFile orbitFile(fileNameIn);
    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&] (UInt arcNo)
    {
      OrbitArc orbit = orbitFile.readArc(arcNo);
      for(UInt i=0; i<orbit.size(); i++)
      {
        const SphericalHarmonics harmonics = gravityfield->sphericalHarmonics(orbit.at(i).time, 1, 1);
        const Vector coeff = harmonics.x(); // [c00, c10, c11, s11]
        const Vector3d cm2ceCorrection = std::sqrt(3.) * harmonics.R() * Vector3d(coeff(2), coeff(3), coeff(1));
        const Rotary3d rot = earthRotation->rotaryMatrix(orbit.at(i).time);;

        if(crf2trf)
        {
          orbit.at(i).position = rot.rotate(orbit.at(i).position) + cm2ceCorrection;
          if(orbit.at(i).velocity.r() > 0)
          {
            const Vector3d omega = rot.rotate(earthRotation->rotaryAxis(orbit.at(i).time));
            orbit.at(i).velocity = rot.rotate(orbit.at(i).velocity) - crossProduct(omega, orbit.at(i).position);
            if(orbit.at(i).acceleration.r() > 0)
              orbit.at(i).acceleration = rot.rotate(orbit.at(i).acceleration)
                                       - crossProduct(omega, crossProduct(omega, orbit.at(i).position))
                                       - 2*crossProduct(omega, orbit.at(i).velocity);
          }
        }
        else
        {
          if(orbit.at(i).velocity.r() > 0)
          {
            const Vector3d omega = rot.rotate(earthRotation->rotaryAxis(orbit.at(i).time));
            if(orbit.at(i).acceleration.r() > 0)
              orbit.at(i).acceleration = rot.inverseRotate(orbit.at(i).acceleration
                                                        + crossProduct(omega, crossProduct(omega, orbit.at(i).position))
                                                        + 2*crossProduct(omega, orbit.at(i).velocity));
            orbit.at(i).velocity = rot.inverseRotate(orbit.at(i).velocity + crossProduct(omega, orbit.at(i).position));
          }
          orbit.at(i).position = rot.inverseRotate(orbit.at(i).position - cm2ceCorrection);
        }
      }
      return orbit;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write rotated orbit to file <"<<fileNameOut<<">"<<Log::endl;
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
