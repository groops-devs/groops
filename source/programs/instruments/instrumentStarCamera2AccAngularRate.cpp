/***********************************************/
/**
* @file instrumentStarCamera2AccAngularRate.cpp
*
* @brief Derivation of angular rates from rotations.
*
* @author Torsten Mayer-Guerr
* @date 2011-04-18
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program derivate from a time series of quaternions
a series of angular rates and angular accelerations.
The derivatives are computed by a polynomial interpolation
with \config{interpolationDegree} of the quaternions.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Derivation of angular rates from rotations.
* @ingroup programsGroup */
class InstrumentStarCamera2AccAngularRate
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentStarCamera2AccAngularRate, SINGLEPROCESS, "derivation of angular rates from rotations", Instrument)

/***********************************************/

void InstrumentStarCamera2AccAngularRate::run(Config &config)
{
  try
  {
    FileName angRateName, angAccName;
    FileName starCameraName;
    UInt     derivationDegree;

    readConfig(config, "outputfileAngularRate", angRateName,      Config::OPTIONAL, "", "[rad/s], VECTOR3D");
    readConfig(config, "outputfileAngularAcc",  angAccName,       Config::OPTIONAL, "", "[rad/s**2], VECTOR3D");
    readConfig(config, "inputfileStarCamera",   starCameraName,   Config::MUSTSET,  "", "");
    readConfig(config, "interpolationDegree",   derivationDegree, Config::DEFAULT,  "2", "derivation by polynomial interpolation of degree n");
    if(isCreateSchema(config)) return;

    // init interpolation coefficients
    // -------------------------------
    Polynomial polynomial(derivationDegree);

    // calculate angular rates & accelerations
    // ---------------------------------------
    logStatus<<"calculate angular rates & accelerations"<<Log::endl;
    InstrumentFile starCameraFile(starCameraName);
    std::list<Arc> angRateArcList, angAccArcList;

    logTimerStart;
    for(UInt arcNo=0; arcNo<starCameraFile.arcCount(); arcNo++)
    {
      logTimerLoop(arcNo, starCameraFile.arcCount());

      StarCameraArc scaArc = starCameraFile.readArc(arcNo);
      const Double  dt     = (scaArc.at(1).time - scaArc.at(0).time).seconds();

      // quaternions
      // -----------
      Matrix q   = scaArc.matrix().column(1,4);
      Matrix dq  = polynomial.derivative(dt, q);
      Matrix ddq = polynomial.derivative2nd(dt, q);

      // angular rate
      // ------------
      Vector3dArc angRateArc;
      for(UInt i=0; i<scaArc.size(); i++)
      {
        Vector3dEpoch angRate;
        angRate.time = scaArc.at(i).time;
        angRate.vector3d.x() = 2*(-q(i,1)*dq(i,0) + q(i,0)*dq(i,1) + q(i,3)*dq(i,2) - q(i,2)*dq(i,3));
        angRate.vector3d.y() = 2*(-q(i,2)*dq(i,0) - q(i,3)*dq(i,1) + q(i,0)*dq(i,2) + q(i,1)*dq(i,3));
        angRate.vector3d.z() = 2*(-q(i,3)*dq(i,0) + q(i,2)*dq(i,1) - q(i,1)*dq(i,2) + q(i,0)*dq(i,3));
        angRateArc.push_back(angRate);
      }
      angRateArcList.push_back(angRateArc);

      // angular accelerations
      // ---------------------
      if(!angAccName.empty())
      {
        Vector3dArc angAccArc;
        for(UInt i=0; i<scaArc.size(); i++)
        {
          Vector3dEpoch angAcc;
          angAcc.time = scaArc.at(i).time;
          angAcc.vector3d.x() = 2*(-q(i,1)*ddq(i,0) + q(i,0)*ddq(i,1) + q(i,3)*ddq(i,2) - q(i,2)*ddq(i,3));
          angAcc.vector3d.y() = 2*(-q(i,2)*ddq(i,0) - q(i,3)*ddq(i,1) + q(i,0)*ddq(i,2) + q(i,1)*ddq(i,3));
          angAcc.vector3d.z() = 2*(-q(i,3)*ddq(i,0) + q(i,2)*ddq(i,1) - q(i,1)*ddq(i,2) + q(i,0)*ddq(i,3));
          angAccArc.push_back(angAcc);
        }
        angAccArcList.push_back(angAccArc);
      }
    }
    logTimerLoopEnd(starCameraFile.arcCount());

    // write data
    // ----------
    if(!angRateName.empty())
    {
      logStatus<<"write angular rate data to file <"<<angRateName<<">"<<Log::endl;
      InstrumentFile::write(angRateName, angRateArcList);
    }
    if(!angAccName.empty())
    {
      logStatus<<"write angular accelerations data to file <"<<angAccName<<">"<<Log::endl;
      InstrumentFile::write(angAccName, angAccArcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
