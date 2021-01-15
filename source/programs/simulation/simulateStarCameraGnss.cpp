/***********************************************/
/**
* @file simulateStarCameraGnss.cpp
*
* @brief Simulates star camera data for a GNSS satellite based on an attitude model.
*
* @author Sebastian Strasser
* @date 2020-11-02
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates \file{star camera}{instrument} measurements at each satellite position
of \configFile{inputfileOrbit}{instrument}.
The resulting rotation matrices rotate from body frame to inertial frame. The body frame refers
to the IGS-specific (not the manufacturer-specific) body frame, as described by
\href{https://doi.org/10.1016/j.asr.2015.06.019}{Montenbruck et al. (2015)}.
The \configFile{inputfileOrbit}{instrument} must contain velocities
(use \program{OrbitAddVelocityAndAcceleration} if needed).

Information about the attitude mode(s) used by the GNSS satellite may be provided via
\configFile{inputfileAttitudeInfo}{instrument}. This file can be created with
\program{GnssAttitudeInfoCreate}. It contains one or more time-dependent entries,
each defining the default attitude mode, the attitude modes used around orbit noon and
midnight, and some parameters required by the various modes.
If no \configFile{inputfileAttitudeInfo}{instrument} is selected, the program defaults
to a nominal yaw-steering attitude model.
A sufficiently high \config{modelingResolution} ensures that the attitude behavior is modeled properly
at all times.

The attitude behavior is defined by the respective mode. Here is a list of the supported
modes with a brief explanation and references:
\begin{itemize}
\item \textbf{nominalYawSteering}:
      Yaw to keep solar panels aligned to Sun (e.g. most GNSS satellites outside eclipse) [1]
\item \textbf{orbitNormal}:
      Keep fixed yaw angle, for example point X-axis in flight direction (e.g. BDS-2G, BDS-3G, QZS-2G) [1]
\item \textbf{catchUpYawSteering}:
      Yaw at maximum yaw rate to catch up to nominal yaw angle (e.g. GPS-* (noon), GPS-IIR (midnight)) [2, 3]
\item \textbf{shadowMaxYawSteeringAndRecovery}:
      Yaw at maximum yaw rate from shadow start to end, recover after shadow (e.g. GPS-IIA (midnight)) [2]
\item \textbf{shadowMaxYawSteeringAndStop}:
      Yaw at maximum yaw rate from shadow start until nominal yaw angle at shadow end is reached,
      then stop (e.g. GLO-M (midnight)) [4]
\item \textbf{shadowConstantYawSteering}:
      Yaw at constant yaw rate from shadow start to end (e.g. GPS-IIF (midnight)) [3]
\item \textbf{centeredMaxYawSteering}:
      Yaw at maximum yaw rate centered around noon/midnight (e.g. QZS-2I, GLO-M (noon)) [4, 8]
\item \textbf{smoothedYawSteering1}:
      Yaw based on an auxiliary Sun vector for a smooth yaw maneuver (e.g. GAL-1) [5]
\item \textbf{smoothedYawSteering2}:
      Yaw based on a modified yaw-steering law for a smooth yaw maneuver (e.g. GAL-2, BDS-3M, BDS-3I) [5, 6]
\item \textbf{betaDependentOrbitNormal}:
      Switch to orbit normal mode if below beta angle threshold (e.g. BDS-2M, BDS-2I, QZS-1) [7, 8]
\end{itemize}

\fig{!hb}{0.9}{gnssAttitudeModes}{fig:gnssAttitudeModes1}{Overview of attitude modes used by GNSS satellites}

See \program{GnssAttitudeInfoCreate} for more details on which satellite uses which attitude modes
and the required parameters for each mode.

References for the attitude modes:
\begin{enumerate}
\item \href{https://doi.org/10.1016/j.asr.2015.06.019}{Montenbruck et al. (2015)}
\item \href{https://doi.org/10.1007/s10291-008-0092-1}{Kouba (2009)}
\item \href{https://doi.org/10.1007/s10291-016-0562-9}{Kuang et al. (2017)}
\item \href{https://doi.org/10.1016/j.asr.2010.09.007}{Dilssner et al. (2011)}
\item \url{https://www.gsc-europa.eu/support-to-developers/galileo-satellite-metadata#3}
\item \href{https://doi.org/10.1007/s10291-018-0783-1}{Wang et al. (2018)}
\item \href{https://doi.org/10.1017/S0373463318000103}{Li et al. (2018)}
\item \url{https://qzss.go.jp/en/technical/qzssinfo/index.html}
\end{enumerate}
)";

/***********************************************/

#include "programs/program.h"
#include "base/kepler.h"
#include "base/polynomial.h"
#include "classes/eclipse/eclipse.h"
#include "classes/ephemerides/ephemerides.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Simulates star camera data for a GNSS satellite based on an attitude model.
* @ingroup programsGroup */
class SimulateStarCameraGnss
{
  enum AttitudeMode
  {
    NOMINAL_YAW_STEERING = 0,
    ORBIT_NORMAL = 1,                         // e.g. BDS-2G, BDS-3G, QZS-2G
    CATCH_UP_YAW_STEERING = 2,                // e.g. GPS-* (noon), GPS-IIR (midnight)
    SHADOW_MAX_YAW_STEERING_AND_RECOVERY = 3, // e.g. GPS-IIA (midnight)
    SHADOW_MAX_YAW_STEERING_AND_STOP = 4,     // e.g. GLO-M (midnight)
    SHADOW_CONSTANT_YAW_STEERING = 5,         // e.g. GPS-IIF (midnight)
    CENTERED_MAX_YAW_STEERING = 6,            // e.g. QZS-2I, GLO-M (noon)
    SMOOTHED_YAW_STEERING_1 = 7,              // e.g. GAL-1
    SMOOTHED_YAW_STEERING_2 = 8,              // e.g. GAL-2, BDS-3M, BDS-3I
    BETA_DEPENDENT_ORBIT_NORMAL = 9           // e.g. BDS-2M, BDS-2I, QZS-1
  };

  class AttitudeInfo
  {
  public:
    Time         timeStart;
    AttitudeMode defaultMode;
    AttitudeMode midnightMode;
    AttitudeMode noonMode;
    Double       maxYawRate;
    Double       yawBias;
    Double       midnightBetaThreshold;
    Double       noonBetaThreshold;
    Double       activationThreshold;
    Double       maxManeuverTime;

    AttitudeInfo(const Time &time, const Vector &data)
    {
      timeStart             = time;
      defaultMode           = AttitudeMode(data(0));
      midnightMode          = AttitudeMode(data(1));
      noonMode              = AttitudeMode(data(2));
      maxYawRate            = DEG2RAD*std::fabs(data(3));
      yawBias               = DEG2RAD*data(4);
      midnightBetaThreshold = DEG2RAD*std::fabs(data(5));
      noonBetaThreshold     = DEG2RAD*std::fabs(data(6));
      activationThreshold   = DEG2RAD*std::fabs(data(7));
      maxManeuverTime       = std::fabs(data(8));
    }
  };

  class Epoch
  {
  public:
    Time     time;
    Double   yawAngle;
    Double   yawRate;
    Double   orbitAngle;
    Double   betaAngle;
    Vector3d pos;
    Vector3d vel;
    Vector3d posSun;

    Epoch() : yawAngle(0), yawRate(0), orbitAngle(0), betaAngle(0) {}
  };

  EphemeridesPtr ephemerides;
  EclipsePtr     eclipse;
  std::vector<Epoch>        epochs;
  std::vector<AttitudeInfo> attitudeInfos;

  // helper methods
  Double       wrapAngle(Double angle) const;  ///< Returns angle wrapped to [-PI, PI).
  Rotary3d     orbitNormal2crf(const Vector3d &posSat, const Vector3d &velSat) const;
  AttitudeInfo getAttitudeInfo(const Time &time) const;
  Epoch        createDefaultEpoch(const Time &time, const Vector3d &posSat, const Vector3d &velSat) const;
  Bool         findShadowBoundaries(UInt idMidnightEpoch, Epoch &shadowStart, Epoch &shadowEnd) const;
  UInt         catchUpYawAngle(const Epoch startEpoch, Double maxYawRate, Bool backwards=FALSE);

  // attitude mode methods
  void modelNominalYawSteering(Epoch &epoch) const;
  void modelOrbitNormal(Epoch &epoch, const AttitudeInfo &attitudeInfo) const;
  UInt modelCatchUpYawSteering(UInt idEpoch, const AttitudeInfo &attitudeInfo);
  UInt modelShadowMaxYawSteeringAndRecovery(UInt idEpoch, const AttitudeInfo &attitudeInfo);
  UInt modelShadowMaxYawSteeringAndStop(UInt idEpoch, const AttitudeInfo &attitudeInfo);
  UInt modelShadowConstantYawSteering(UInt idEpoch);
  UInt modelCenteredMaxYawSteering(UInt idEpoch, const AttitudeInfo &attitudeInfo);
  UInt modelSmoothedYawSteering1(UInt idEpoch);
  UInt modelSmoothedYawSteering2(UInt idEpoch, const AttitudeInfo &attitudeInfo);
  UInt modelBetaDependentOrbitNormal(UInt idEpoch, const AttitudeInfo &attitudeInfo);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateStarCameraGnss, SINGLEPROCESS, "Simulates star camera data for a GNSS satellite based on an attitude model.", Simulation, Gnss, Instrument)

/***********************************************/

void SimulateStarCameraGnss::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameStarCamera, fileNameOrbit, fileNameAttitudeInfo;
    Double   modelingResolution;
    UInt     interpolationDegree;

    readConfig(config, "outputfileStarCamera",  fileNameStarCamera,   Config::MUSTSET,  "",  "rotation from body frame to CRF");
    readConfig(config, "inputfileOrbit",        fileNameOrbit,        Config::MUSTSET,  "",  "attitude is modeled based on this orbit");
    readConfig(config, "inputfileAttitudeInfo", fileNameAttitudeInfo, Config::OPTIONAL, "",  "attitude modes used by the satellite and respective parameters");
    readConfig(config, "interpolationDegree",   interpolationDegree,  Config::DEFAULT,  "7", "polynomial degree for orbit interpolation");
    readConfig(config, "modelingResolution",    modelingResolution,   Config::DEFAULT,  "1", "[s] resolution for attitude model evaluation");
    readConfig(config, "ephemerides",           ephemerides,          Config::MUSTSET,  "",  "");
    readConfig(config, "eclipse",               eclipse,              Config::MUSTSET,  "",  "model to determine if satellite is in Earth's shadow");
    if(isCreateSchema(config)) return;

    logStatus<<"read orbit file <"<<fileNameOrbit<<">"<<Log::endl;
    OrbitArc orbitArc = InstrumentFile::read(fileNameOrbit);
    if(orbitArc.size() && orbitArc.at(0).velocity.r()==0.)
      throw(Exception("orbit does not contain velocity data"));

    if(fileNameAttitudeInfo.empty())
    {
      logStatus<<"no attitude info provided, using nominal yaw-steering attitude"<<Log::endl;
      attitudeInfos.push_back(AttitudeInfo(Time(), Vector(9)));
    }
    else
    {
      logStatus<<"read attitude info from file <"<<fileNameAttitudeInfo<<">"<<Log::endl;
      MiscValuesArc arc = InstrumentFile::read(fileNameAttitudeInfo);
      for(auto epoch : arc)
        attitudeInfos.push_back(AttitudeInfo(epoch.time, epoch.data()));
    }

    // increase sampling and compute beta angle, orbit angle, and yaw angle/rate from default attitude mode for all epochs
    const std::vector<Time> timesOrbit = orbitArc.times();
    const Time deltaTime = seconds2time(std::fabs(modelingResolution));
    std::vector<Time> times;
    for(UInt i=0; timesOrbit.front()+i*deltaTime<=timesOrbit.back(); i++)
      times.push_back(timesOrbit.front() + i*deltaTime);
    Polynomial polynomial(interpolationDegree);
    {
      Matrix positionVelocity = polynomial.interpolate(times, timesOrbit, orbitArc.matrix().column(1, 6), 1);
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        epochs.push_back(createDefaultEpoch(times.at(idEpoch), Vector3d(positionVelocity.slice(idEpoch, 0, 1, 3)), Vector3d(positionVelocity.slice(idEpoch, 3, 1, 3))));
    }

    // Kepler-integrate forwards and backwards in time to consider maneuvers at orbit boundaries
    auto keplerExtend = [&](const Time &timeStart, const Vector3d &posStart, const Vector3d &velStart, Double integrationStep, Double integrationLimit)
    {
      Kepler kepler(timeStart, posStart, velStart);
      for(UInt i=1; (i-1)*std::fabs(integrationStep)<=std::fabs(integrationLimit); i++)
      {
        Time time = timeStart + seconds2time(i*integrationStep);
        Vector3d pos, vel;
        kepler.orbit(time, pos, vel);
        epochs.push_back(createDefaultEpoch(time, pos, vel));
        times.push_back(time);
      }
    };
    keplerExtend(orbitArc.front().time, orbitArc.front().position, orbitArc.front().velocity, -deltaTime.seconds(), getAttitudeInfo(orbitArc.front().time).maxManeuverTime);
    keplerExtend(orbitArc.back().time,  orbitArc.back().position,  orbitArc.back().velocity,   deltaTime.seconds(), getAttitudeInfo(orbitArc.back().time).maxManeuverTime);
    std::sort(epochs.begin(), epochs.end(), [](auto &e1, auto &e2){ return e1.time<e2.time; });
    std::sort(times.begin(), times.end());

    // model transitions caused by change of default attitude mode
    AttitudeMode previousDefaultMode = getAttitudeInfo(epochs.front().time).defaultMode;
    for(UInt idEpoch=1; idEpoch<epochs.size(); idEpoch++)
    {
      const AttitudeInfo attitudeInfo = getAttitudeInfo(epochs.at(idEpoch).time);
      if(attitudeInfo.defaultMode!=previousDefaultMode)
      {
        if(attitudeInfo.maxYawRate>0.)
        {
          epochs.at(idEpoch-1).yawRate = wrapAngle(epochs.at(idEpoch).yawAngle-epochs.at(idEpoch-1).yawAngle)>=0 ? attitudeInfo.maxYawRate : -attitudeInfo.maxYawRate;
          catchUpYawAngle(epochs.at(idEpoch-1), attitudeInfo.maxYawRate);
        }
        else
          logWarning<<"default attitude mode change at "<<epochs.at(idEpoch-1).time.dateTimeStr()<<" but unable to model transition because maxYawRate data is missing"<<Log::endl;
      }
      previousDefaultMode = attitudeInfo.defaultMode;
    }

    // apply specific attitude modes for noon and midnight
    for(UInt idEpoch=0; idEpoch<epochs.size(); idEpoch++)
    {
      const AttitudeInfo attitudeInfo = getAttitudeInfo(epochs.at(idEpoch).time);
      const AttitudeMode attitudeMode = std::fabs(epochs.at(idEpoch).orbitAngle)<=PI/2 ? attitudeInfo.midnightMode : attitudeInfo.noonMode;
      if(attitudeMode==attitudeInfo.defaultMode)
        continue; // same as default mode => nothing to adjust

      switch(attitudeMode)
      {
        case NOMINAL_YAW_STEERING:                 modelNominalYawSteering(epochs.at(idEpoch));                           break;
        case ORBIT_NORMAL:                         modelOrbitNormal(epochs.at(idEpoch), attitudeInfo);                    break;
        case CATCH_UP_YAW_STEERING:                idEpoch = modelCatchUpYawSteering(idEpoch, attitudeInfo);              break;
        case SHADOW_MAX_YAW_STEERING_AND_RECOVERY: idEpoch = modelShadowMaxYawSteeringAndRecovery(idEpoch, attitudeInfo); break;
        case SHADOW_MAX_YAW_STEERING_AND_STOP:     idEpoch = modelShadowMaxYawSteeringAndStop(idEpoch, attitudeInfo);     break;
        case SHADOW_CONSTANT_YAW_STEERING:         idEpoch = modelShadowConstantYawSteering(idEpoch);                     break;
        case CENTERED_MAX_YAW_STEERING:            idEpoch = modelCenteredMaxYawSteering(idEpoch, attitudeInfo);          break;
        case SMOOTHED_YAW_STEERING_1:              idEpoch = modelSmoothedYawSteering1(idEpoch);                          break;
        case SMOOTHED_YAW_STEERING_2:              idEpoch = modelSmoothedYawSteering2(idEpoch, attitudeInfo);            break;
        case BETA_DEPENDENT_ORBIT_NORMAL:          idEpoch = modelBetaDependentOrbitNormal(idEpoch, attitudeInfo);        break;
        default: throw(Exception(attitudeMode%"unknown attitude mode: %i"s));
      }
    }

    // Compute output star camera data
    Matrix quaternions(times.size(), 4);
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      copy(rotaryZ(-Angle(epochs.at(idEpoch).yawAngle)).quaternion().trans(), quaternions.row(idEpoch));
      if(idEpoch>0 && inner(quaternions.row(idEpoch), quaternions.row(idEpoch-1))<0.)
        quaternions.row(idEpoch) *= -1; // ensure same sign for correct interpolation
    }
    quaternions = polynomial.interpolate(timesOrbit, times, quaternions, 1);
    StarCameraArc starCameraArc;
    for(UInt idEpoch=0; idEpoch<timesOrbit.size(); idEpoch++)
    {
      StarCameraEpoch starCameraEpoch;
      starCameraEpoch.time   = timesOrbit.at(idEpoch);
      starCameraEpoch.rotary = orbitNormal2crf(orbitArc.at(idEpoch).position, orbitArc.at(idEpoch).velocity) * Rotary3d(quaternions.row(idEpoch).trans());
      starCameraArc.push_back(starCameraEpoch);
    }

    logStatus<<"write star camera file <"<<fileNameStarCamera<<">"<<Log::endl;
    InstrumentFile::write(fileNameStarCamera, starCameraArc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SimulateStarCameraGnss::wrapAngle(Double angle) const
{
  while(angle>= PI) angle -= 2*PI;
  while(angle< -PI) angle += 2*PI;

  return angle;
}

/***********************************************/

Rotary3d SimulateStarCameraGnss::orbitNormal2crf(const Vector3d &posSat, const Vector3d &velSat) const
{
  try
  {
    Vector3d r = normalize(posSat);
    Vector3d n = normalize(crossProduct(posSat, velSat));

    return Rotary3d(crossProduct(-n, -r), -n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SimulateStarCameraGnss::AttitudeInfo SimulateStarCameraGnss::getAttitudeInfo(const Time &time) const
{
  try
  {
    auto iter = std::find_if(attitudeInfos.rbegin(), attitudeInfos.rend(), [&](auto info){ return info.timeStart<=time; });
    if(iter==attitudeInfos.rend())
      throw(Exception("no attitude modes found for "+time.dateTimeStr()));
    return *iter;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SimulateStarCameraGnss::Epoch SimulateStarCameraGnss::createDefaultEpoch(const Time &time, const Vector3d &posSat, const Vector3d &velSat) const
{
  try
  {
    Epoch epoch;
    epoch.time   = time;
    epoch.pos    = posSat;
    epoch.vel    = velSat;
    epoch.posSun = ephemerides->position(epoch.time, Ephemerides::SUN);

    // Argument of latitude of satellite and sun
    const Vector3d z  = normalize(crossProduct(epoch.pos, epoch.vel));
    const Vector3d x  = normalize(crossProduct(Vector3d(0,0,1), z));
    const Vector3d y  = crossProduct(z, x);
    const Double   u  = std::atan2(inner(epoch.pos, y), inner(epoch.pos, x));
    const Double   u0 = std::atan2(inner(epoch.posSun, y), inner(epoch.posSun, x));

    epoch.orbitAngle = wrapAngle(u - u0 + PI);
    epoch.betaAngle  = std::acos(inner(-z, epoch.posSun)/epoch.posSun.r()) - PI/2;

    const AttitudeInfo attitudeInfo = getAttitudeInfo(time);
    if(attitudeInfo.defaultMode==NOMINAL_YAW_STEERING)
      modelNominalYawSteering(epoch);
    else if(attitudeInfo.defaultMode==ORBIT_NORMAL)
      modelOrbitNormal(epoch, attitudeInfo);
    else
      throw(Exception("only nominal yaw-steering mode or orbit normal mode are supported as default attitude modes"));

    return epoch;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool SimulateStarCameraGnss::findShadowBoundaries(UInt idMidnightEpoch, Epoch &shadowStart, Epoch &shadowEnd) const
{
  try
  {
    if(eclipse->factor(epochs.at(idMidnightEpoch).time, epochs.at(idMidnightEpoch).pos, ephemerides)>0.)
      return FALSE; // satellite never entered full shadow

    auto isOutsideShadow = [&](auto &epoch){ return eclipse->factor(epoch.time, epoch.pos, ephemerides)>0.; };

    auto iterEnd = std::find_if(epochs.begin()+idMidnightEpoch+1, epochs.end(), isOutsideShadow);
    if(iterEnd==epochs.end())
      return FALSE;

    auto iterStart = std::find_if(std::make_reverse_iterator(epochs.begin()+idMidnightEpoch), epochs.rend(), isOutsideShadow);
    if(iterStart==epochs.rend())
      return FALSE;

    shadowStart = *(iterStart-1);
    shadowEnd   = *(iterEnd-1);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt SimulateStarCameraGnss::catchUpYawAngle(const Epoch startEpoch, Double maxYawRate, Bool backwards)
{
  try
  {
    const Double startYawRate = (startEpoch.yawRate>=0 ? maxYawRate : -maxYawRate);
    Double previousYawAngleDiff = 0.;
    auto hasCaughtUpYaw = [&](UInt idEpoch)
    {
      Double yawAngle     = wrapAngle(startEpoch.yawAngle + startYawRate*(epochs.at(idEpoch).time-startEpoch.time).seconds());
      Double yawAngleDiff = wrapAngle(epochs.at(idEpoch).yawAngle - yawAngle);

      // Stop once true yaw angle catches up with nominal yaw angle
      if(((previousYawAngleDiff>0. && yawAngleDiff<0.) || (previousYawAngleDiff<0. && yawAngleDiff>0.)) && std::fabs(yawAngleDiff-previousYawAngleDiff)<PI)
        return TRUE;

      previousYawAngleDiff        = yawAngleDiff;
      epochs.at(idEpoch).yawAngle = yawAngle;
      epochs.at(idEpoch).yawRate  = startYawRate;

      return FALSE;
    };

    UInt idStartEpoch = std::distance(epochs.begin(), std::find_if(epochs.begin(), epochs.end(), [&](auto &epoch){ return epoch.time>startEpoch.time; }));
    if(backwards)
    {
      for(UInt idEpoch=idStartEpoch+1; idEpoch-->0;)
        if(hasCaughtUpYaw(idEpoch))
          return idEpoch;
    }
    else
    {
      for(UInt idEpoch=idStartEpoch; idEpoch<epochs.size(); idEpoch++)
        if(hasCaughtUpYaw(idEpoch))
          return idEpoch;
    }

    return backwards ? 0 : epochs.size();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SimulateStarCameraGnss::modelNominalYawSteering(Epoch &epoch) const
{
  try
  {
    epoch.yawAngle = wrapAngle(std::atan2(-std::tan(epoch.betaAngle), std::sin(epoch.orbitAngle)));
    epoch.yawRate  = epoch.vel.r()/epoch.pos.r() * std::tan(epoch.betaAngle)*std::cos(epoch.orbitAngle) /
                     (std::sin(epoch.orbitAngle)*std::sin(epoch.orbitAngle) + std::tan(epoch.betaAngle)*std::tan(epoch.betaAngle));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SimulateStarCameraGnss::modelOrbitNormal(Epoch &epoch, const AttitudeInfo &attitudeInfo) const
{
  try
  {
    epoch.yawAngle = wrapAngle(attitudeInfo.yawBias);
    epoch.yawRate  = 0.;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt SimulateStarCameraGnss::modelCatchUpYawSteering(UInt idEpoch, const AttitudeInfo &attitudeInfo)
{
  try
  {
    if(attitudeInfo.maxYawRate==0.)
      throw(Exception("required data missing for model: CATCH_UP_YAW_STEERING (maxYawRate)"));

    // check if maximum yaw rate is exceeded between current and next epoch
    if(idEpoch+1<epochs.size() && std::fabs(epochs.at(idEpoch).yawRate)<=attitudeInfo.maxYawRate && std::fabs(epochs.at(idEpoch+1).yawRate)>attitudeInfo.maxYawRate)
    {
      Epoch startEpoch = epochs.at(idEpoch+1);

      // consider anomalous turns at low beta angle for satellites with yaw bias
      if(attitudeInfo.yawBias!=0. && std::fabs(startEpoch.betaAngle)<std::fabs(attitudeInfo.yawBias) && attitudeInfo.yawBias*startEpoch.betaAngle>0.)
        startEpoch.yawRate = (startEpoch.betaAngle-attitudeInfo.yawBias)>=0. ? -attitudeInfo.maxYawRate : attitudeInfo.maxYawRate;

      return catchUpYawAngle(startEpoch, attitudeInfo.maxYawRate);
    }

    return idEpoch;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt SimulateStarCameraGnss::modelShadowMaxYawSteeringAndRecovery(UInt idEpoch, const AttitudeInfo &attitudeInfo)
{
  try
  {
    if(attitudeInfo.maxYawRate==0.)
      throw(Exception("required data missing for model: SHADOW_MAX_YAW_STEERING_AND_RECOVERY (maxYawRate)"));

    // check if satellite crosses midnight between current and next epoch
    if(std::fabs(epochs.at(idEpoch).orbitAngle)<=PI/2 && idEpoch+1<epochs.size() && epochs.at(idEpoch).orbitAngle<0. && epochs.at(idEpoch+1).orbitAngle>=0)
    {
      Epoch shadowStart, shadowEnd;
      if(!findShadowBoundaries(idEpoch, shadowStart, shadowEnd))
        return idEpoch; // cannot model maneuver without knowing shadow start and end

      // shadow maneuver
      Double startYawRate = 0.;
      if(attitudeInfo.yawBias>0.)
        startYawRate = attitudeInfo.maxYawRate;
      else if(attitudeInfo.yawBias<0.)
        startYawRate = -attitudeInfo.maxYawRate;
      else
        startYawRate = (shadowStart.yawRate>=0. ? attitudeInfo.maxYawRate : -attitudeInfo.maxYawRate);
      for(auto &&epoch : epochs)
        if(epoch.time>=shadowStart.time && epoch.time<=shadowEnd.time)
        {
          epoch.yawAngle = wrapAngle(shadowStart.yawAngle + startYawRate*(epoch.time-shadowStart.time).seconds());
          epoch.yawRate = startYawRate;
        }

      // post-shadow recovery maneuver
      const Double shadowEndNominalYawAngle = shadowEnd.yawAngle;
      shadowEnd.yawAngle = wrapAngle(shadowStart.yawAngle + startYawRate*(shadowEnd.time-shadowStart.time).seconds());
      shadowEnd.yawRate  = wrapAngle(shadowEnd.yawAngle-shadowEndNominalYawAngle)>=0 ? -attitudeInfo.maxYawRate : attitudeInfo.maxYawRate;
      return catchUpYawAngle(shadowEnd, attitudeInfo.maxYawRate);
    }

    return idEpoch;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt SimulateStarCameraGnss::modelShadowMaxYawSteeringAndStop(UInt idEpoch, const AttitudeInfo &attitudeInfo)
{
  try
  {
    if(attitudeInfo.maxYawRate==0.)
      throw(Exception("required data missing for model: SHADOW_MAX_YAW_STEERING_AND_STOP (maxYawRate)"));

    // check if satellite crosses midnight between current and next epoch
    if(std::fabs(epochs.at(idEpoch).orbitAngle)<=PI/2 && idEpoch+1<epochs.size() && epochs.at(idEpoch).orbitAngle<0. && epochs.at(idEpoch+1).orbitAngle>=0)
    {
      Epoch shadowStart, shadowEnd;
      if(!findShadowBoundaries(idEpoch, shadowStart, shadowEnd))
        return idEpoch; // cannot model maneuver without knowing shadow start and end

      // shadow maneuver
      for(UInt i=0; i<epochs.size(); i++)
        if(epochs.at(i).time>shadowStart.time && epochs.at(i).time<=shadowEnd.time)
        {
          epochs.at(i).yawAngle = shadowEnd.yawAngle;
          epochs.at(i).yawRate = 0.;
          idEpoch = i;
        }
      catchUpYawAngle(shadowStart, attitudeInfo.maxYawRate);
    }

    return idEpoch;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt SimulateStarCameraGnss::modelShadowConstantYawSteering(UInt idEpoch)
{
  try
  {
    // check if satellite crosses midnight between current and next epoch
    if(std::fabs(epochs.at(idEpoch).orbitAngle)<=PI/2 && idEpoch+1<epochs.size() && epochs.at(idEpoch).orbitAngle<0. && epochs.at(idEpoch+1).orbitAngle>=0)
    {
      Epoch shadowStart, shadowEnd;
      if(!findShadowBoundaries(idEpoch, shadowStart, shadowEnd))
        return idEpoch; // cannot model maneuver without knowing shadow start and end

      // shadow maneuver
      const Double startYawRate = wrapAngle(shadowEnd.yawAngle-shadowStart.yawAngle)/(shadowEnd.time-shadowStart.time).seconds();
      for(UInt i=0; i<epochs.size(); i++)
        if(epochs.at(i).time>=shadowStart.time && epochs.at(i).time<=shadowEnd.time)
        {
          epochs.at(i).yawAngle = wrapAngle(shadowStart.yawAngle + startYawRate*(epochs.at(i).time-shadowStart.time).seconds());
          epochs.at(i).yawRate  = startYawRate;
          idEpoch = i;
        }
    }

    return idEpoch;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt SimulateStarCameraGnss::modelCenteredMaxYawSteering(UInt idEpoch, const AttitudeInfo &attitudeInfo)
{
  try
  {
    if(attitudeInfo.maxYawRate==0.)
      throw(Exception("required data missing for model: CENTERED_MAX_YAW_STEERING (maxYawRate)"));

    if(std::fabs(epochs.at(idEpoch).orbitAngle)<=PI/2) // satellite is in midnight half of the orbit
    {
      Bool hasCrossedMidnight = idEpoch+1<epochs.size() && wrapAngle(epochs.at(idEpoch).orbitAngle)<0. &&  wrapAngle(epochs.at(idEpoch+1).orbitAngle)>=0;
      Bool isModeActive = attitudeInfo.midnightBetaThreshold==0. || std::abs(epochs.at(idEpoch).betaAngle)<attitudeInfo.midnightBetaThreshold;
      if(!hasCrossedMidnight || !isModeActive)
        return idEpoch;
    }
    else // satellite is in noon half of the orbit
    {
      Bool hasCrossedNoon = idEpoch+1<epochs.size() && wrapAngle(epochs.at(idEpoch).orbitAngle+PI)<0. &&  wrapAngle(epochs.at(idEpoch+1).orbitAngle+PI)>=0;
      Bool isModeActive = attitudeInfo.noonBetaThreshold==0. || std::abs(epochs.at(idEpoch).betaAngle)<attitudeInfo.noonBetaThreshold;
      if(!hasCrossedNoon || !isModeActive)
        return idEpoch;
    }

    catchUpYawAngle(epochs.at(idEpoch+1), attitudeInfo.maxYawRate, TRUE/*backwards*/);
    return catchUpYawAngle(epochs.at(idEpoch+1), attitudeInfo.maxYawRate);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt SimulateStarCameraGnss::modelSmoothedYawSteering1(UInt idEpoch)
{
  try
  {
    const Double sinBetaX = std::sin(15.0*DEG2RAD);
    const Double sinBetaY = std::sin(2.0*DEG2RAD);

    auto computeS0 = [&](const Epoch &epoch) { return orbitNormal2crf(epoch.pos, epoch.vel).inverseRotate(normalize(epoch.posSun-epoch.pos)); };
    auto isInAuxiliaryRegion = [&](const Vector3d &S0) { return std::fabs(S0.x())<sinBetaX && std::fabs(S0.y())<sinBetaY; };

    if(idEpoch+1<epochs.size() && isInAuxiliaryRegion(computeS0(epochs.at(idEpoch+1))) && !isInAuxiliaryRegion(computeS0(epochs.at(idEpoch))))
    {
      const Double Gamma = (computeS0(epochs.at(idEpoch+1)).y()<0 ? -1 : 1);
      for(UInt i=idEpoch+2; i<epochs.size(); i++)
      {
        const Vector3d S0 = computeS0(epochs.at(i));
        if(!isInAuxiliaryRegion(S0))
          return i;

        const Double SHy = 0.5*(sinBetaY*Gamma + S0.y()) + 0.5*(sinBetaY*Gamma - S0.y()) * std::cos(PI*std::fabs(S0.x())/sinBetaX);
        const Vector3d SH(S0.x(), SHy, std::sqrt(1. - S0.x()*S0.x() - SHy*SHy)*(S0.z()<0 ? -1 : 1));
        const Double denom = std::sqrt(1. - SH.z()*SH.z());
        epochs.at(i).yawAngle = std::atan2(SH.y()/denom, SH.x()/denom);
      }
    }

    return idEpoch;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt SimulateStarCameraGnss::modelSmoothedYawSteering2(UInt idEpoch, const AttitudeInfo &attitudeInfo)
{
  try
  {
    if(attitudeInfo.activationThreshold==0. || attitudeInfo.maxManeuverTime==0. ||
       (std::fabs(epochs.at(idEpoch).orbitAngle<=PI/2 ? attitudeInfo.midnightBetaThreshold : attitudeInfo.noonBetaThreshold))==0.)
      throw(Exception("required data missing for model: SMOOTHED_YAW_STEERING_2 (betaThreshold, activationThreshold, maxManeuverTime)"));

    auto isUsingModifiedLaw = [&](const Epoch &epoch)
    {
      Vector3d n = crossProduct(epoch.pos, epoch.vel);
      Double epsilon = std::acos(inner(normalize(epoch.pos), normalize(crossProduct(n, crossProduct(n, epoch.posSun)))));
      if(epsilon>PI/2)
        epsilon = PI - epsilon;

      if(std::fabs(epochs.at(idEpoch).orbitAngle)<=PI/2) // satellite is in midnight half of the orbit
        return std::fabs(epochs.at(idEpoch).betaAngle)<attitudeInfo.midnightBetaThreshold && std::fabs(epsilon)<attitudeInfo.activationThreshold;
      else // satellite is in noon half of the orbit
        return std::fabs(epochs.at(idEpoch).betaAngle)<attitudeInfo.noonBetaThreshold && std::fabs(epsilon)<attitudeInfo.activationThreshold;
    };

    if(idEpoch+1<epochs.size() && isUsingModifiedLaw(epochs.at(idEpoch+1)) && !isUsingModifiedLaw(epochs.at(idEpoch)))
      for(UInt i = idEpoch+2; i<epochs.size(); i++)
      {
        if(!isUsingModifiedLaw(epochs.at(i)))
          return idEpoch;

        const Double offset = (epochs.at(idEpoch+1).yawAngle<0 ? -PI/2 : PI/2);
        const Double tMod   = (epochs.at(i).time - epochs.at(idEpoch+1).time).seconds();
        epochs.at(i).yawAngle = offset + (epochs.at(idEpoch+1).yawAngle - offset) * std::cos(2*PI/attitudeInfo.maxManeuverTime * tMod);
      }

    return idEpoch;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt SimulateStarCameraGnss::modelBetaDependentOrbitNormal(UInt idEpoch, const AttitudeInfo &attitudeInfo)
{
  try
  {
    if(attitudeInfo.maxYawRate==0. ||attitudeInfo.activationThreshold==0. || attitudeInfo.maxManeuverTime==0. ||
       (std::fabs(epochs.at(idEpoch).orbitAngle<=PI/2 ? attitudeInfo.midnightBetaThreshold : attitudeInfo.noonBetaThreshold))==0.)
      throw(Exception("required data missing for model: BETA_DEPENDENT_ORBIT_NORMAL (maxYawRate, betaThreshold, activationThreshold, maxManeuverTime)"));

    auto isBelowBetaThreshold = [&](auto &epoch)
    {
      return std::fabs(epoch.betaAngle)<(std::fabs(epoch.orbitAngle)<=PI/2 ? attitudeInfo.midnightBetaThreshold : attitudeInfo.noonBetaThreshold);
    };
    auto isBelowYawThreshold = [&](auto &epoch) { return std::fabs(wrapAngle(epoch.yawAngle+attitudeInfo.yawBias))<=attitudeInfo.activationThreshold; };

    if(idEpoch+1<epochs.size() && isBelowBetaThreshold(epochs.at(idEpoch)) && !isBelowYawThreshold(epochs.at(idEpoch)) && isBelowYawThreshold(epochs.at(idEpoch+1)))
    {
      UInt idEpochEnd = idEpoch+1;
      for(UInt i=idEpoch+1; i<epochs.size(); i++)
      {
        if(!isBelowBetaThreshold(epochs.at(i)))
          break;
        if(isBelowYawThreshold(epochs.at(i)))
          idEpochEnd = i;
      }
      if((epochs.at(idEpochEnd).time-epochs.at(idEpoch+1).time).seconds()*attitudeInfo.maxYawRate<attitudeInfo.activationThreshold)
        return idEpoch; // period in orbit normal mode would be shorter than period required for transition

      const Double yawRateEntry = (wrapAngle(epochs.at(idEpoch).yawAngle-attitudeInfo.yawBias)>=0.)    ? -attitudeInfo.maxYawRate :  attitudeInfo.maxYawRate;
      const Double yawRateExit  = (wrapAngle(epochs.at(idEpochEnd).yawAngle-attitudeInfo.yawBias)>=0.) ?  attitudeInfo.maxYawRate : -attitudeInfo.maxYawRate;
      for(UInt i=idEpoch+1; i<=idEpochEnd; i++)
        modelOrbitNormal(epochs.at(i), attitudeInfo);
      epochs.at(idEpoch).yawRate = yawRateEntry;
      catchUpYawAngle(epochs.at(idEpoch), attitudeInfo.maxYawRate);
      epochs.at(idEpochEnd).yawRate = yawRateExit;
      return catchUpYawAngle(epochs.at(idEpochEnd), attitudeInfo.maxYawRate);
    }

    return idEpoch;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
