/***********************************************/
/**
* @file slrObservation.cpp
*
* @brief Code & Phase observations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "slr/slrSatellite.h"
#include "slr/slrStation.h"
#include "slr/slrObservation.h"

/***********************************************/

static Bool positionVelocityTime(const SlrStation &station, const SlrSatellite &satellite,
                                 const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, const Time &timeTrans,
                                 Vector3d &posTrans, Vector3d &velTrans, Time &timeBounce, Vector3d &posBounce, Vector3d &velBounce,
                                 Time &timeRecv, Vector3d &posRecv)
{
  try
  {
    // station position at transmit time
    posTrans = rotationCrf2Trf(timeTrans).inverseRotate(station.position(timeTrans));
    if(std::isnan(posTrans.quadsum()))
      return FALSE;

    // satellite position at bounce time
    posBounce = satellite.position(timeTrans);
    if(std::isnan(posBounce.quadsum()))
      return FALSE;
    Vector3d posOld;
    for(UInt i=0; (i<10) && ((posBounce-posOld).r() > 0.0001); i++) // iteration
    {
      timeBounce = timeTrans + seconds2time((posBounce-posTrans).r()/LIGHT_VELOCITY);
      posOld     = posBounce;
      posBounce  = satellite.position(timeBounce) + satellite.centerOfMass2Reflector(timeBounce, normalize(posBounce-posTrans));
    }
    velBounce = satellite.velocity(timeBounce);

    // station position at receive time
    timeRecv = timeBounce + seconds2time((posBounce-posTrans).r()/LIGHT_VELOCITY);
    posRecv  = rotationCrf2Trf(timeRecv).inverseRotate(station.position(timeRecv));
    for(UInt i=0; (i<10) && ((posRecv-posOld).r() > 0.0001); i++) // iteration
    {
      timeRecv = timeBounce + seconds2time((posRecv-posBounce).r()/LIGHT_VELOCITY);
      posOld   = posRecv;
      posRecv  = rotationCrf2Trf(timeRecv).inverseRotate(station.position(timeRecv));
    }
    velTrans = 1./(timeRecv-timeTrans).seconds() * (posRecv-posTrans);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Bool SlrObservation::init(const SlrStation &station, const SlrSatellite &satellite,
                          const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, Angle elevationCutOff)
{
  try
  {
    if((!station.useable()) || (!satellite.useable()))
      return FALSE;

    const UInt obsCount = observations.rows();
    Vector useable(obsCount);
    Vector sigmas0Tmp(obsCount);
    std::vector<Time> timesTransCorrected = station.correctedTimes(timesTrans);

    for(UInt i=0; i<obsCount; i++)
    {
      Time timeBounce, timeRecv;
      Vector3d posTrans, posBounce, posRecv;
      Vector3d velTrans, velBounce;
      if(!positionVelocityTime(station, satellite, rotationCrf2Trf, timesTransCorrected.at(i),
                               posTrans, velTrans, timeBounce, posBounce, velBounce, timeRecv, posRecv))
        continue;

      // direction to satellite in local frame (north, east, up)
      const Vector3d k         = station.global2localFrame(timesTransCorrected.at(i)).transform(rotationCrf2Trf(timesTransCorrected.at(i)).rotate(normalize(posBounce-posTrans)));
      const Angle    azimut    = k.lambda();
      const Angle    elevation = k.phi();

      sigmas0Tmp(i) = station.accuracy(timesTransCorrected.at(i),
                                       observations(i)-0.5*((posBounce-posTrans).r() + (posRecv-posBounce).r()),
                                       sigmas0(i), redundancies(i), laserWavelength(i), azimut, elevation);
      useable(i) = (k.phi() >= elevationCutOff) && (sigmas0Tmp(i) > 0) && !std::isnan(posRecv.quadsum());
    }

    // Remove unused observations
    const UInt obsCountNew = sum(useable);
    auto timesTmp      = timesTrans;
    auto obsTmp        = observations;
    auto wavelengthTmp = laserWavelength;
    timesTrans.resize(obsCountNew);
    observations       = Vector(obsCountNew);
    residuals          = Vector(obsCountNew);
    redundancies       = Vector(obsCountNew);
    sigmas0            = Vector(obsCountNew);
    sigmas             = Vector(obsCountNew);
    laserWavelength    = Vector(obsCountNew);
    UInt idx = 0;
    for(UInt i=0; i<useable.size(); i++)
      if(useable(i))
      {
        timesTrans.at(idx)   = timesTmp.at(i);
        observations(idx)    = obsTmp(i);
        sigmas(idx) = sigmas0(idx) = sigmas0Tmp(i);
        laserWavelength(idx) = wavelengthTmp(i);
        idx++;
      }

    return (obsCountNew > 0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrObservation::setHomogenizedResiduals(const_MatrixSliceRef residuals_, const_MatrixSliceRef redundancies_)
{
  try
  {
    residuals    = residuals_;
    redundancies = redundancies_;
    for(UInt i=0; i<residuals.rows(); i++)
      residuals(i) *= sigmas(i);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void SlrObservationEquation::compute(const SlrObservation &observation, const SlrStation &station_, const SlrSatellite &satellite_,
                                     const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                                     const std::function<void(SlrObservationEquation &eqn)> &reduceModels, Bool homogenize)
{
  try
  {

    type        = RANGE;
    station     = &station_;
    satellite   = &satellite_;

    const UInt obsCount = observation.observations.rows();
    sigmas      = observation.sigmas;
    sigmas0     = observation.sigmas0;
    l           = observation.observations;
    A           = Matrix(obsCount, 8);
    index.resize(obsCount);
    std::iota(index.begin(), index.end(), 0);
    count       = std::vector<UInt>(obsCount, 1);
    timesStat   = station->correctedTimes(observation.timesTrans);
    timesSat.resize(obsCount);
    posStat.resize(obsCount);
    posSat.resize(obsCount);
    azimutStat.resize(obsCount);
    elevationStat.resize(obsCount);
    laserWavelength = observation.laserWavelength;

    for(UInt i=0; i<obsCount; i++)
    {
      Time     timesRecv;
      Vector3d posRecv, velStat, velSat;
      positionVelocityTime(station_, satellite_, rotationCrf2Trf, timesStat.at(i),
                           posStat.at(i), velStat, timesSat.at(i), posSat.at(i), velSat, timesRecv, posRecv);

      // reduced observations
      l(i) -= 0.5*((posSat.at(i)-posStat.at(i)).r() + (posRecv-posSat.at(i)).r());

      // curved space-time
      const Double r1  = posStat.at(i).r();
      const Double r2  = posSat.at(i).r();
      const Double r12 = (posStat.at(i) - posSat.at(i)).r();
      l(i) += 2*DEFAULT_GM/std::pow(LIGHT_VELOCITY,2)*log((r1+r2+r12)/(r1+r2-r12));

      // geometry of roundtrip
      // ---------------------
      const Vector3d k1     = normalize(posSat.at(i) - posStat.at(i)); // unit vector stat -> sat
      const Vector3d k2     = normalize(posRecv - posSat.at(i));       // unit vector sat -> stat
      const Double   d1Stat = inner(k1, velStat)/LIGHT_VELOCITY;
      const Double   d1Sat  = inner(k1, velSat) /LIGHT_VELOCITY;
      const Double   d2Stat = inner(k2, velStat)/LIGHT_VELOCITY;
      const Double   d2Sat  = inner(k2, velSat) /LIGHT_VELOCITY;

      // direction to satellite in local frame (north, east, up)
      const Vector3d kLocal = station_.global2localFrame(timesStat.at(i)).transform(rotationCrf2Trf(timesStat.at(i)).rotate(k1));
      azimutStat.at(i)    = kLocal.lambda();
      elevationStat.at(i) = kLocal.phi();

      // design matrix
      // -------------
      A(i, idxPosStat+0) = 0.5*(k1.x() * (-1+d1Sat) + k2.x() * (+1+d1Sat+d1Stat+d2Stat));  // station coord x
      A(i, idxPosStat+1) = 0.5*(k1.y() * (-1+d1Sat) + k2.y() * (+1+d1Sat+d1Stat+d2Stat));  // station coord y
      A(i, idxPosStat+2) = 0.5*(k1.z() * (-1+d1Sat) + k2.z() * (+1+d1Sat+d1Stat+d2Stat));  // station coord z
      A(i, idxPosSat+0)  = 0.5*(k1.x() * (+1+d1Sat) + k2.x() * (-1+d1Sat+d1Stat+d2Stat));  // satellite coord x
      A(i, idxPosSat+1)  = 0.5*(k1.y() * (+1+d1Sat) + k2.y() * (-1+d1Sat+d1Stat+d2Stat));  // satellite coord y
      A(i, idxPosSat+2)  = 0.5*(k1.z() * (+1+d1Sat) + k2.z() * (-1+d1Sat+d1Stat+d2Stat));  // satellite coord z
      A(i, idxTime)      = 0.5*LIGHT_VELOCITY*(d1Sat-d1Stat+d2Stat-d2Sat);                 // with respect to transmit time
      A(i, idxRange)     = 1.0;                                                            // one way range correction (troposphere, ...)
    }

    if(reduceModels)
      reduceModels(*this);

    // homogenize
    // -----------
    if(homogenize)
      for(UInt i=0; i<obsCount; i++)
      {
        if(l.size()) l.row(i) *= 1/sigmas(i);
        if(A.size()) A.row(i) *= 1/sigmas(i);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
