/***********************************************/
/**
* @file slr.cpp
*
* @brief global navigation satellite system.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr.h"
#include "slr/slrObservation.h"
#include "slr/slrDesignMatrix.h"
#include "slr/slrSatellite.h"
#include "slr/slrStation.h"
#include "slr/slrParametrization/slrParametrization.h"
#include "slr/slrSatelliteGenerator/slrSatelliteGenerator.h"
#include "slr/slrStationGenerator/slrStationGenerator.h"

/***********************************************/

void Slr::init(const std::vector<Time> &times, SlrSatelliteGeneratorPtr satelliteGenerator, SlrStationGeneratorPtr stationGenerator,
               EarthRotationPtr earthRotation, SlrParametrizationPtr parametrization)
{
  try
  {
    // init earth rotation
    // -------------------
    this->times = times;
    polynomialEop.init(times, 3);
    eop = Matrix(times.size(), 8); // Matrix eop columns: xp, yp, sp, deltaUT, LOD, X, Y, S
    for(UInt i=0; i<times.size(); i++)
      earthRotation->earthOrientationParameter(times.at(i), eop(i,0), eop(i,1), eop(i,2), eop(i,3), eop(i,4), eop(i,5), eop(i,6), eop(i,7));
    // UT1-UTC => UT1-GPS (avoid leap seconds jumps for interpolation)
    for(UInt i=0; i<times.size(); i++)
      eop(i,3) -= (times.at(i)-timeGPS2UTC(times.at(i))).seconds();

    funcRotationCrf2Trf = std::bind(&Slr::rotationCrf2Trf, this, std::placeholders::_1);

    // init satellites
    // -----------------
    satellites = satelliteGenerator->satellites(times);
    for(UInt idSat=0; idSat<satellites.size(); idSat++)
      satellites.at(idSat)->id_ = idSat;

    // init stations
    // --------------
    stations = stationGenerator->stations(times, satellites, earthRotation);
    for(UInt idStat=0; idStat<stations.size(); idStat++)
      stations.at(idStat)->id_ = idStat;

    // init parametrization
    // --------------------
    this->parametrization = parametrization;
    if(parametrization)
    {
      parametrization->init(this, parametrization->getParametrizationGravity());
      funcReduceModels = std::bind(&SlrParametrization::observationCorrections, parametrization, std::placeholders::_1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d Slr::rotationCrf2Trf(const Time &time) const
{
  try
  {
    Matrix eop = polynomialEop.interpolate({time}, this->eop);
    const Double xp      = eop(0, 0);
    const Double yp      = eop(0, 1);
    const Double sp      = eop(0, 2);
    const Double deltaUT = eop(0, 3) + (time-timeGPS2UTC(time)).seconds();
    const Double X       = eop(0, 5);
    const Double Y       = eop(0, 6);
    const Double S       = eop(0, 7);

    const Double ERA = Planets::ERA(timeGPS2UTC(time) + seconds2time(deltaUT));
    const Double r2  = X*X + Y*Y;
    const Double E   = (r2!=0.) ? std::atan2(Y, X) : 0.;
    const Double D   = std::atan(std::sqrt(r2/(1-r2)));

    return  rotaryX(Angle(-yp)) * rotaryY(Angle(-xp)) *
            rotaryZ(Angle(sp+ERA-S-E)) *
            rotaryY(Angle(D)) * rotaryZ(Angle(E));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Slr::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    logStatus<<"setup parameters"<<Log::endl;
    normalEquationInfo.initNewParameterNames();
    if(parametrization)
      parametrization->initParameter(normalEquationInfo);
    normalEquationInfo.calculateIndex();
    if(normalEquationInfo.parameterCount())
    {
      logInfo<<"+ ======="<<Log::endl;
      logInfo<<normalEquationInfo.parameterCount()%"%9i parameters in "s<<normalEquationInfo.blockCount()<<" normal equation matrix blocks"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector Slr::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo) const
{
  try
  {
    return parametrization ? parametrization->aprioriParameter(normalEquationInfo) : Vector();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Slr::basicObservationEquations(const SlrNormalEquationInfo &/*normalEquationInfo*/, UInt idStat, UInt idSat, UInt idPass, SlrObservationEquation &eqn) const
{
  try
  {
    eqn.compute(*stations.at(idStat)->observations.at(idSat).at(idPass),
                *stations.at(idStat), *satellites.at(idSat), funcRotationCrf2Trf, funcReduceModels, TRUE);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Slr::designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    if(eqn.l.rows() && parametrization)
      parametrization->designMatrix(normalEquationInfo, eqn, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Slr::constraints(const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(parametrization)
      parametrization->constraints(normalEquationInfo, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Slr::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz)
{
  try
  {
    return parametrization ? parametrization->updateParameter(normalEquationInfo, x, Wz) : 0.;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Slr::updateCovariance(const SlrNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
{
  try
  {
    if(parametrization)
      parametrization->updateCovariance(normalEquationInfo, covariance);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Slr::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix)
{
  try
  {
    if(parametrization)
      parametrization->writeResults(normalEquationInfo, suffix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Byte> Slr::selectStations(PlatformSelectorPtr selector)
{
  try
  {
    std::vector<const Platform*> platforms(stations.size(), nullptr);
    for(UInt idStat=0; idStat<stations.size(); idStat++)
      if(stations.at(idStat)->useable())
        platforms.at(idStat) = &stations.at(idStat)->platform;
    return selector->select(times.front(), times.back(), platforms);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Byte> Slr::selectSatellites(PlatformSelectorPtr selector)
{
  try
  {
    std::vector<const Platform*> platforms(satellites.size(), nullptr);
    for(UInt idSat=0; idSat<satellites.size(); idSat++)
      if(satellites.at(idSat)->useable())
        platforms.at(idSat) = &satellites.at(idSat)->platform;
    return selector->select(times.front(), times.back(), platforms);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Bool Slr::InfoParameterChange::update(Double change)
{
  try
  {
    count++;
    rms += change*change;
    if(std::fabs(change) > std::fabs(maxChange))
    {
      maxChange = change;
      return TRUE;
    }
    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Slr::InfoParameterChange::print(Double convertToMeter, Double &maxChangeTotal)
{
  try
  {
    rms = std::sqrt(rms/count);
    if(!info.empty())
    {
      std::string space(5-std::min(unit.size(), std::size_t(4)), ' ');
      logInfo<<"  rms ="<<rms%"%7.1f "s<<unit<<","<<space<<"max ="<<maxChange%"%8.1f "s<<unit<<","<<space<<info<<Log::endl;
    }

    maxChangeTotal = std::max(maxChangeTotal, std::fabs(convertToMeter*maxChange));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
