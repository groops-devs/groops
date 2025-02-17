/***********************************************/
/**
* @file troposphereSaastamoinen.cpp
*
* @brief Saastamoinen troposphere model.
* @see Troposphere
*
* @author André Hauschild
* @date 2023-04-27
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "troposphere.h"
#include "troposphereSaastamoinen.h"

/***********************************************/

TroposphereSaastamoinen::TroposphereSaastamoinen(Config &config)
{
  try
  {
    readConfig(config, "inputfileGpt", fileNameGpt, Config::MUSTSET, "{groopsDataDir}/troposphere/gpt3_grid1deg.dat",  "gridded GPT data");
    readConfig(config, "aHeight",      a_ht,        Config::DEFAULT, "2.53e-5", "parameter a (height correction)");
    readConfig(config, "bHeight",      b_ht,        Config::DEFAULT, "5.49e-3", "parameter b (height correction)");
    readConfig(config, "cHeight",      c_ht,        Config::DEFAULT, "1.14e-3", "parameter c (height correction)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TroposphereSaastamoinen::init(const std::vector<std::string> &/*stationNames*/, const std::vector<Vector3d> &stationPositions)
{
  try
  {
    initEmpiricalCoefficients(fileNameGpt, stationPositions);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TroposphereSaastamoinen::slantDelay(UInt stationId, const Time &time, Double frequency, Angle azimuth, Angle elevation) const
{
  try
  {
    computeEmpiricalCoefficients(time);
    const Double nmfh = mappingFunctionHydrostatic(stationId,time,frequency,azimuth,elevation);
    const Double nmfw = mappingFunctionWet        (stationId,time,frequency,azimuth,elevation);

    return nmfh*zhd(stationId) + nmfw*zwd(stationId);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TroposphereSaastamoinen::mappingFunctionHydrostatic(UInt stationId, const Time &time, Double /*frequency*/, Angle /*azimuth*/, Angle elevation) const
{
  try
  {
    computeEmpiricalCoefficients(time);
    const Double sinE  = sin(elevation);
    return mappingFunction(sinE, ah(stationId), bh(stationId), ch(stationId))
           + (1./sinE - mappingFunction(sinE, a_ht, b_ht, c_ht)) * height(stationId)*0.001;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TroposphereSaastamoinen::mappingFunctionWet(UInt stationId, const Time &time, Double /*frequency*/, Angle /*azimuth*/, Angle elevation) const
{
  try
  {
    computeEmpiricalCoefficients(time);
    return mappingFunction(sin(elevation), aw(stationId), bw(stationId), cw(stationId));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TroposphereSaastamoinen::mappingFunctionGradient(UInt /*stationId*/, const Time &/*time*/, Double /*frequency*/, Angle azimuth, Angle elevation, Double &dx, Double &dy) const
{
  try
  {
    const Double mfgw = 1./(std::sin(elevation)*std::tan(elevation) + 0.0031); // hydrostatic gradient mapping function [Chen and Herring, 1997] (unitless)
    dx = mfgw * std::cos(azimuth);
    dy = mfgw * std::sin(azimuth);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TroposphereSaastamoinen::getAprioriValues(UInt stationId, const Time &time, Double /*frequency*/, Double &zenithDryDelay, Double &zenithWetDelay, Double &gradientDryNorth,
                                               Double &gradientWetNorth, Double &gradientDryEast, Double &gradientWetEast, Double &aDry, Double &aWet) const
{
  try
  {
    computeEmpiricalCoefficients(time);
    zenithDryDelay   = zhd(stationId);
    zenithWetDelay   = zwd(stationId);
    gradientDryNorth = gnh(stationId);
    gradientWetNorth = gnw(stationId);
    gradientDryEast  = geh(stationId);
    gradientWetEast  = gew(stationId);
    aDry             = ah(stationId);
    aWet             = aw(stationId);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
