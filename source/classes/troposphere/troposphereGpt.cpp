/***********************************************/
/**
* @file troposphereGpt.cpp
*
* @brief GPT empirical troposphere model.
* @see Troposphere
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2017-02-23
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "troposphere.h"
#include "troposphereGpt.h"

/***********************************************/

TroposphereGpt::TroposphereGpt(Config &config)
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

void TroposphereGpt::init(const std::vector<Vector3d> &stationPositions)
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

Double TroposphereGpt::slantDelay(const Time &time, UInt stationId, Angle azimuth, Angle elevation) const
{
  try
  {
    computeEmpiricalCoefficients(time);

    const Double sinE = std::sin(elevation);
    const Double vmfh = mappingFunction(sinE, ah(stationId), bh(stationId), ch(stationId))
                      + (1./sinE - mappingFunction(sinE, a_ht, b_ht, c_ht)) * height(stationId)*0.001;
    const Double vmfw = mappingFunction(sinE, aw(stationId), bw(stationId), cw(stationId));
    const Double mfgh = 1. / (sinE*std::tan(elevation) + 0.0031); // hydrostatic gradient mapping function [Chen and Herring, 1997]
    const Double mfgw = 1. / (sinE*std::tan(elevation) + 0.0007); // wet -"-

    return vmfh*zhd(stationId) + vmfw*zwd(stationId)
        + (mfgh*gnh(stationId) + mfgw*gnw(stationId)) * std::cos(azimuth)
        + (mfgh*geh(stationId) + mfgw*gew(stationId)) * std::sin(azimuth);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TroposphereGpt::mappingFunctionHydrostatic(const Time &time, UInt stationId, Angle /*azimuth*/, Angle elevation) const
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

Double TroposphereGpt::mappingFunctionWet(const Time &time, UInt stationId, Angle /*azimuth*/, Angle elevation) const
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

void TroposphereGpt::mappingFunctionGradient(const Time &/*time*/, UInt /*stationId*/, Angle azimuth, Angle elevation, Double &dx, Double &dy) const
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

void TroposphereGpt::getAprioriValues(const Time &time, UInt stationId, Double &zenithDryDelay, Double &zenithWetDelay, Double &gradientDryNorth,
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
