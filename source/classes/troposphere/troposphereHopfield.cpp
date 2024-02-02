/***********************************************/
/**
* @file troposphereHopfield.cpp
*
* @brief Hopfield troposphere model.
* @see Troposphere
*
* @author André Hauschild
* @date 2023-04-27
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "troposphere.h"
#include "troposphereHopfield.h"

/***********************************************/

TroposphereHopfield::TroposphereHopfield(Config &config)
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

void TroposphereHopfield::init(const std::vector<Vector3d> &stationPositions)
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

Double TroposphereHopfield::slantDelay(const Time &time, UInt stationId, Angle azimuth, Angle elevation) const
{
  try
  {
    /*computeEmpiricalCoefficients(time);*/
    const Double nmfh = mappingFunctionHydrostatic(time,stationId,azimuth,elevation);
    const Double nmfw = mappingFunctionWet        (time,stationId,azimuth,elevation);

    const Double T =  291.1;
    const Double P = 1010.0;
    const Double e =   10.4;

    // Dry and wet zenith delay [m]

    /* EPOS implementation */
    const Double zhd = (77.6e-6 * (-613.3768/T + 148.98) * P)/5.0;
    const Double zwd = (77.6e-6 * 11000.0 * 4810.0 * e/pow(T,2))/5.0;

    /*
     * according to Hoffman-Wellenhof 2008, p. 132
    const Double zhd = 77.64e-6/5.0*met.P/met.T*(40136+148.72*(met.T-273.16));
    const Double zwd = 11000e-6/5.0*met.e/pow(met.T,2)*(-12.96*met.T+3.718e5);
    */

    return nmfh*zhd + nmfw*zwd;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TroposphereHopfield::mappingFunctionHydrostatic(const Time &/*time*/, UInt /*stationId*/, Angle /*azimuth*/, Angle elevation) const
{
  try
  {
    /*computeEmpiricalCoefficients(time);*/
    return 1.0/sin(sqrt(pow(elevation,2)+pow(2.5*DEG2RAD,2)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TroposphereHopfield::mappingFunctionWet(const Time &/*time*/, UInt /*stationId*/, Angle /*azimuth*/, Angle elevation) const
{
  try
  {
    /*computeEmpiricalCoefficients(time);*/
    return 1.0/sin(sqrt(pow(elevation,2)+pow(1.5*DEG2RAD,2)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TroposphereHopfield::mappingFunctionGradient(const Time &/*time*/, UInt /*stationId*/, Angle azimuth, Angle elevation, Double &dx, Double &dy) const
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

void TroposphereHopfield::getAprioriValues(const Time &time, UInt stationId, Double &zenithDryDelay, Double &zenithWetDelay, Double &gradientDryNorth,
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
