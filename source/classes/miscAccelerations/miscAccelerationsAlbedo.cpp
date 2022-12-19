/***********************************************/
/**
* @file miscAccelerationsAlbedo.cpp
*
* @brief DEPRECATED. Use radiationPressure instead.
* @see MiscAccelerations
*
* @author Torsten Mayer-Guerr
* @date 2014-21-10
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "files/fileGriddedData.h"
#include "files/fileSatelliteModel.h"
#include "classes/miscAccelerations/miscAccelerations.h"
#include "classes/miscAccelerations/miscAccelerationsRadiationPressure.h"
#include "classes/miscAccelerations/miscAccelerationsAlbedo.h"

/***********************************************/

MiscAccelerationsAlbedo::MiscAccelerationsAlbedo(Config &config)
{
  try
  {
    FileName    fileNameReflectivity, fileNameEmissivity;

    renameDeprecatedConfig(config, "inputfileEmissity", "inputfileEmissivity", date2time(2020, 8, 20));

    readConfig(config, "inputfileReflectivity", fileNameReflectivity, Config::OPTIONAL, "{groopsDataDir}/albedo/earth_ceres_reflectivity_grid2.5.dat", "");
    readConfig(config, "inputfileEmissivity",   fileNameEmissivity,   Config::OPTIONAL, "{groopsDataDir}/albedo/earth_ceres_emissivity_grid2.5.dat", "");
    readConfig(config, "solarflux",             solarflux,            Config::DEFAULT,  "1367", "solar flux constant in 1 AU [W/m**2]");
    readConfig(config, "factor",                factor,               Config::DEFAULT,  "1.0",  "the result is multiplied by this factor, set -1 to subtract the field");
    if(isCreateSchema(config)) return;

    // read reflectivity
    if(!fileNameReflectivity.empty())
    {
      GriddedData grid;
      readFileGriddedData(fileNameReflectivity, grid);
      points       = grid.points;
      areas        = grid.areas;
      reflectivity = grid.values;
    }

    // read emissivity
    if(!fileNameEmissivity.empty())
    {
      GriddedData grid;
      readFileGriddedData(fileNameEmissivity, grid);
      if(points.size() && (points.size() != grid.points.size()))
        throw(Exception("reflectivity and emissivity must be given on same grid"));
      points   = grid.points;
      areas    = grid.areas;
      emissivity = grid.values;
    }

    // convert area from unit sphere
    for(UInt i=0; i<points.size(); i++)
      areas.at(i) *= pow(points.at(i).r(), 2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d MiscAccelerationsAlbedo::acceleration(SatelliteModelPtr satellite, const Time &time,
                                               const Vector3d &position, const Vector3d &/*velocity*/,
                                               const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides)
{
  try
  {
    if(!satellite)
      throw(Exception("No satellite model given"));
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    Vector3d posSat = rotSat.inverseRotate(position);
    Vector3d posSun = rotSat.inverseRotate(ephemerides->position(time, Ephemerides::SUN));

    const Double AU           = 149597870700.0;
    const Double distanceSun  = posSun.normalize();
    const Double s0           = solarflux/LIGHT_VELOCITY*pow(AU/distanceSun,2);

    // gridded data given monthly?
    UInt index1 = 0;
    UInt index2 = 0;
    if((reflectivity.size()==12) || (emissivity.size()==12))
    {
      UInt   year, month, day, hour, minute;
      Double second;
      time.date(year, month, day, hour, minute, second);
      index1 = std::min(month-1, reflectivity.size()-1);
      index2 = std::min(month-1, emissivity.size()-1);
    }

    Vector3d F; // force in SRF
    Vector   absorbedPressure(satellite->surfaces.size()); // power for each satellite surface
    for(UInt i=0; i<points.size(); i++)
    {
      Vector3d posEarth     = rotSat.inverseRotate(rotEarth.inverseRotate(points.at(i)));
      Vector3d direction    = posSat-posEarth;
      const Double distance = direction.normalize();
      posEarth.normalize();

      // Cosine of angle of reflected radiation
      const Double cosReflexion = inner(direction, posEarth);
      if(cosReflexion<=0) // element visible?
        continue;
      // Cosine of angle of incident radiation
      const Double cosIncident = inner(posSun, posEarth);

      // Reflected Irradiance
      Double eReflectivity = 0;
      if(reflectivity.size() && (cosIncident>0))
        eReflectivity = reflectivity.at(index1).at(i)/(PI*distance*distance)*cosIncident*s0*cosReflexion*areas.at(i);

      // Emitted Irradiance
      Double eEmittance = 0;
      if(emissivity.size())
        eEmittance = emissivity.at(index2).at(i)/(4*PI*distance*distance)*s0*cosReflexion*areas.at(i);

      F += MiscAccelerationsRadiationPressure::force(satellite, direction, eReflectivity, eEmittance, absorbedPressure);
    }

    // compute thermal radition
    for(UInt i=0; i<satellite->surfaces.size(); i++)
    {
      auto &surface = satellite->surfaces.at(i);
      if(surface.specificHeatCapacity < 0)
        F -= (2./3.) * surface.area * absorbedPressure(i) * surface.normal; // spontaneous reemission of absorbed radiation
    }

    return (factor/satellite->mass) * rotEarth.rotate(rotSat.rotate(F));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
