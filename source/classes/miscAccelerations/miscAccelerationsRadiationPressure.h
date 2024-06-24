/***********************************************/
/**
* @file miscAccelerationsRadiationPressure.h
*
* @brief Solar and Earth radiation pressure model.
* @see MiscAccelerations
*
* @author Torsten Mayer-Guerr
* @date 2022-08-23
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONSRADIATIONPRESSURE__
#define __GROOPS_MISCACCELERATIONSRADIATIONPRESSURE__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerationsRadiationPressure = R"(
\subsection{RadiationPressure}\label{miscAccelerationsType:RadiationPressure}
This class computes acceleration acting on a satellite caused by Solar and Earth radiation pressure
and thermal radiation.

Solar radiation pressure: The solar constant at 1~AU can be set via \config{solarFlux}.
The \config{factorSolarRadation} can be used to scale the computed acceleration of the direct solar radiation.

Earth radiation pressure:
Input are a time series of gridded albedo values (unitless) as \configFile{inputfileAlbedoTimeSeries}{griddedDataTimeSeries}
and a time series of gridded longwave flux (W/m$^2$) as \configFile{inputfileLongwaveFluxTimeSeries}{griddedDataTimeSeries}.
Both files are optional and if not specified, the respective effect on the acceleration is not computed.
The \config{factorEarthRadation} can be used to scale the computed acceleration of the earth radiation.

The thermal radiation (TRP) of the satellite itself is either computed as direct re-emission or
based on the actual temperature of the satellite surfaces, depending on the seetings of the
\file{satellite macro model}{satelliteModel}. The second one uses a transient temperature model
with a temporal differential equation which disallows parallel computing.
The \config{factorThermalRadiation} can be used to scale the computed acceleration of the TRP.

The algorithms are described in:

Woeske et. al. (2019), GRACE accelerometer calibration by high precision non-gravitational force modeling,
Advances in Space Research, \url{https://doi.org/10.1016/j.asr.2018.10.025}.
)";
#endif

/***********************************************/

#include "files/fileGriddedDataTimeSeries.h"
#include "classes/eclipse/eclipse.h"
#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief Solar and Earth radiation pressure model.
 * @ingroup miscAccelerationsGroup
 * @see MiscAccelerations */
class MiscAccelerationsRadiationPressure : public MiscAccelerationsBase
{
  class State
  {
  public:
    SatelliteModelPtr satellite;
    Time              time;
    Vector            absorbedPressure;  // for each surface
    Bool              needTemperature;
    Vector            temperature;       // for each surface
  };

  Double                      solarflux;
  EclipsePtr                  eclipse;
  Bool                        hasAlbedo, hasFlux;
  InFileGriddedDataTimeSeries fileAlbedoTimeSeries, fileFluxTimeSeries;
  std::vector<Vector3d>       points;
  std::vector<Double>         areas;
  Double                      factorSun, factorEarth, factorThermal;
  std::vector<State>          states;

public:
  MiscAccelerationsRadiationPressure(Config &config);

  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) override;

  static Vector3d force(SatelliteModelPtr satellite, const Vector3d &direction, Double visible, Double infrared, Vector &absorbedPressure);
};

/***********************************************/

inline MiscAccelerationsRadiationPressure::MiscAccelerationsRadiationPressure(Config &config)
{
  try
  {
    FileName fileNameInAlbedo, fileNameInFlux;

    readConfig(config, "solarflux",                       solarflux,        Config::DEFAULT,  "1367", "solar flux constant in 1 AU [W/m^2]");
    readConfig(config, "eclipse",                         eclipse,          Config::MUSTSET,  "",     "");
    readConfig(config, "inputfileAlbedoTimeSeries",       fileNameInAlbedo, Config::OPTIONAL, "{groopsDataDir}/earthRadiationPressure/CERES_SYN1deg_albedo_climatology.dat",       "GriddedDataTimeSeries of albedo values (unitless)");
    readConfig(config, "inputfileLongwaveFluxTimeSeries", fileNameInFlux,   Config::OPTIONAL, "{groopsDataDir}/earthRadiationPressure/CERES_SYN1deg_longwaveFlux_climatology.dat", "GriddedDataTimeSeries of longwave flux values [W/m^2]");
    readConfig(config, "factorSolarRadation",             factorSun,        Config::DEFAULT,  "1.0",  "Solar radiation pressure is multiplied by this factor");
    readConfig(config, "factorEarthRadation",             factorEarth,      Config::DEFAULT,  "1.0",  "Earth radiation preussure is multiplied by this factor");
    readConfig(config, "factorThermalRadiation",          factorThermal,    Config::DEFAULT,  "1.0",  "Thermal (re-)radiation is multiplied by this factor");
    if(isCreateSchema(config)) return;

    hasAlbedo = FALSE;
    if(!fileNameInAlbedo.empty())
    {
      fileAlbedoTimeSeries.open(fileNameInAlbedo);
      GriddedData grid = fileAlbedoTimeSeries.grid();
      points    = grid.points;
      areas     = grid.areas;
      hasAlbedo = TRUE;
    }

    hasFlux = FALSE;
    if(!fileNameInFlux.empty())
    {
      fileFluxTimeSeries.open(fileNameInFlux);
      GriddedData grid = fileFluxTimeSeries.grid();
      if(points.size() && (points.size() != grid.points.size()))  // grid is already set up from albedo file
        throw(Exception("Albedo and longwave flux must be given on same grid"));
      points  = grid.points;
      areas   = grid.areas;
      hasFlux = TRUE;
    }

    // convert area from unit sphere
    for(UInt i=0; i<points.size(); i++)
      areas.at(i) *= std::pow(points.at(i).r(), 2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsRadiationPressure::acceleration(SatelliteModelPtr satellite, const Time &time,
                                                                 const Vector3d &position, const Vector3d &/*velocity*/,
                                                                 const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides)
{
  try
  {
    if(!satellite)
      throw(Exception("No satellite model given"));
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    Vector absorbedPressure(satellite->surfaces.size()); // power for each satellite surface

    // from sun to satellite in SRF
    Vector3d posSat = rotSat.inverseRotate(position);
    Vector3d posSun = rotSat.inverseRotate(ephemerides->position(time, Ephemerides::SUN));
    constexpr Double AU = 149597870700.0;
    const Double ny = eclipse ? eclipse->factor(time, position, ephemerides) : 1; // shadowing from Earth and Moon

    // Solar Radiation Pressure (SRP)
    // ------------------------------
    Vector3d     directionSunSat = posSat-posSun;
    const Double distanceSunSat  = directionSunSat.normalize();
    Vector3d F = factorSun * force(satellite, directionSunSat, ny*solarflux/LIGHT_VELOCITY*std::pow(AU/distanceSunSat, 2), 0., absorbedPressure);

    // Earth Radiation Pressure (ERP)
    // ------------------------------
    if((factorEarth || factorThermal) && (hasAlbedo || hasFlux))
    {
      Vector dataAlbedo, dataFlux;
      if(hasAlbedo) dataAlbedo = fileAlbedoTimeSeries.data(time);
      if(hasFlux)   dataFlux   = fileFluxTimeSeries.data(time);

      for(UInt i=0; i<points.size(); i++)
      {
        Vector3d posEarth             = rotSat.inverseRotate(rotEarth.inverseRotate(points.at(i))); // in SRF
        Vector3d directionEarthSat    = posSat-posEarth;
        const Double distanceEarthSat = directionEarthSat.normalize();

        // Cosine of angle of reflected radiation
        const Double cosReflexion = inner(directionEarthSat, normalize(posEarth));
        if(cosReflexion <= 0) // element visible?
          continue;
        const Double f = cosReflexion*areas.at(i)/(PI*distanceEarthSat*distanceEarthSat);

        // Reflected Irradiance
        Double eReflectivity = 0;
        if(hasAlbedo)
        {
          Vector3d directionSunEarth    = posEarth-posSun;
          const Double distanceSunEarth = directionSunEarth.normalize();
          const Double cosIncident      = -inner(directionSunEarth, normalize(posEarth)); // Cosine of angle of incident radiation
          if(cosIncident > 0)
            eReflectivity = f * dataAlbedo(i) * cosIncident * solarflux/LIGHT_VELOCITY * std::pow(AU/distanceSunEarth, 2);
        }

        // Emitted Irradiance
        Double eEmittance = 0;
        if(hasFlux)
          eEmittance = f * dataFlux(i)/LIGHT_VELOCITY;

        F += factorEarth * force(satellite, directionEarthSat, eReflectivity, eEmittance, absorbedPressure);
      }
    } // if((factorEarth || factorThermal) && (hasAlbedo || hasFlux)

    // Thermal Radiation Pressure (TRP)
    // --------------------------------
    if(factorThermal && absorbedPressure.size() && quadsum(absorbedPressure))
    {
      constexpr Double sigma = 5.670374419e-8; // Stefan–Boltzmann constant [W/m^2/K^4]

      // set temperature to be in equillibrium with input power
      auto computeEquillibriumTemperature = [](State &state, const Time &time, const Vector &absorbedPressure)
      {
        auto &surfaces = state.satellite->surfaces;
        state.time             = time;
        state.absorbedPressure = absorbedPressure;
        state.temperature      = Vector(surfaces.size());
        for(UInt i=0; i<surfaces.size(); i++)
          if(surfaces.at(i).specificHeatCapacity > 0)
            state.temperature(i) = std::pow(LIGHT_VELOCITY/sigma/surfaces.at(i).absorptionInfrared * std::max(state.absorbedPressure(i), 0.), 1./4.);
      };

      // find current state of satellite or create new
      auto state = std::find_if(states.begin(), states.end(), [&](const State &s) {return s.satellite == satellite;});
      if(state == states.end())
      {
        // first state -> set temperature to be in equillibrium with input power
        state = states.insert(states.end(), State());
        state->satellite       = satellite;
        state->needTemperature = std::any_of(satellite->surfaces.begin(), satellite->surfaces.end(), [](const auto &s) {return s.specificHeatCapacity > 0;});
        if(state->needTemperature)
          computeEquillibriumTemperature(*state, time, absorbedPressure);
      }

      // Update temperature
      if(state->needTemperature)
      {
        if((state->time > time) || ((time - state->time).seconds() > 15*60))
        {
          logWarning<<time.dateTimeStr()<<" thermal radiation pressure: temperature state time "<<state->time.dateTimeStr()<<" is not within 15 min before -> restart with equillibrium temperature"<<Log::endl;
          logWarning<<"maybe disable parallel computing?"<<Log::endl;
          computeEquillibriumTemperature(*state, time, absorbedPressure);
        }
        else
        {
          constexpr Double dt = 1.; // integration step size [seconds]
          while(state->time < time)
          {
            axpy(dt/(time-state->time).seconds(), absorbedPressure - state->absorbedPressure, state->absorbedPressure);
            for(UInt i=0; i<satellite->surfaces.size(); i++)
              if(satellite->surfaces.at(i).specificHeatCapacity > 0)
              {
                auto &surface = satellite->surfaces.at(i);
                state->temperature(i) += dt/surface.specificHeatCapacity *
                                        (LIGHT_VELOCITY*state->absorbedPressure(i) - surface.absorptionInfrared * sigma * std::pow(state->temperature(i), 4));

              }
            state->time += seconds2time(dt);
          }
        }
      }

      // compute thermal radition
      for(UInt i=0; i<satellite->surfaces.size(); i++)
      {
        auto &surface = satellite->surfaces.at(i);
        if(surface.specificHeatCapacity > 0)
          F -= factorThermal * (2./3.*sigma/LIGHT_VELOCITY) * surface.area * surface.absorptionInfrared * std::pow(state->temperature(i), 4) * surface.normal;
        else if(surface.specificHeatCapacity < 0)
          F -= factorThermal * (2./3.) * surface.area * absorbedPressure(i) * surface.normal; // spontaneous reemission of absorbed radiation
      }
    } // if(factorThermal)

    return (1./satellite->mass) * rotEarth.rotate(rotSat.rotate(F));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsRadiationPressure::force(SatelliteModelPtr satellite, const Vector3d &direction, Double visible, Double infrared, Vector &absorbedPressure)
{
  try
  {
    Vector3d F;
    for(UInt i=0; i<satellite->surfaces.size(); i++)
    {
      auto &surface = satellite->surfaces.at(i);

      // absorptionX + diffusionX + reflexionX = 1.
      const Double reflexion  = surface.reflexionVisible  * visible + surface.reflexionInfrared  * infrared;
      const Double diffusion  = surface.diffusionVisible  * visible + surface.diffusionInfrared  * infrared;
      const Double absorption = surface.absorptionVisible * visible + surface.absorptionInfrared * infrared;

      if((surface.type == SatelliteModel::Surface::PLATE) || (surface.type == SatelliteModel::Surface::GRACESHADOW))
      {
        const Double cosPhi = -inner(direction, surface.normal);
        if(cosPhi < 1e-10)
          continue;

        // hard-coded GRACE self-shadowing of bottom plate and inner sides
        Double fraction = 1.;
        if(surface.type == SatelliteModel::Surface::GRACESHADOW)
        {
          constexpr Double sin50 = 0.766044443119; // std::sin(DEG2RAD*50.);
          constexpr Double cos50 = 0.642787609687; // std::cos(DEG2RAD*50.);
          constexpr Double a     = 3122.5;         // length of bottom side in x-direction
          constexpr Double b     = 1643.52084;     // width of the bottom side in y-direction
          constexpr Double s     = 73.10281*2;     // length of inner plate = 73.1 mm
          const Double alpha = DEG2RAD*50;
          const Double c     = std::sqrt(std::pow(s,2)+std::pow(b,2)-2*s*b*std::cos(PI-alpha));

          if(surface.normal.z() > 0.999) // bottom plate
          {
            // projected lower point of side plate into bottom side plane along light direction (scaled)
            Double sx = s/a*sin50 * std::fabs(direction.x()/direction.z());
            Double sy = s/b*sin50 * std::fabs(direction.y()/direction.z()) - s/b*cos50;
            if(sx > 1) {sy /= sx; sx = 1.;}
            if(sy > 1) {sx /= sy; sy = 1.;}
            fraction = 1.-std::max(sy, 0.)*(1.-0.5*sx);
          }
          else // inner sides
          {
            Double sy1=0, sy2=1;
            if(std::fabs(direction.z()) > 1e-10)
            {
              // -------------- sy1 --------------
              const Double lambda1 = std::atan2(direction.z(),std::fabs(direction.y()));
              if(lambda1 < 1e-10)             // direction is comming from underneath the satellite
              {
                const Double gamma1 = PI-alpha-std::fabs(lambda1);
                sy1 = -(std::sin(std::fabs(lambda1))*b/std::sin(gamma1))/s;
              }
              else                            // direction is comming from above the satellite
              {
                const Double gamma1 = alpha-lambda1;
                sy1 = (std::sin(lambda1)*b/std::sin(gamma1))/s;
              }
              // -------------- sy2 --------------
              const Double lambda2 = std::atan2(-direction.z(),std::fabs(-direction.y()));
              const Double epsilon = std::asin(sin50*s/c);
              if(lambda2 <= epsilon)
              {
                const Double gamma2 = PI-(PI-epsilon-alpha)-(epsilon-lambda2);
                sy2 = (c*std::sin(epsilon-lambda2)/std::sin(gamma2))/s;
              }
              else
              {
                const Double gamma2 = PI-(epsilon+alpha)-(lambda2-epsilon);
                sy2 = -(c*std::sin(lambda2-epsilon)/std::sin(gamma2))/s;
              }
            }

            const Double dx = std::fabs(direction.x()/(surface.normal.z()*direction.z()+surface.normal.y()*direction.y()));
            constexpr Double h = b*sin50+s*0.984807753012;           //(sin80) height of outer point of other side above inner side
            Double sx1 = b/a*sin50 * dx;                             // bottom
            Double sx2 = h/a * dx;                                   // other inner side

            if(sx1 > 1) {sy1 /= sx1; sx1 = 1.;}
            if(sy1 > 1) {sx1 /= sy1; sy1 = 1.;}
            if(sy1 < 0) {sx1 -= (sx2-sx1)/sy2*sy1; sy1 = 0.;} // final
            if(sx2 > 1) {sy2 = sy1 + (sy2-sy1)*(1-sx1)/(sx2-sx1); sx2 = 1.;}
            if(sy2 > 1) {sx2 = sx1 + (sx2-sx1)*(1-sy1)/(sy2-sy1); sy2 = 1.;}
            fraction = 1. - std::max(sy1, 0.)*(1.-0.5*sx1)- std::max(sy2-sy1, 0.)*((1.-sx1) - 0.5*(sx2-sx1));
          }
        }

        // Source: Montenbruck et al. (2015) Enhanced solar radiation pressure modeling for Galileo satellites. (Equ.5) DOI 10.1007/s00190-014-0774-0
        F += fraction*surface.area*cosPhi*(absorption+diffusion)*direction - fraction*surface.area*cosPhi*(2./3.*diffusion+2.*cosPhi*reflexion)*surface.normal;
        if(absorbedPressure.size())
          absorbedPressure(i) += fraction*cosPhi*absorption;
      }
      else if(surface.type == SatelliteModel::Surface::CYLINDER)
      {
        // Source: Rodriguez Solano (2014) Impact of non-conservative force modeling on GNSS satellite orbits and global solutions. PhD thesis
        const Double cosPhi = -inner(direction, surface.normal);
        if(cosPhi < 1e-10)
          continue;
        F += surface.area*cosPhi*(absorption+diffusion)*direction - surface.area*cosPhi*(PI/6.*diffusion+4./3.*cosPhi*reflexion)*surface.normal;
        if(absorbedPressure.size())
          absorbedPressure(i) += cosPhi*absorption;
      }
      else if(surface.type == SatelliteModel::Surface::SPHERE)
      {
        // Source: Sosnica Krzysztof: Determination of Precise Satellite Orbits and Geodetic Parameters using Satellite Laser Ranging
        const Double CR = (1.0+4./9.*surface.diffusionVisible);
        F += CR * surface.area * direction * (visible + infrared);
       }
    } // for(i=surface)

    return F;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
