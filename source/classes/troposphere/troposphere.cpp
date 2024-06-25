/***********************************************/
/**
* @file troposphere.cpp
*
* @brief Signal delay in the atmosphere.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-04-03
*/
/***********************************************/

#define DOCSTRING_Troposphere

#include "base/import.h"
#include "config/configRegister.h"
#include "files/fileGriddedData.h"
#include "troposphereGpt.h"
#include "troposphereMendesAndPavlis.h"
#include "troposphereSaastamoinen.h"
#include "troposphereHopfield.h"
#include "troposphereViennaMapping.h"
#include "troposphere.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Troposphere, "troposphereType",
                      TroposphereViennaMapping,
                      TroposphereGpt,
                      TroposphereMendesAndPavlis,
                      TroposphereSaastamoinen,
                      TroposphereHopfield)

GROOPS_READCONFIG_CLASS(Troposphere, "troposphereType")

/***********************************************/

TropospherePtr Troposphere::create(Config &config, const std::string &name)
{
  try
  {
    TropospherePtr troposphere;
    std::string type;

    readConfigChoice(config, name, type, Config::MUSTSET, "", "signal delay in the atmosphere");
    if(readConfigChoiceElement(config, "viennaMapping",   type, "Vienna Mapping Function"))
      troposphere = TropospherePtr(new TroposphereViennaMapping(config));
    if(readConfigChoiceElement(config, "gpt",             type, "GPT empirical troposphere model"))
      troposphere = TropospherePtr(new TroposphereGpt(config));
    if(readConfigChoiceElement(config, "mendesAndPavlis", type, "SLR troposphere model by Mendes and Pavlis, 2004"))
      troposphere = TropospherePtr(new TroposphereMendesAndPavlis(config));
    if(readConfigChoiceElement(config, "saastamoinen",    type, "Saastamoinen troposphere model"))
      troposphere = TropospherePtr(new TroposphereSaastamoinen(config));
    if(readConfigChoiceElement(config, "hopfield",        type, "Hopfield troposphere model (simplified)"))
      troposphere = TropospherePtr(new TroposphereHopfield(config));
    endChoice(config);

    return troposphere;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void Troposphere::initEmpiricalCoefficients(const FileName &fileNameGpt, const std::vector<Vector3d> &stationPositions)
{
  try
  {
    GriddedDataRectangular grid;
    readFileGriddedData(fileNameGpt, grid);

    const UInt stationCount = stationPositions.size();
    longitude = Vector(stationCount);
    latitude  = Vector(stationCount);
    height    = Vector(stationCount);
    coeff     = Matrix(grid.values.size(), stationCount);

    for(UInt stationId=0; stationId<stationCount; stationId++)
    {
      Angle  L, B;
      grid.ellipsoid(stationPositions.at(stationId), L, B, height(stationId));
      longitude(stationId) = L;
      latitude(stationId)  = B;

      // interpolation factors and indices
      Double tx   = (grid.latitudes.front()-latitude(stationId)) * (grid.latitudes.size()-1)/(grid.latitudes.front()-grid.latitudes.back());
      Double ty   = std::fmod(longitude(stationId)-grid.longitudes.front()+4*PI, 2*PI)/(2*PI) * grid.longitudes.size();
      UInt   idxX = std::min(static_cast<UInt>(std::max(std::floor(tx), 0.)), grid.latitudes.size()-2);
      UInt   idxY = static_cast<UInt>(std::floor(ty));
      tx -= idxX; // latitudinal  interpolation factor [0..1]
      ty -= idxY; // longitudinal interpolation factor [0..1]

      for(UInt k=0; k<grid.values.size(); k++)
      {
        coeff(k, stationId) = (1-tx) * (1-ty) * grid.values.at(k)(idxX,   (idxY)  %grid.longitudes.size())
                            + (1-tx) * (ty)   * grid.values.at(k)(idxX,   (idxY+1)%grid.longitudes.size())
                            + (tx)   * (1-ty) * grid.values.at(k)(idxX+1, (idxY)  %grid.longitudes.size())
                            + (tx)   * (ty)   * grid.values.at(k)(idxX+1, (idxY+1)%grid.longitudes.size());
      }
    }

    timeRef = Time();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Troposphere::computeEmpiricalCoefficients(const Time &time) const
{
  try
  {
    if(timeRef.mjdInt() == time.mjdInt())
      return; // computing empirical coefficients once per day is sufficient, since they only have annual and semiannual variations

    timeRef = time;

    const Double t = (time.mjd()-J2000)/365.25;
    Vector fourier(5);
    fourier(0) = 1.;               // constant
    fourier(1) = std::cos(2*PI*t); // cos annual
    fourier(2) = std::sin(2*PI*t); // sin annual
    fourier(3) = std::cos(4*PI*t); // cos semiannual
    fourier(4) = std::sin(4*PI*t); // sin semiannual

    topo = ah = aw = bh = bw = ch = cw = p = T = Q = dT = la = Tm = zhd = zwd = gnh = geh = gnw = gew = Vector(longitude.rows());
    for(UInt stationId=0; stationId<longitude.size(); stationId++)
    {
      topo(stationId)  = coeff(0, stationId);                                                // [m]

      ah(stationId)    = inner(fourier, coeff.slice(1+ 0*5,  stationId, fourier.rows(), 1));
      aw(stationId)    = inner(fourier, coeff.slice(1+ 1*5,  stationId, fourier.rows(), 1));
      bh(stationId)    = inner(fourier, coeff.slice(1+ 2*5,  stationId, fourier.rows(), 1));
      bw(stationId)    = inner(fourier, coeff.slice(1+ 3*5,  stationId, fourier.rows(), 1));
      ch(stationId)    = inner(fourier, coeff.slice(1+ 4*5,  stationId, fourier.rows(), 1));
      cw(stationId)    = inner(fourier, coeff.slice(1+ 5*5,  stationId, fourier.rows(), 1));

      p(stationId)     = inner(fourier, coeff.slice(1+ 6*5,  stationId, fourier.rows(), 1)); // [Pa]
      T(stationId)     = inner(fourier, coeff.slice(1+ 7*5,  stationId, fourier.rows(), 1)); // [Kelvin]
      Q(stationId)     = inner(fourier, coeff.slice(1+ 8*5,  stationId, fourier.rows(), 1)); // [kg/kg]
      dT(stationId)    = inner(fourier, coeff.slice(1+ 9*5,  stationId, fourier.rows(), 1)); // tempLapseRate [Kelvin/m]
      la(stationId)    = inner(fourier, coeff.slice(1+10*5,  stationId, fourier.rows(), 1)); // waterWaporDecreaseFactor []
      Tm(stationId)    = inner(fourier, coeff.slice(1+11*5,  stationId, fourier.rows(), 1)); // waterWaporMeanTemperature [Kelvin]

      gnh(stationId)   = inner(fourier, coeff.slice(1+12*5, stationId, fourier.rows(), 1)); // [m]
      geh(stationId)   = inner(fourier, coeff.slice(1+13*5, stationId, fourier.rows(), 1)); // [m]
      gnw(stationId)   = inner(fourier, coeff.slice(1+14*5, stationId, fourier.rows(), 1)); // [m]
      gew(stationId)   = inner(fourier, coeff.slice(1+15*5, stationId, fourier.rows(), 1)); // [m]

      constexpr Double g    = 9.80665;                  // mean gravity in [m/s^2]
      constexpr Double dMtr = 2.8965e-2;                // molar mass of dry air in [kg/mol]
      constexpr Double Rg   = 8.3143;                   // universal gas constant in [J/K/mol]
      constexpr Double k1   = 77.604;                   // coefficient [K/hPa]
      constexpr Double k2   = 64.79;                    // coefficient [K/hPa]
      constexpr Double k2p  = k2 - k1*18.0152/28.9644;  // coefficient [K/hPa]
      constexpr Double k3   = 377600.;                  // coefficient [KK/hPa]

      const Double h   = height(stationId) - topo(stationId);            // reduced height
      const Double Tv = T(stationId) * (1. + 0.6077*Q(stationId));       // virtual temperature in [Kelvin]
      const Double P  = 0.01 * p(stationId) * std::exp(-g*dMtr/Rg/Tv*h); // pressure in [hPa]

      // zenith hydrostatic delay (Saastamoinen, 1972) and zenith wet delay (Askne and Nordius, 1987)
      zhd(stationId) = 0.0022768 * P / (1. - 0.00266*std::cos(2*latitude(stationId)) - 0.00000028*height(stationId));
      zwd(stationId) = 1e-6*Rg/dMtr/g * (k2p + k3/Tm(stationId)) * P * Q(stationId)/(0.622 + 0.378*Q(stationId))/(1+la(stationId));
      zwd(stationId) *= std::exp(-h/2000); // extrapolate ZWD to station height
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
