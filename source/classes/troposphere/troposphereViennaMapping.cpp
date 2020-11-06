/***********************************************/
/**
* @file troposphereViennaMapping.cpp
*
* @brief Vienna Mapping Functions.
* @see Troposphere
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-04-03
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "troposphere.h"
#include "troposphereViennaMapping.h"

/***********************************************/

TroposphereViennaMapping::TroposphereViennaMapping(Config &config)
{
  try
  {
    readConfig(config, "inputfileVmfCoefficients", fileNamesCoefficients, Config::MUSTSET, "{groopsDataDir}/troposphere/", "ah, aw, zhd, zwd coefficients");
    readConfig(config, "inputfileGpt",             fileNameGpt,           Config::MUSTSET, "{groopsDataDir}/troposphere/gpt3_grid1deg.dat",  "gridded GPT data");
    readConfig(config, "aHeight",                  a_ht,                  Config::DEFAULT, "2.53e-5", "parameter a (height correction)");
    readConfig(config, "bHeight",                  b_ht,                  Config::DEFAULT, "5.49e-3", "parameter b (height correction)");
    readConfig(config, "cHeight",                  c_ht,                  Config::DEFAULT, "1.14e-3", "parameter c (height correction)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TroposphereViennaMapping::init(const std::vector<Vector3d> &stationPositions)
{
  try
  {
    initEmpiricalCoefficients(fileNameGpt, stationPositions);

    GriddedData grid;
    // count epochs
    times.clear();
    for(const FileName &fileName : fileNamesCoefficients)
    {
      InFileGriddedDataTimeSeries file(fileName);
      times.insert(times.end(), file.times().begin(), file.times().end());
      if(!grid.points.size())
        grid = file.grid();
      if(!std::equal(grid.points.begin(), grid.points.end(), file.grid().points.begin(), file.grid().points.end(),
                     [](const Vector3d &p1, const Vector3d &p2) {return (p1-p2).r() < 0.1;}))
        throw(Exception(fileName.str()+": VMF grids differ"));
      if(file.dataCount() != 8)
        throw(Exception(fileName.str()+": data columns "+file.dataCount()%"%i != 8"s));
      if(file.splineDegree() != 1)
        throw(Exception(fileName.str()+": splineDegree "+file.splineDegree()%"%i != 1"s));
    }

    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end()); // remove duplicates
    sampling  = medianSampling(times).mjd();
    if(!isRegular(times))
      throw(Exception("VMF data epochs are uncontinuous"));

    // test grid
    std::vector<std::vector<UInt>>   index(stationPositions.size());
    std::vector<std::vector<Double>> factor(stationPositions.size());
    GriddedDataRectangular gridRectangular;
    if(gridRectangular.init(grid))
    {
      // bilinear interpolation factors and indices
      for(UInt stationId=0; stationId<stationPositions.size(); stationId++)
      {
        Double tx   = (gridRectangular.latitudes.front()-latitude(stationId)) * (gridRectangular.latitudes.size()-1)/(gridRectangular.latitudes.front()-gridRectangular.latitudes.back());
        Double ty   = std::fmod(longitude(stationId)-gridRectangular.longitudes.front()+4*PI, 2*PI)/(2*PI) * gridRectangular.longitudes.size();
        UInt   idxX = std::min(static_cast<UInt>(std::max(std::floor(tx), 0.)), gridRectangular.latitudes.size()-2);
        UInt   idxY = static_cast<UInt>(std::floor(ty));
        tx -= idxX; // latitudinal  interpolation factor [0..1]
        ty -= idxY; // longitudinal interpolation factor [0..1]

        index.at(stationId).push_back(gridRectangular.longitudes.size()*(idxX  ) + (idxY)  %gridRectangular.longitudes.size());
        index.at(stationId).push_back(gridRectangular.longitudes.size()*(idxX  ) + (idxY+1)%gridRectangular.longitudes.size());
        index.at(stationId).push_back(gridRectangular.longitudes.size()*(idxX+1) + (idxY)  %gridRectangular.longitudes.size());
        index.at(stationId).push_back(gridRectangular.longitudes.size()*(idxX+1) + (idxY+1)%gridRectangular.longitudes.size());

        factor.at(stationId).push_back((1-tx) * (1-ty));
        factor.at(stationId).push_back((1-tx) * (ty)  );
        factor.at(stationId).push_back((tx)   * (1-ty));
        factor.at(stationId).push_back((tx)   * (ty)  );
      }
    }
    else
    {
      // find closest station
      for(UInt stationId=0; stationId<stationPositions.size(); stationId++)
      {
        std::vector<Double> distances(grid.points.size());
        for(UInt i=0; i<distances.size(); i++)
          distances.at(i) = (stationPositions.at(stationId)-grid.points.at(i)).r();
        auto iter = std::min_element(distances.begin(), distances.end());
        if(iter == distances.end() || *iter > 10e3)
        {
          Angle  lon, lat;
          Double h;
          grid.ellipsoid(stationPositions.at(stationId), lon, lat, h);
          throw(Exception("no troposphere data for station id "+stationId%"%i (L="s+(lon*RAD2DEG)%"%f, B="s+(lat*RAD2DEG)%"%f)"s));
        }
        index.at(stationId).push_back(std::distance(distances.begin(), iter));
        factor.at(stationId).push_back(1.0);
      }
    }

    ah  = Matrix(times.size(), stationPositions.size());
    aw  = Matrix(times.size(), stationPositions.size());
    zhd = Matrix(times.size(), stationPositions.size());
    zwd = Matrix(times.size(), stationPositions.size());
    gnh = Matrix(times.size(), stationPositions.size());
    geh = Matrix(times.size(), stationPositions.size());
    gnw = Matrix(times.size(), stationPositions.size());
    gew = Matrix(times.size(), stationPositions.size());

    for(const FileName &fileName : fileNamesCoefficients)
    {
      InFileGriddedDataTimeSeries file(fileName);

      // loop over all epochs in VMF troposphere data file
      for(UInt k=0; k<file.times().size(); k++)
      {
        const UInt   idEpoch = std::distance(times.begin(), std::lower_bound(times.begin(), times.end(), file.times().at(k)));
        const Matrix data    = file.data(k);

        for(UInt stationId=0; stationId<stationPositions.size(); stationId++)
        {
          for(UInt i=0; i<index.at(stationId).size(); i++)
          {
            ah (idEpoch, stationId) += factor.at(stationId).at(i) * data(index.at(stationId).at(i), 0);
            aw (idEpoch, stationId) += factor.at(stationId).at(i) * data(index.at(stationId).at(i), 1);
            zhd(idEpoch, stationId) += factor.at(stationId).at(i) * data(index.at(stationId).at(i), 2);
            zwd(idEpoch, stationId) += factor.at(stationId).at(i) * data(index.at(stationId).at(i), 3);
            gnh(idEpoch, stationId) += factor.at(stationId).at(i) * data(index.at(stationId).at(i), 4);
            geh(idEpoch, stationId) += factor.at(stationId).at(i) * data(index.at(stationId).at(i), 5);
            gnw(idEpoch, stationId) += factor.at(stationId).at(i) * data(index.at(stationId).at(i), 6);
            gew(idEpoch, stationId) += factor.at(stationId).at(i) * data(index.at(stationId).at(i), 7);
          }

          // interpolate to height
          if(index.at(stationId).size() > 1)
          {
            constexpr Double g    = 9.80665;   // mean gravity in [m/s^2]
            constexpr Double dMtr = 2.8965e-2; // molar mass of dry air in [kg/mol]
            constexpr Double Rg   = 8.3143;    // universal gas constant in [J/K/mol]

            computeEmpiricalCoefficients(times.at(idEpoch));

            const Double h  = height(stationId) - topo(stationId);
            const Double Tv = T(stationId) * (1. + 0.6077*Q(stationId));  // virtual temperature in [Kelvin]
//             const Double Tv = (T(stationId) + dT(stationId) * hRed) * (1. + 0.378*Q(stationId)/(0.622 + 0.378*Q(stationId))); // virtual temperature in [Kelvin]
            zhd(idEpoch, stationId) *= std::exp(-g*dMtr/Rg/Tv*h);
            zwd(idEpoch, stationId) *= std::exp(-h/2000);
          }
        }
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TroposphereViennaMapping::findIndex(const Time &time, UInt &idx, Double &tau) const
{
  try
  {
    // interpolate coefficients
    tau = (time-times.front()).mjd()/sampling;
    if(tau < 0)
    {
      if((time-times.front()).seconds()<-1.0) // possible clock error
        throw(Exception("time out of range: "+time.dateTimeStr()));
      tau = 0;
      idx = 0;
      return;
    }
    idx = static_cast<UInt>(std::floor(tau));
    tau -= idx;
    if(idx == times.size()-1)
    {
      idx -= 1;
      tau += 1.;
    }
    if(idx >= times.size()-1)
      throw(Exception("time out of range: "+time.dateTimeStr()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TroposphereViennaMapping::slantDelay(const Time &time, UInt stationId, Angle azimuth, Angle elevation) const
{
  try
  {
    UInt   idx;
    Double tau;
    findIndex(time, idx, tau);
    computeEmpiricalCoefficients(time);

    const Double _ah  = (1-tau) * ah (idx, stationId) + tau * ah (idx+1, stationId);
    const Double _aw  = (1-tau) * aw (idx, stationId) + tau * aw (idx+1, stationId);
    const Double _zhd = (1-tau) * zhd(idx, stationId) + tau * zhd(idx+1, stationId);
    const Double _zwd = (1-tau) * zwd(idx, stationId) + tau * zwd(idx+1, stationId);
    const Double _gnh = (1-tau) * gnh(idx, stationId) + tau * gnh(idx+1, stationId);
    const Double _geh = (1-tau) * gnw(idx, stationId) + tau * gnw(idx+1, stationId);
    const Double _gnw = (1-tau) * geh(idx, stationId) + tau * geh(idx+1, stationId);
    const Double _gew = (1-tau) * gew(idx, stationId) + tau * gew(idx+1, stationId);

    const Double sinE = std::sin(elevation);
    const Double vmfh = mappingFunction(sinE, _ah, bh(stationId), ch(stationId))
                      + (1./sinE - mappingFunction(sinE, a_ht, b_ht, c_ht)) * height(stationId)*0.001;
    const Double vmfw = mappingFunction(sinE, _aw, bw(stationId), cw(stationId));
    const Double mfgh = 1. / (sinE*std::tan(elevation) + 0.0031); // hydrostatic gradient mapping function [Chen and Herring, 1997]
    const Double mfgw = 1. / (sinE*std::tan(elevation) + 0.0007); // wet -"-

    return vmfh*_zhd + vmfw*_zwd + (mfgh*_gnh + mfgw*_gnw) * std::cos(azimuth) + (mfgh*_geh + mfgw*_gew) * std::sin(azimuth);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TroposphereViennaMapping::mappingFunctionHydrostatic(const Time &time, UInt stationId, Angle /*azimuth*/, Angle elevation) const
{
  try
  {
    UInt   idx;
    Double tau;
    findIndex(time, idx, tau);
    computeEmpiricalCoefficients(time);

    const Double sinE  = std::sin(elevation);

    return mappingFunction(sinE, (1-tau) * ah(idx, stationId) + tau * ah(idx+1, stationId), bh(stationId), ch(stationId))
           + (1./sinE - mappingFunction(sinE, a_ht, b_ht, c_ht)) * height(stationId)*0.001;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TroposphereViennaMapping::mappingFunctionWet(const Time &time, UInt stationId, Angle /*azimuth*/, Angle elevation) const
{
  try
  {
    UInt   idx;
    Double tau;
    findIndex(time, idx, tau);
    computeEmpiricalCoefficients(time);

    return mappingFunction(std::sin(elevation), (1-tau) * aw(idx, stationId) + tau * aw(idx+1, stationId), bw(stationId), cw(stationId));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TroposphereViennaMapping::mappingFunctionGradient(const Time &/*time*/, UInt /*stationId*/, Angle azimuth, Angle elevation, Double &dx, Double &dy) const
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

void TroposphereViennaMapping::getAprioriValues(const Time &time, UInt stationId, Double &zenithDryDelay, Double &zenithWetDelay,
                                                Double &gradientDryNorth, Double &gradientWetNorth, Double &gradientDryEast, Double &gradientWetEast,
                                                Double &aDry, Double &aWet) const
{
  try
  {
    UInt   idx;
    Double tau;
    findIndex(time, idx, tau);

    aDry             = (1-tau) * ah (idx, stationId) + tau * ah (idx+1, stationId);
    aWet             = (1-tau) * aw (idx, stationId) + tau * aw (idx+1, stationId);
    zenithDryDelay   = (1-tau) * zhd(idx, stationId) + tau * zhd(idx+1, stationId);
    zenithWetDelay   = (1-tau) * zwd(idx, stationId) + tau * zwd(idx+1, stationId);
    gradientDryNorth = (1-tau) * gnh(idx, stationId) + tau * gnh(idx+1, stationId);
    gradientWetNorth = (1-tau) * gnw(idx, stationId) + tau * gnw(idx+1, stationId);
    gradientDryEast  = (1-tau) * geh(idx, stationId) + tau * geh(idx+1, stationId);
    gradientWetEast  = (1-tau) * gew(idx, stationId) + tau * gew(idx+1, stationId);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
