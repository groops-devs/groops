/***********************************************/
/**
* @file GriddedDataTimeSeries2Era5Flux.cpp
*
* @brief Convert ERA5 accumulated TOA radiation into reflected shortwave and emitted longwave flux.
*
* @author Alexander Pamler
* @date 2026-10-12
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts a \file{griddedDataTimeSeries}{griddedDataTimeSeries} of ERA5 accumulated top-of-atmosphere
radiation (e.g. from \program{NetCdf2GriddedDataTimeSeries}) into time series of reflected shortwave flux
and emitted longwave flux (W/m$^2$) as input for the radiation pressure model.

The input must contain TOA incident solar radiation (tisr), top net solar radiation (tsr) and
top net thermal radiation (ttr) in J/m$^2$, accumulated over \config{accumulationInterval}
and time tagged at the end of each accumulation interval (ERA5 valid_time), with spline degree 1.

Reflected shortwave flux = (tisr - tsr)/accumulationInterval,
emitted longwave flux    = -ttr/accumulationInterval.

The values are averaged area-weighted onto a regular cell-centered grid with \config{resolution}
at ellipsoidal \config{height}. The output is written with spline degree 0 (constant within each
accumulation interval), the epochs are the interval boundaries.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedDataTimeSeries.h"

/***** CLASS ***********************************/

/** @brief Convert ERA5 accumulated TOA radiation into reflected shortwave and emitted longwave flux.
* @ingroup programsGroup */
class GriddedDataTimeSeries2Era5Flux
{
  struct Contribution
  {
    UInt   index;  // output cell
    Double weight; // overlap area on unit sphere
  };

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedDataTimeSeries2Era5Flux, SINGLEPROCESS, "Convert ERA5 accumulated TOA radiation to shortwave/longwave flux", Grid, TimeSeries)

/***********************************************/

void GriddedDataTimeSeries2Era5Flux::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutSw, fileNameOutLw, fileNameIn;
    UInt     columnTisr, columnTsr, columnTtr;
    Double   accumulation, inputResolution, resolution, height;

    readConfig(config, "outputfileShortwaveFlux",        fileNameOutSw,   Config::MUSTSET, "",      "reflected shortwave flux [W/m^2]");
    readConfig(config, "outputfileLongwaveFlux",         fileNameOutLw,   Config::MUSTSET, "",      "emitted longwave flux [W/m^2]");
    readConfig(config, "inputfileGriddedDataTimeSeries", fileNameIn,      Config::MUSTSET, "",      "ERA5 accumulations [J/m^2]");
    readConfig(config, "columnTisr",                     columnTisr,      Config::DEFAULT, "0",     "data column of TOA incident solar radiation");
    readConfig(config, "columnTsr",                      columnTsr,       Config::DEFAULT, "1",     "data column of top net solar radiation");
    readConfig(config, "columnTtr",                      columnTtr,       Config::DEFAULT, "2",     "data column of top net thermal radiation");
    readConfig(config, "accumulationInterval",           accumulation,    Config::DEFAULT, "3600",  "[seconds] accumulation period of the input");
    readConfig(config, "inputResolution",                inputResolution, Config::DEFAULT, "0.25",  "[degree] grid spacing of the input");
    readConfig(config, "resolution",                     resolution,      Config::DEFAULT, "1.0",   "[degree] cell size of the output grid");
    readConfig(config, "height",                         height,          Config::DEFAULT, "30000", "[m] ellipsoidal height of the output grid points");
    if(isCreateSchema(config)) return;

    // ============================================

    logStatus<<"read input <"<<fileNameIn<<">"<<Log::endl;
    InFileGriddedDataTimeSeries fileIn(fileNameIn);
    const GriddedData        gridIn  = fileIn.grid();
    const std::vector<Time> &timesIn = fileIn.times();

    if(fileIn.splineDegree() != 1)
      throw(Exception("expected spline degree 1 (one node per epoch), got "+fileIn.splineDegree()%"%i"s));
    if(std::max({columnTisr, columnTsr, columnTtr}) >= fileIn.dataCount())
      throw(Exception("data column index exceeds data count "+fileIn.dataCount()%"%i"s));
    for(UInt i=1; i<timesIn.size(); i++)
      if(std::fabs((timesIn.at(i)-timesIn.at(i-1)).seconds() - accumulation) > 1e-3)
        throw(Exception("epochs must be continuous with spacing accumulationInterval, gap at "+timesIn.at(i).dateTimeStr()));

    // ============================================

    // output grid: regular, cell-centered, north to south, west to east
    // ----------------------------------------------------------------
    const UInt rows = static_cast<UInt>(std::round(180./resolution));
    const UInt cols = static_cast<UInt>(std::round(360./resolution));

    GriddedData gridOut;
    gridOut.ellipsoid = gridIn.ellipsoid;
    for(UInt i=0; i<rows; i++)
    {
      const Double phiTop = 90.-i*resolution;
      const Double phiBot = phiTop-resolution;
      const Double area   = DEG2RAD*resolution * (std::sin(DEG2RAD*phiTop)-std::sin(DEG2RAD*phiBot)); // unit sphere
      for(UInt k=0; k<cols; k++)
      {
        const Double lambda = -180.+(k+0.5)*resolution;
        gridOut.points.push_back(gridOut.ellipsoid(Angle(DEG2RAD*lambda), Angle(DEG2RAD*(phiTop-0.5*resolution)), height));
        gridOut.areas.push_back(area);
      }
    }

    // overlap of input cells with output cells (computed once, grid is fixed)
    // -----------------------------------------------------------------------
    logStatus<<"compute overlap of "<<gridIn.points.size()<<" input cells with "<<gridOut.points.size()<<" output cells"<<Log::endl;
    std::vector<std::vector<Contribution>> contributions(gridIn.points.size());
    for(UInt idx=0; idx<gridIn.points.size(); idx++)
    {
      Angle L, B;
      Double h;
      gridIn.ellipsoid(gridIn.points.at(idx), L, B, h);
      const Double lambda = std::round(RAD2DEG*L/inputResolution)*inputResolution;
      const Double phi    = std::round(RAD2DEG*B/inputResolution)*inputResolution;

      // input cell bounds, poles clipped
      const Double phiTop = std::min( 90., phi+0.5*inputResolution);
      const Double phiBot = std::max(-90., phi-0.5*inputResolution);
      Double lamWest = lambda-0.5*inputResolution;
      lamWest -= 360.*std::floor((lamWest+180.)/360.); // -> [-180, 180)
      const Double lamEast = lamWest+inputResolution;  // may exceed 180 -> wraps

      const UInt iStart = static_cast<UInt>(std::max(0., std::floor((90.-phiTop)/resolution)));
      const UInt iEnd   = std::min(rows-1, static_cast<UInt>(std::floor((90.-phiBot)/resolution)));
      const UInt kStart = static_cast<UInt>(std::floor((lamWest+180.)/resolution));
      const UInt kEnd   = static_cast<UInt>(std::floor((lamEast+180.)/resolution));

      for(UInt i=iStart; i<=iEnd; i++)
      {
        const Double rowTop = 90.-i*resolution;
        const Double rowBot = rowTop-resolution;
        const Double dSinPhi = std::sin(DEG2RAD*std::min(phiTop, rowTop)) - std::sin(DEG2RAD*std::max(phiBot, rowBot));
        if(dSinPhi <= 0)
          continue;
        for(UInt k=kStart; k<=kEnd; k++)
        {
          const Double colWest = -180.+k*resolution;
          const Double dLambda  = std::min(lamEast, colWest+resolution) - std::max(lamWest, colWest);
          if(dLambda <= 0)
            continue;
          contributions.at(idx).push_back({i*cols + (k%cols), DEG2RAD*dLambda*dSinPhi});
        }
      }
    }

    // ============================================

    // convert epoch by epoch
    // ----------------------
    const UInt nodeCount = fileIn.nodeCount();
    std::vector<Matrix> dataSw(nodeCount), dataLw(nodeCount);
    UInt   emptyCells = 0;
    Double meanSw = 0, meanLw = 0;
    const Double totalArea = std::accumulate(gridOut.areas.begin(), gridOut.areas.end(), 0.);

    logStatus<<"convert "<<nodeCount<<" epochs"<<Log::endl;
    for(UInt idNode=0; idNode<nodeCount; idNode++)
    {
      const Matrix in = fileIn.data(idNode);

      Vector sumSw(gridOut.points.size()), sumLw(gridOut.points.size());
      Vector weightSw(gridOut.points.size()), weightLw(gridOut.points.size());
      for(UInt idx=0; idx<in.rows(); idx++)
      {
        const Double valueSw = (in(idx, columnTisr) - in(idx, columnTsr))/accumulation;
        const Double valueLw = -in(idx, columnTtr)/accumulation;
        for(const auto &c : contributions.at(idx))
        {
          if(std::isfinite(valueSw)) {sumSw(c.index) += c.weight*valueSw; weightSw(c.index) += c.weight;}
          if(std::isfinite(valueLw)) {sumLw(c.index) += c.weight*valueLw; weightLw(c.index) += c.weight;}
        }
      }

      dataSw.at(idNode) = Matrix(gridOut.points.size(), 1);
      dataLw.at(idNode) = Matrix(gridOut.points.size(), 1);
      for(UInt i=0; i<gridOut.points.size(); i++)
      {
        if((weightSw(i) <= 0) || (weightLw(i) <= 0))
          emptyCells++;
        dataSw.at(idNode)(i, 0) = (weightSw(i) > 0) ? sumSw(i)/weightSw(i) : 0.;
        dataLw.at(idNode)(i, 0) = (weightLw(i) > 0) ? sumLw(i)/weightLw(i) : 0.;
        meanSw += gridOut.areas.at(i)*dataSw.at(idNode)(i, 0);
        meanLw += gridOut.areas.at(i)*dataLw.at(idNode)(i, 0);
      }

      if((idNode+1)%24 == 0 || idNode+1 == nodeCount)
        logInfo<<"  "<<idNode+1<<" of "<<nodeCount<<" epochs"<<Log::endl;
    }

    // sanity check: global means should be roughly SW ~ 100 W/m^2, LW ~ 240 W/m^2
    logInfo<<"global mean reflected shortwave flux: "<<meanSw/totalArea/nodeCount<<" W/m^2"<<Log::endl;
    logInfo<<"global mean emitted longwave flux:    "<<meanLw/totalArea/nodeCount<<" W/m^2"<<Log::endl;
    if(emptyCells)
      logWarning<<emptyCells<<" output cells without valid input (set to 0)"<<Log::endl;

    // ============================================

    // epochs: interval boundaries (spline degree 0 -> nodeCount+1 epochs)
    // node k = accumulation interval ending at timesIn[k]
    std::vector<Time> timesOut;
    timesOut.push_back(timesIn.front() - seconds2time(accumulation));
    timesOut.insert(timesOut.end(), timesIn.begin(), timesIn.end());

    logStatus<<"write shortwave flux to <"<<fileNameOutSw<<">"<<Log::endl;
    writeFileGriddedDataTimeSeries(fileNameOutSw, 0, timesOut, gridOut, dataSw);
    logStatus<<"write longwave flux to <"<<fileNameOutLw<<">"<<Log::endl;
    writeFileGriddedDataTimeSeries(fileNameOutLw, 0, timesOut, gridOut, dataLw);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/