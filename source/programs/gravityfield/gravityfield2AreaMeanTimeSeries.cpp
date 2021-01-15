/***********************************************/
/**
* @file gravityfield2AreaMeanTimeSeries.cpp
*
* @brief Generates a time series as mean values over an area from a time variable gravity field.
*
* @author Torsten Mayer-Guerr
* @date 2006-09-24
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes a time series of time variable
\configClass{gravityfield}{gravityfieldType} functionals averaged over a given area,
e.g. equivalent water heights in the amazon basin.  The type of functional
(e.g gravity anomalies or geoid heights) can be choosen with \configClass{kernel}{kernelType}.
The average is performed at each time step by a weigthed average over all \configClass{grid}{gridType}
points where the weight is the associated area at each point. If \config{removeMean} is set
the temporal mean is removed from the time series. To speed up the computation
the gravity field can be converted to spherical harmonics before the computation
with \config{convertToHarmonics}.

Additionally the root mean square of the values in the area at each time step
can is computed if \config{compueRms} is set.

Additionally the accuracy of the value at each time step can be computed if \config{compueSigma} is set.

The \configFile{outputfileTimeSeries}{instrument} is an instrument file with one, two, or three data columns.
First data column contains the computed functionals and the following columns contain the RMS and the accuracies (optionally).

To visualize the results use \program{PlotGraph}.

\fig{!hb}{0.8}{gravityfield2AreaMeanTimeSeries}{fig:gravityfield2AreaMeanTimeSeries}{Amazon basin: Area mean values of ITSG-Grace2016 monthly solutions with error bars from computed sigma.}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/grid/grid.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Generates a time series as mean values over an area from a time variable gravity field.
* @ingroup programsGroup */
class Gravityfield2AreaMeanTimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2AreaMeanTimeSeries, PARALLEL, "generates a time series as mean values over an area from a time variable gravity field", Gravityfield, TimeSeries)

/***********************************************/

void Gravityfield2AreaMeanTimeSeries::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName        fileNameOut;
    GridPtr         grid;
    TimeSeriesPtr   timeSeries;
    KernelPtr       kernel;
    GravityfieldPtr gravityfield;
    Bool            convertToHarmonics, computeRms, computeSigma, removeMean, multiplyWithArea;

    readConfig(config, "outputfileTimeSeries", fileNameOut,        Config::MUSTSET,  "", "");
    readConfig(config, "grid",                 grid,               Config::MUSTSET,  "", "");
    readConfig(config, "timeSeries",           timeSeries,         Config::MUSTSET,  "", "");
    readConfig(config, "kernel",               kernel,             Config::MUSTSET,  "", "");
    readConfig(config, "gravityfield",         gravityfield,       Config::MUSTSET,  "", "");
    readConfig(config, "convertToHarmonics",   convertToHarmonics, Config::DEFAULT,  "1", "gravityfield is converted to spherical harmonics before evaluation, may accelerate the computation");
    readConfig(config, "multiplyWithArea",     multiplyWithArea,   Config::DEFAULT,  "0", "multiply time series with total area (useful for mass estimates)");
    readConfig(config, "removeMean",           removeMean,         Config::DEFAULT,  "0", "remove the temporal mean of the series");
    readConfig(config, "computeRms",           computeRms,         Config::DEFAULT,  "0", "additional rms each time step");
    readConfig(config, "computeSigma",         computeSigma,       Config::DEFAULT,  "0", "additional error bars at each time step");
    if(isCreateSchema(config)) return;

    std::vector<Time>     times  = timeSeries->times();
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();

    // create grid
    // -----------
    Double totalArea = 0.0;
    for(UInt k=0; k<areas.size(); k++)
      totalArea += areas.at(k);

    for(UInt k=0; k<areas.size(); k++)
      areas.at(k) *= multiplyWithArea ? std::pow(points.at(k).r(), 2) : 1.0/totalArea;

    logInfo<<"  area:             "<<totalArea/(4*PI)*100<<"% of Earth's surface ("<<totalArea*pow(DEFAULT_R/1000,2)<<" km^2)"<<Log::endl;
    logInfo<<"  number of points: "<<points.size()<<Log::endl;

    // Create values on grid
    // ---------------------
    const UInt count = times.size();
    std::vector<Double> value(count);
    std::vector<Double> rms(count, 0.);
    std::vector<Double> sigma(count, 1.0);

    // slow general case
    // -----------------
    if(!convertToHarmonics)
    {
      logStatus<<"create gravity functionals"<<Log::endl;
      Parallel::forEach(value, [&](UInt i)
      {
        Double mean = 0;
        for(UInt k=0; k<points.size(); k++)
          mean += areas.at(k) * gravityfield->field(times.at(i), points.at(k), *kernel);
        return mean;
      }, comm);

      if(computeRms)
      {
        logStatus<<"compute RMS"<<Log::endl;
        Parallel::forEach(rms, [&](UInt i)
        {
          Double rms = 0;
          for(UInt k=0; k<points.size(); k++)
            rms += areas.at(k) * pow(gravityfield->field(times.at(i), points.at(k), *kernel),2);
          return sqrt(rms);
        }, comm);
      }

      if(computeSigma)
      {
        logStatus<<"compute accuracy"<<Log::endl;
        Vector B = areas; // B = linear function from spherical harmonics to block mean
        Parallel::forEach(sigma, [&](UInt i) {return sqrt(inner(B, gravityfield->variance(times.at(i), points, *kernel) * B));}, comm);
      }
    }

    // fast spherical harmonics case
    // -----------------------------
    if(convertToHarmonics)
    {
      logStatus<<"create gravity functionals"<<Log::endl;
      Vector B;
      Parallel::forEach(value, [&](UInt i)
      {
        const SphericalHarmonics harm = gravityfield->sphericalHarmonics(times.at(i));
        const Vector x = harm.x();
        if(B.rows() < x.rows())
          B = MiscGriddedData::synthesisSphericalHarmonicsMatrix(harm.maxDegree(), harm.GM(), harm.R(), points, kernel, harm.isInterior()).trans() * Vector(areas);
        return inner(B.row(0, x.rows()), x);
      }, comm);

      if(computeRms)
      {
        logStatus<<"compute RMS"<<Log::endl;
        Matrix B;
        Parallel::forEach(rms, [&](UInt i)
        {
          const SphericalHarmonics harm = gravityfield->sphericalHarmonics(times.at(i));
          const Vector x = harm.x();
          if(B.columns()<x.rows())
          {
            B = MiscGriddedData::synthesisSphericalHarmonicsMatrix(harm.maxDegree(), harm.GM(), harm.R(), points, kernel, harm.isInterior());
            for(UInt k=0; k<points.size(); k++)
              B.row(k) *= std::sqrt(areas.at(k));
          }
          return std::sqrt(quadsum(B.column(0, x.rows()) * x));
        }, comm);
      }

      if(computeSigma)
      {
        logStatus<<"compute accuracy"<<Log::endl;
logWarning<<"Variance propagation from spherical harmonics to AMV may be not correct implemented. Set convertToHarmonics to FALSE."<<Log::endl;
        Vector B;
        Parallel::forEach(sigma, [&](UInt i)
        {
          const SphericalHarmonics harm = gravityfield->sphericalHarmonics(times.at(i));
          const Matrix C = gravityfield->sphericalHarmonicsCovariance(times.at(i));
          if(B.rows() < C.rows())
            B = MiscGriddedData::synthesisSphericalHarmonicsMatrix(harm.maxDegree(), harm.GM(), harm.R(), points, kernel, harm.isInterior()).trans() * Vector(areas);
          // full covariance matrix
          if(C.getType() == Matrix::SYMMETRIC)
            return sqrt(inner(B.row(0, C.rows()), C * B.row(0, C.rows())));
          // only diagonal matrix
          Double sum = 0;
          for(UInt i=0; i<C.rows(); i++)
            sum += B(i)*C(i,0)*B(i);
          return sqrt(sum);
        }, comm);
      }
    }

    if(!Parallel::isMaster(comm)) return;

    // remove (weigthed) temporal mean
    // -------------------------------
    if(removeMean)
    {
      logStatus<<"remove mean"<<Log::endl;
      Double mean   = 0.0;
      Double weight = 0.0;
      for(UInt i=0; i<count; i++)
      {
        mean   += value.at(i)/pow(sigma.at(i),2);
        weight += 1./pow(sigma.at(i),2);
      }
      mean /= weight;
      for(UInt i=0; i<count; i++)
        value.at(i) -= mean;
    }

    // sort into matrix
    // ----------------
    Matrix A(count, 2+computeRms+computeSigma);
    for(UInt i=0; i<count; i++)
    {
      UInt idx = 0;
      A(i, idx++) = times.at(i).mjd();
      A(i, idx++) = value.at(i);
      if(computeRms)
        A(i, idx++) = rms.at(i);
      if(computeSigma)
        A(i, idx++) = sigma.at(i);
    }

    // write file
    // ----------
    if(!fileNameOut.empty())
    {
      logStatus<<"write time series to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, Arc(times, A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
