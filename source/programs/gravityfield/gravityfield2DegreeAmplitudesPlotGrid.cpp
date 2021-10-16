/***********************************************/
/**
* @file gravityfield2DegreeAmplitudesPlotGrid.cpp
*
* @brief Time series of degree amplitudes from a time variable gravityfield.
*
* @author Torsten Mayer-Guerr
* @date 2016-01-22
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes a \configClass{timeSeries}{timeSeriesType}
of a time variable \configClass{gravityfield}{gravityfieldType} and saves it as degree amplitudes.
The expansion is limited in the range between \config{minDegree} and \config{maxDegree} inclusivly
\begin{equation}
  \sigma_n = \frac{GM}{R}\left(\frac{R}{r}\right)^{n+1}k_n\sqrt{\sum_{m=0}^n c_{nm}^2+s_{nm}^2}.
\end{equation}

The \configFile{outputfileTimeSeries}{matrix} is a matrix with
every row containing the time, degree, degree amplitude, and the formal error.

To visualize the results use \program{PlotGraph}).

See also \program{Gravityfield2DegreeAmplitudes}.

\fig{!hb}{0.5}{gravityfield2DegreeAmplitudesPlotGrid}{fig:gravityfield2DegreeAmplitudesPlotGrid}{Degree amplitudes of monthly ITSG-Grace2016 solutions relative to GOCO05s.}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/kernel/kernel.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Time series of degree amplitudes from a time variable gravityfield.
* @ingroup programsGroup */
class Gravityfield2DegreeAmplitudesPlotGrid
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2DegreeAmplitudesPlotGrid, SINGLEPROCESS, "Time series of degree amplitudes from a time variable gravityfield.", Gravityfield, TimeSeries)
GROOPS_RENAMED_PROGRAM(Gravityfield2DegreeAmplitudeTimeSeries, Gravityfield2DegreeAmplitudesPlotGrid, date2time(2020, 07, 21))

/***********************************************/

void Gravityfield2DegreeAmplitudesPlotGrid::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName        outName;
    FileName        outNameDegree;
    GravityfieldPtr gravityfield;
    KernelPtr       kernel;
    TimeSeriesPtr   timeSeries;
    Angle           gap;
    UInt            minDegree, maxDegree;
    Double          evalRadius = NAN_EXPR;
    Double          GM, R;

    readConfig(config, "outputfileTimeSeries", outName,      Config::MUSTSET,  "",   "each row: mjd, degree, amplitude, formal error");
    readConfig(config, "gravityfield",         gravityfield, Config::MUSTSET,  "",   "");
    readConfig(config, "kernel",               kernel,       Config::MUSTSET,  "",   "");
    readConfig(config, "timeSeries",           timeSeries,   Config::MUSTSET,  "",   "");
    readConfig(config, "evaluationRadius",     evalRadius,   Config::OPTIONAL, "",   "evaluate the gravity field at this radius (default: evaluate at surface");
    readConfig(config, "polarGap",             gap,          Config::DEFAULT,  "0.", "exclude polar regions (aperture angle in degrees)");
    readConfig(config, "minDegree",            minDegree,    Config::MUSTSET,  "0",  "minimal degree");
    readConfig(config, "maxDegree",            maxDegree,    Config::MUSTSET,  "",   "maximal degree");
    readConfig(config, "GM",                   GM,           Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                    R,            Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    if(std::isnan(evalRadius)) evalRadius = R;
    const std::vector<Time> times = timeSeries->times();
    const Vector kn = kernel->inverseCoefficients(Vector3d(0, 0, evalRadius), maxDegree);

    // time series of degree amplituds as matrix
    //------------------------------------------
    Matrix A(times.size()*(maxDegree+1-minDegree), 4); // mjd, degree, amplitude, formal error
    UInt idx = 0;
    for(UInt i=0; i<times.size(); i++)
    {
      SphericalHarmonics harm = gravityfield->sphericalHarmonics(times.at(i), maxDegree, minDegree, GM, R);
      for(UInt n=minDegree; n<=maxDegree; n++)
      {
        const UInt   minOrder   = static_cast<UInt>(gap*static_cast<Double>(n)+0.5); // Sneeuw
        const Double areaFactor = (minOrder > 0) ? ((2.*n+1.)/(2.*n+2.-2.*minOrder)) : (1.0);
        const Double factor     = areaFactor * std::pow(harm.GM()/harm.R() * std::pow(harm.R()/evalRadius, n+1) * kn(n), 2);

        Double ampl = 0;
        for(UInt m=minOrder; m<=n; m++)
          ampl += factor * (std::pow(harm.cnm()(n,m), 2) + std::pow(harm.snm()(n,m), 2));

        Double sigma = 0;
        if(harm.sigma2cnm().size())
          for(UInt m=minOrder; m<=n; m++)
            sigma += factor * (harm.sigma2cnm()(n,m) + harm.sigma2snm()(n,m));

        A(idx, 0) = times.at(i).mjd();
        A(idx, 1) = n;
        A(idx, 2) = std::sqrt(ampl);
        A(idx, 3) = std::sqrt(sigma);
        idx++;
      }
    }

    // write file
    //-----------
    logStatus<<"write time series of degree amplitudes to file <"<<outName<<">"<<Log::endl;
    writeFileMatrix(outName, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
