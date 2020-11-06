/***********************************************/
/**
* @file gravityfield2DegreeAmplitudes.cpp
*
* @brief Computes degree amplitudes of a gravity field.
*
* @author Andreas Kvas
* @date 2020-01-18
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes degree amplitudes from a \configClass{gravityfield}{gravityfieldType}
and saves them to a \file{matrix}{matrix} file with three columns: the degree, the degree amplitude, and the formal errors.

The coefficients can be converted to different functionals with \configClass{kernel}{kernelType}.
The gravity field can be evaluated at different altitudes by specifying \config{evaluationRadius}.
Polar regions can be excluded by setting \config{polarGap}.
If set the expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusivly.
The coefficients are related to the reference radius~\config{R}
and the Earth gravitational constant \config{GM}.

See also \program{PotentialCoefficients2DegreeAmplitudes}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Computes degree amplitudes of a gravity field.
* @ingroup programsGroup */
class Gravityfield2DegreeAmplitudes
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2DegreeAmplitudes, SINGLEPROCESS, "computes degree amplitudes of a gravity field", Gravityfield, PotentialCoefficients)

/***********************************************/

void Gravityfield2DegreeAmplitudes::run(Config &config)
{
  try
  {
    enum DegreeType {RMS, CUMMULATE, MEDIAN};
    DegreeType      degreeType = RMS;
    FileName        fileNameCoeff;
    UInt            minDegree, maxDegree = INFINITYDEGREE;
    Time            time;
    Double          GM, R;
    Double          evalRadius = NAN_EXPR;
    GravityfieldPtr gravityfield;
    KernelPtr       kernel;
    Angle           gap;

    readConfig(config, "outputfileMatrix", fileNameCoeff, Config::MUSTSET,  "", "three column matrix with degree, signal amplitude, formal error");
    readConfig(config, "gravityfield",     gravityfield,  Config::MUSTSET,  "", "");
    readConfig(config, "kernel",           kernel,        Config::MUSTSET,  "", "");
    std::string choice;
    if(readConfigChoice(config, "type", choice, Config::MUSTSET, "", "type of variances"))
    {
      if(readConfigChoiceElement(config, "rms",          choice, "degree amplitudes (square root of degree variances)")) degreeType = RMS;
      if(readConfigChoiceElement(config, "accumulation", choice, "cumulate variances over degrees"))                     degreeType = CUMMULATE;
      if(readConfigChoiceElement(config, "median",       choice, "meadian of absolute values per degree"))               degreeType = MEDIAN;
      endChoice(config);
    }
    readConfig(config, "time",             time,         Config::OPTIONAL, "", "at this time the gravity field will be evaluated");
    readConfig(config, "evaluationRadius", evalRadius,   Config::OPTIONAL, "", "evaluate the gravity field at this radius (default: evaluate at surface");
    readConfig(config, "polarGap",         gap,          Config::DEFAULT,  "0.0", "exclude polar regions (aperture angle in degrees)");
    readConfig(config, "minDegree",        minDegree,    Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",        maxDegree,    Config::OPTIONAL, "",  "");
    readConfig(config, "GM",               GM,           Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                R,            Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    // Check input
    // -----------
    if(maxDegree<minDegree)
      throw(Exception("Maximum degree lower than minimum degree (" + maxDegree%"%i"s + " vs. " + minDegree%"%i"s + ")."));

    if(std::isnan(evalRadius)) evalRadius = R;

    // Create potential coefficients
    // -----------------------------
    logStatus<<"create spherical harmonics"<<Log::endl;
    SphericalHarmonics harm = gravityfield->sphericalHarmonics(time, maxDegree, minDegree, GM, R);
    maxDegree = harm.maxDegree();
    const Vector kn     = kernel->inverseCoefficients(Vector3d(0, 0, evalRadius), maxDegree, harm.isInterior());
    const Bool hasSigma = harm.sigma2cnm().size() || harm.sigma2snm().size();

    auto vectorMedian = [](std::vector<Double> &data)
    {
      std::partial_sort(data.begin(), data.begin()+data.size()/2+1, data.end());
      return (data.size()%2) ? data.at(data.size()/2) : (0.5*(data.at(data.size()/2-1)+data.at(data.size()/2)));
    };

    auto vectorAccumulate = [](std::vector<Double> &data)
    {
      std::for_each(data.begin(), data.end(), [](Double &v) { v = v*v; });
      std::partial_sum(data.begin(), data.end(), data.begin());
      std::for_each(data.begin(), data.end(), [](Double &v) { v = std::sqrt(v); });
    };

    std::vector<Double> x, y, sigma;
    for(UInt n=0; n<=maxDegree; n++)
    {
      const UInt   minOrder   = static_cast<UInt>(gap*static_cast<Double>(n)+0.5); // Sneeuw
      const Double areaFactor = (minOrder>0) ? ((2.*n+1.)/(2.*n+2.-2.*minOrder)) : (1.0);
      const Double factor     = areaFactor * std::pow(harm.GM()/harm.R() * std::pow(harm.R()/evalRadius, n+1) * kn(n), 2);

      std::vector<Double> coefficients, formalErrors;
      for(UInt m=minOrder; m<=n; m++)
      {
        coefficients.push_back(factor * std::pow(harm.cnm()(n, m),2));
        if(hasSigma) formalErrors.push_back(factor * harm.sigma2cnm()(n, m));
        if(m > 0)
        {
          coefficients.push_back(factor * std::pow(harm.snm()(n, m),2));
          if(hasSigma) formalErrors.push_back(factor * harm.sigma2snm()(n, m));
        }
      }

      // degree variances
      // ----------------
      x.push_back(n);
      if( (degreeType == RMS) || (degreeType == CUMMULATE) )
      {
        y.push_back( std::sqrt(std::accumulate(coefficients.begin(), coefficients.end(), 0.0)) );
        if(hasSigma) sigma.push_back( std::sqrt(std::accumulate(formalErrors.begin(), formalErrors.end(), 0.0)) );
      }
      else if(degreeType == MEDIAN)
      {
        y.push_back( std::sqrt( (2*n+1) * vectorMedian(coefficients)) );
        if(hasSigma) sigma.push_back( std::sqrt( (2*n+1) * vectorMedian(formalErrors)) );
      }
    } // for(n)

    if(degreeType == CUMMULATE)
    {
      vectorAccumulate(y);
      if(hasSigma) vectorAccumulate(sigma);
    }

    Matrix data(x.size(), 3, NAN_EXPR);
    copy(Vector(x), data.column(0));
    copy(Vector(y), data.column(1));
    if(hasSigma) copy(Vector(sigma), data.column(2));

    // write
    // -----
    logStatus<<"write degree amplitudes to file <"<<fileNameCoeff<<">"<<Log::endl;
    writeFileMatrix(fileNameCoeff, data);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
