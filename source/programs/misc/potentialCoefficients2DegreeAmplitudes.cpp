/***********************************************/
/**
* @file potentialCoefficients2DegreeAmplitudes.cpp
*
* @brief Degree amplitudes of potential coefficients.
*
* @author Torsten Mayer-Guerr
* @date 2020-04-28
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes degree amplitudes from
\file{potentialCoefficients files}{potentialCoefficients}
and saves them to a \file{matrix}{matrix} file.

The coefficients can be filtered with \configClass{filter}{sphericalHarmonicsFilterType} and converted
to different functionals with \configClass{kernel}{kernelType}. The gravity field can be evaluated at
different altitudes by specifying \config{evaluationRadius}. Polar regions can be excluded
by setting \config{polarGap}. If set the expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusivly. The coefficients are related to the reference radius~\config{R}
and the Earth gravitational constant \config{GM}.

The \configFile{outputfileMatrix}{matrix} contains in the first 3 columns the degree, the degree amplitude, and
the formal errors. For each additional \configFile{inputfilePotentialCoefficients}{potentialCoefficients} three columns
are appended: the degree amplitude, the formal errors, and the difference to the first file.

For example the data columns for 4 \configFile{inputfilePotentialCoefficients}{potentialCoefficients} are
\begin{itemize}
\item degree=\verb|data0|
\item PotentialCoefficients0: signal=\verb|data1|, error=\verb|data2|,
\item PotentialCoefficients1: signal=\verb|data3|, error=\verb|data4|,  difference=\verb|data5|,
\item PotentialCoefficients2: signal=\verb|data6|, error=\verb|data7|,  difference=\verb|data8|,
\item PotentialCoefficients3: signal=\verb|data9|, error=\verb|data10|, difference=\verb|data11|.
\end{itemize}

See also \program{Gravityfield2DegreeAmplitudes}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/system.h"
#include "files/fileMatrix.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/kernel/kernel.h"
#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilter.h"

/***** CLASS ***********************************/

/** @brief Degree amplitudes of potential coefficients.
* @ingroup programsGroup */
class PotentialCoefficients2DegreeAmplitudes
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(PotentialCoefficients2DegreeAmplitudes, SINGLEPROCESS, "degree amplitudes of potential coefficients", Misc, Gravityfield, PotentialCoefficients)

/***********************************************/

void PotentialCoefficients2DegreeAmplitudes::run(Config &config)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNamesIn;
    enum DegreeType       {RMS, CUMMULATE};
    KernelPtr             kernel;
    SphericalHarmonicsFilterPtr filter;
    DegreeType            degreeType = RMS;
    Double                r = NAN_EXPR;
    Angle                 gap(0.0);
    UInt                  minDegree, maxDegree = INFINITYDEGREE;
    Double                GM, R;
    std::string           choice;

    readConfig(config, "outputfileMatrix",               fileNameOut, Config::MUSTSET, "", "matrix with degree, signal amplitude, formal error, differences");
    readConfig(config, "inputfilePotentialCoefficients", fileNamesIn, Config::MUSTSET,  "{groopsDataDir}/potential/", "");
    readConfig(config, "kernel",                         kernel,      Config::MUSTSET, "geoidHeight", "");
    readConfig(config, "filter",                         filter,      Config::OPTIONAL, "", "filter the coefficients");
    if(readConfigChoice(config, "type", choice, Config::MUSTSET, "", "type of variances"))
    {
      if(readConfigChoiceElement(config, "rms",          choice, "degree amplitudes (square root of degree variances)")) degreeType = RMS;
      if(readConfigChoiceElement(config, "accumulation", choice, "cumulate variances over degrees"))                     degreeType = CUMMULATE;
      endChoice(config);
    }
    readConfig(config, "evaluationRadius", r,         Config::OPTIONAL, "", "evaluate the gravity field at this radius (default: evaluate at surface");
    readConfig(config, "polarGap",         gap,       Config::DEFAULT,  "0.0", "exclude polar regions (aperture angle in degrees)");
    readConfig(config, "minDegree",        minDegree, Config::DEFAULT,  "2", "");
    readConfig(config, "maxDegree",        maxDegree, Config::OPTIONAL, "",  "");
    readConfig(config, "GM",               GM,        Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                R,         Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;


    std::vector<SphericalHarmonics> harm(fileNamesIn.size());
    for(UInt i=0; i<fileNamesIn.size(); i++)
      if(System::exists(fileNamesIn.at(i)))
      {
        logStatus<<"read potential coefficients <"<<fileNamesIn.at(i)<<">"<<Log::endl;
        readFileSphericalHarmonics(fileNamesIn.at(i), harm.at(i));
        if(maxDegree == INFINITYDEGREE) maxDegree = harm.at(i).maxDegree();
        maxDegree = std::max(maxDegree, harm.at(i).maxDegree());
      }
      else
        logWarning<<"file <"<<fileNamesIn.at(i)<<"> not exist. continue"<<Log::endl;

    if(std::isnan(r)) r = R;


    logStatus<<"compute degree amplitudes"<<Log::endl;
    Matrix A(maxDegree+1-minDegree, 3*harm.size(), NAN_EXPR);
    UInt idx = 1;
    for(UInt i=0; i<harm.size(); i++)
    {
      harm.at(i) = harm.at(i).get(INFINITYDEGREE, minDegree, GM, R);
      if(filter)
        harm.at(i) = filter->filter(harm.at(i));

      // use standard deviation instead of variances
      for(UInt n=0; n<harm.at(i).sigma2cnm().rows(); n++)
        for(UInt m=0; m<=n; m++)
        {
          harm.at(i).sigma2cnm()(n,m) = std::sqrt(harm.at(i).sigma2cnm()(n,m));
          harm.at(i).sigma2snm()(n,m) = std::sqrt(harm.at(i).sigma2snm()(n,m));
        }

      const Vector kn = kernel->inverseCoefficients(Vector3d(r, 0, 0), maxDegree);

      // degree variances
      // ----------------
      auto variances = [&](const_MatrixSliceRef cnm, const_MatrixSliceRef snm)
      {
        Vector variance(maxDegree+1-minDegree, NAN_EXPR);
        for(UInt n=minDegree; n<std::min(maxDegree+1, cnm.rows()); n++)
        {
          const UInt   minOrder   = static_cast<UInt>(gap*static_cast<Double>(n)+0.5); // Sneeuw
          const Double areaFactor = (minOrder>0) ? ((2.*n+1.)/(2.*n+2.-2.*minOrder)) : (1.0);
          const Double factor     = areaFactor * std::pow(GM/R * std::pow(R/r, n+1) * kn(n), 2);
          variance(n-minDegree)   = factor * (quadsum(cnm.slice(n, minOrder, 1, n+1-minOrder)) + quadsum(snm.slice(n, minOrder, 1, n+1-minOrder)));
        } // for(n)
        return variance;
      };
      // ----------------

      copy(variances(harm.at(i).cnm(),       harm.at(i).snm()),       A.column(idx++)); // signal
      copy(variances(harm.at(i).sigma2cnm(), harm.at(i).sigma2snm()), A.column(idx++)); // formal errors
      if(i > 0)                                                                         // differences
      {
        const UInt size = std::min(harm.at(i).maxDegree(), harm.at(0).maxDegree()) + 1;
        copy(variances(harm.at(i).cnm().slice(0,0,size,size)-harm.at(0).cnm().slice(0,0,size,size),
                       harm.at(i).snm().slice(0,0,size,size)-harm.at(0).snm().slice(0,0,size,size)), A.column(idx++));
      }
    }

    if(degreeType == CUMMULATE)
      for(UInt n=1; n<A.rows(); n++)
        axpy(1., A.row(n-1), A.row(n));

    for(UInt n=0; n<A.rows(); n++)
      for(UInt i=0; i<A.columns(); i++)
        A(n,i) = std::sqrt(A(n,i));

    for(UInt n=0; n<A.rows(); n++)
      A(n,0) = n+minDegree;

    // write
    // -----
    logStatus<<"write degree amplitudes to file <"<<fileNameOut<<">"<<Log::endl;
    writeFileMatrix(fileNameOut, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
