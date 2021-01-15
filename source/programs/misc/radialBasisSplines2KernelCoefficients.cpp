/***********************************************/
/**
* @file radialBasisSplines2KernelCoefficients.cpp
*
* @brief Kernel/Covariance-function from gravity field accuracies, if not given from signal (cnm,snm), if not given from kaulas rule.
*
* @author Torsten Mayer-Guerr
* @date 2008-08-08
*
*/
/***********************************************/
// Latex documentation

#define DOCSTRING docstring
static const char *docstring = R"(
This program calculates the coefficients $k_n$ of a \configClass{kernel:coefficients}{kernelType:coefficients} according to
\begin{equation}
  k_n = \frac{GM}{4\pi R}\frac{\sigma_n}{\sqrt{2n+1}}.
\end{equation}
from a given \configClass{gravityfield}{gravityfieldType},
with \config{R} and \config{GM} describing the reference radius and the geocentric constant, respectively.
The $\sigma_n$
stand for the gravity field accuracies (from degree \config{minDegree} to \config{maxDegree}), if they are given.
If no accuracies are provided, the $\sigma_n$
represent the square root of the degree variances of the gravity field.
If \config{maxDegree} excedes the maximum degree given by \configClass{gravityfield}{gravityfieldType},
the higher degrees are complemented by Kaula's rule
The output of the coefficients is given in the file  \configFile{outputfileCoefficients}{matrix}.
)";

/***********************************************/

#include "base/import.h"
#include "files/fileMatrix.h"
#include "classes/gravityfield/gravityfield.h"
#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief Kernel/Covariance-function from gravity field accuracies, if not given from signal (cnm,snm), if not given from kaulas rule.
* @ingroup programsGroup */
class RadialBasisSplines2KernelCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(RadialBasisSplines2KernelCoefficients, SINGLEPROCESS, "Kernel/Covariance-function from gravity field accuracies, if not given from signal (cnm,snm), if not given from kaulas rule", Misc)
GROOPS_RENAMED_PROGRAM(KernelDegreeVariances, RadialBasisSplines2KernelCoefficients, date2time(2020, 11, 9))

/***********************************************/

void RadialBasisSplines2KernelCoefficients::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName        fileNameCoeff;
    GravityfieldPtr gravityfield;
    UInt            minDegree, maxDegree = INFINITYDEGREE;
    Double          GM, R;
    Double          kaulaPower, kaulaFactor;

    readConfig(config, "outputfileCoefficients", fileNameCoeff, Config::MUSTSET,  "",  "");
    readConfig(config, "gravityfield",           gravityfield,  Config::OPTIONAL, "",  "use sigmas, if not given use signal (cnm,snm), if not given use kaulas rule");
    readConfig(config, "minDegree",              minDegree,     Config::MUSTSET,  "2", "");
    readConfig(config, "maxDegree",              maxDegree,     Config::OPTIONAL, "",  "");
    readConfig(config, "GM",                     GM,            Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                      R,             Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "kaulaPower",             kaulaPower,    Config::DEFAULT,  "2",    "sigma = kaulaFactor/degree^kaulaPower");
    readConfig(config, "kaulaFactor",            kaulaFactor,   Config::DEFAULT,  "1e-5", "sigma = kaulaFactor/degree^kaulaPower");
    if(isCreateSchema(config)) return;

    logStatus<<"use accuracies, if not given use signal, if not given use kaulas rule"<<Log::endl;
    Vector coeff;
    if(gravityfield)
    {
      // Use variances
      SphericalHarmonics field = gravityfield->sphericalHarmonics(Time(), maxDegree, minDegree, GM, R);
      maxDegree  = field.maxDegree();
      coeff = Vector(maxDegree+1);
      for(UInt n=minDegree; n<=maxDegree; n++)
      {
        Double sum = 0.;
        if(field.sigma2cnm().size())
          for(UInt m=0; m<=n; m++)
            sum += field.sigma2cnm()(n,m) + field.sigma2snm()(n,m);
        // If no variances, use signal instead
        if(sum==0)
          for(UInt m=0; m<=n; m++)
            sum += std::pow(field.cnm()(n,m), 2) + std::pow(field.snm()(n,m), 2);
        coeff(n) = GM/R * std::sqrt(sum/(4*PI)/(2.*n+1.));
      }
    }

    if(coeff.rows() == 0)
    {
      if(maxDegree == INFINITYDEGREE)
        throw(Exception("maxDegree or gravityfield must be given"));
      coeff = Vector(maxDegree+1);
    }

    // Fill the rest with kaula
    for(UInt n=maxDegree+1; n-->minDegree;)
    {
      if(coeff(n)!=0)
        break;
      const Double sum = (2*n+1) * std::pow(kaulaFactor/std::pow(n, kaulaPower), 2);
      coeff(n) = GM/R * std::sqrt(sum/(4*PI)/(2.*n+1.));
    }

    logStatus<<"write coefficient vector to file <"<<fileNameCoeff<<">"<<Log::endl;
    writeFileMatrix(fileNameCoeff, coeff);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
