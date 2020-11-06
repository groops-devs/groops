/***********************************************/
/**
* @file kaula2SigmaPotentialCoefficients.cpp
*
* @brief Create variances according kaulas rule of thumb.
*
* @author Torsten Mayer-Guerr
* @date 2008-08-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create signal standard deviations of potential coefficients according Kaula's rule of thumb
\begin{equation}
  \sigma_n = \frac{f}{n^p},
\end{equation}
with the degree $n$, the \config{factor} $f$, and the \config{power} $p$.

The standard deviations are written as formal errors of
 \configFile{outputfilePotentialCoefficients}{potentialCoefficients}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileSphericalHarmonics.h"

/***** CLASS ***********************************/

/** @brief Create variances according kaulas rule of thumb.
* The coefficients are writen as SphericalHarmonics.
* @ingroup programsGroup */
class Kaula2SigmaPotentialCoefficients
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Kaula2SigmaPotentialCoefficients, SINGLEPROCESS, "Create variances according kaulas rule of thumb", Misc, PotentialCoefficients)

/***********************************************/

void Kaula2SigmaPotentialCoefficients::run(Config &config)
{
  try
  {
    FileName potentialCoefficientsName;
    UInt     minDegree, maxDegree;
    Double   GM, R;
    Double   power, factor;

    readConfig(config, "outputfilePotentialCoefficients", potentialCoefficientsName, Config::MUSTSET,  "", "");
    readConfig(config, "minDegree", minDegree, Config::DEFAULT,  "2", "");
    readConfig(config, "maxDegree", maxDegree, Config::MUSTSET,  "",  "");
    readConfig(config, "GM",        GM,        Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",         R,         Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "power",     power,     Config::DEFAULT,  "2",    "sigma = factor/degree^power");
    readConfig(config, "factor",    factor,    Config::DEFAULT,  "1e-5", "sigma = factor/degree^power");
    if(isCreateSchema(config)) return;

    logStatus<<"create potential coefficients"<<Log::endl;
    Matrix cnm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      const Double v = factor/std::pow(n, power);
      sigma2cnm(n,0) = v*v;
      for(UInt m=1; m<=n; m++)
        sigma2cnm(n,m) = sigma2snm(n,m) = v*v;
    }

    logStatus<<"writing potential coefficients to file <"<<potentialCoefficientsName<<">"<<Log::endl;
    writeFileSphericalHarmonics(potentialCoefficientsName, SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
