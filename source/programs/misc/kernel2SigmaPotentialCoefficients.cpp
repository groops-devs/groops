/***********************************************/
/**
* @file kernel2SigmaPotentialCoefficients.cpp
*
* @brief Create variances of spherical harmonics by convolution a kernel with white noise.
* The coefficients are writen as SphericalHarmonics.
*
* @author Torsten Mayer-Guerr
* @date 2008-08-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create variances of spherical harmonics by convolution a kernel with white noise,
e.g. to display filter coefficients of a Gaussian filter.
The coefficients are written as formal errors of \configFile{outputfilePotentialCoefficients}{potentialCoefficients}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/sphericalHarmonics.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Create variances of spherical harmonics by convolution a kernel with white noise.
* The coefficients are writen as SphericalHarmonics.
* @ingroup programsGroup */
class Kernel2SigmaPotentialCoefficients
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Kernel2SigmaPotentialCoefficients, SINGLEPROCESS, "create variances of spherical harmonics by convolution a kernel with white noise", Misc, PotentialCoefficients)

/***********************************************/

void Kernel2SigmaPotentialCoefficients::run(Config &config)
{
  try
  {
    FileName  potentialCoefficientsName;
    KernelPtr kernel;
    UInt      minDegree, maxDegree = INFINITYDEGREE;
    Double    GM, R;
    Double    factor;

    readConfig(config, "outputfilePotentialCoefficients", potentialCoefficientsName, Config::MUSTSET, "", "");
    readConfig(config, "kernel",    kernel,    Config::MUSTSET,  "",  "");
    readConfig(config, "minDegree", minDegree, Config::DEFAULT,  "2", "");
    readConfig(config, "maxDegree", maxDegree, Config::OPTIONAL, "",  "");
    readConfig(config, "GM",        GM,        Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",         R,         Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "factor",    factor,    Config::DEFAULT,  "1", "");
    if(isCreateSchema(config)) return;

    logStatus<<"create potential coefficients"<<Log::endl;
    Matrix cnm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    const Vector kn = kernel->coefficients(Vector3d(0,0,R), maxDegree);
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      const Double v = std::pow(factor*R/GM*kn(n), 2)/(2.*n+1.);
      sigma2cnm(n,0) = v;
      for(UInt m=1; m<=n; m++)
        sigma2cnm(n,m) = sigma2snm(n,m) = v;
    }

    logStatus<<"writing potential coefficients to file"<<Log::endl;
    writeFileSphericalHarmonics(potentialCoefficientsName, SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
