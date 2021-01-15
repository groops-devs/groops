/***********************************************/
/**
* @file grs2PotentialCoefficients.cpp
*
* @brief Spherical harmonics from Geodetic Reference System (GRS).
*
* @author Torsten Mayer-Guerr
* @date 2004-12-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program creates potential coefficients from the defining constants
of a Geodetic Reference System (GRS). The potential coeffiencts excludes the centrifugal part.
The default values creates the GRS80.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileSphericalHarmonics.h"

/***** CLASS ***********************************/

/** @brief Spherical harmonics from Geodetic Reference System (GRS).
* @ingroup programsGroup */
class Grs2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Grs2PotentialCoefficients, SINGLEPROCESS, "spherical harmonics from Geodetic Reference System (GRS)", Misc, PotentialCoefficients)

/***********************************************/

void Grs2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outName;
    Double   GM, R, J2, omega;
    UInt     degree;

    readConfig(config, "outputfilePotentialCoefficients", outName, Config::MUSTSET,  "", "");
    readConfig(config, "maxDegree", degree, Config::MUSTSET, "10",           "");
    readConfig(config, "GM",        GM,     Config::DEFAULT, "3.986005e+14", "Geocentric gravitational constant");
    readConfig(config, "R",         R,      Config::DEFAULT, STRING_DEFAULT_GRS80_a,  "reference radius");
    readConfig(config, "J2",        J2,     Config::DEFAULT, "108263e-8",    "Dynamical form factor");
    readConfig(config, "omega",     omega,  Config::DEFAULT, "7292115e-11",  "Angular velocity of rotation");
    if(isCreateSchema(config)) return;

    logStatus<<"create spherical harmonics"<<Log::endl;
    Matrix cnm(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    Double e_sec, q0;
    Double e_old = 1.0;
    Double e     = 0.5;

    do
    {
      e_sec = e * sqrt(1-e*e)/(1-e*e);
      q0    = (1.0 + (3.0/(e_sec*e_sec)))*atan(e_sec)-(3.0/e_sec);
      e_old = e;
      e     = sqrt(3.0*J2 + (4.0*omega*omega*R*R*R * e*e*e) / (15.0*GM*q0));
    }
    while( fabs(e-e_old) > 1.0e-13);

    cnm(0,0) = 1.0;
    for(UInt n=2; n<=degree; n+=2)
    {
      UInt m = n/2;
      Double Jn = pow(-1.,m+1)*((3*pow(e*e,m))/((2.*m+1.)*(2.*m+3.)))*(1.-m+5.*m*(J2/(e*e)));
      cnm(n,0) = -Jn/sqrt(2.*n+1.);
    }

    logStatus<<"writing potential coefficients to file <"<<outName<<">"<<Log::endl;
    writeFileSphericalHarmonics(outName, SphericalHarmonics(GM, R, cnm, snm));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
