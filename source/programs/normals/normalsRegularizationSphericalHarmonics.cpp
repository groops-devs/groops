/***********************************************/
/**
* @file normalsRegularizationSphericalHarmonics.cpp
*
* @brief Diagonal regularization matrix from gravity field accuracies, if not given from signal (cnm,snm), if not given from kaulas rule.
* The inverse accuracies @f$1/\sigma_n^2@f$ are used as weights in the regularization matrix.
* The diagonal is saved as Vector.
*
* @author Torsten Mayer-Guerr
* @date 2008-08-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Diagonal regularization matrix from gravity field accuracies,
if not given from signal (cnm,snm), if not given from kaulas rule.
The inverse accuracies $1/\sigma_n^2$ are used as weights in the regularization matrix.
The diagonal is saved as Vector.

The corresponding pseudo observations can be computed with \program{Gravityfield2SphericalHarmonicsVector}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Diagonal regularization matrix from gravity field accuracies, if not given from signal (cnm,snm), if not given from kaulas rule.
* The inverse accuracies @f$1/\sigma_n^2@f$ are used as weights in the regularization matrix.
* The diagonal is saved as Vector.
* @ingroup programsGroup */
class NormalsRegularizationSphericalHarmonics
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsRegularizationSphericalHarmonics, SINGLEPROCESS, "diagonal regularization matrix from gravity field accuracies, if not given from signal (cnm,snm), if not given from kaulas rule", NormalEquation)

/***********************************************/

void NormalsRegularizationSphericalHarmonics::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName        fileNameDiagonal;
    UInt            minRegularizationDegree = 0;
    UInt            maxRegularizationDegree = INFINITYDEGREE;
    UInt            minDegree, maxDegree;
    GravityfieldPtr gravityfield;
    SphericalHarmonicsNumberingPtr numbering;
    Double          GM, R;
    Bool            makeIsotropic;
    Double          kaulaPower, kaulaFactor;

    readConfig(config, "outputfileDiagonalmatrix", fileNameDiagonal,        Config::MUSTSET,  "",  "");
    readConfig(config, "gravityfield",             gravityfield,            Config::OPTIONAL, "",  "use sigmas, if not given use signal (cnm,snm), if not given use kaulas rule");
    readConfig(config, "minRegularizationDegree",  minRegularizationDegree, Config::OPTIONAL, "",  "");
    readConfig(config, "maxRegularizationDegree",  maxRegularizationDegree, Config::OPTIONAL, "",  "");
    readConfig(config, "minDegree",                minDegree,               Config::MUSTSET,  "2", "");
    readConfig(config, "maxDegree",                maxDegree,               Config::MUSTSET,  "",  "");
    readConfig(config, "numbering",                numbering,               Config::MUSTSET,  "",  "numbering scheme for regul matrix");
    readConfig(config, "GM",                       GM,                      Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                        R,                       Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "makeIsotropic",            makeIsotropic,           Config::DEFAULT,  "0",    "");
    readConfig(config, "kaulaPower",               kaulaPower,              Config::DEFAULT,  "2",    "sigma = kaulaFactor*degree**kaulaPower");
    readConfig(config, "kaulaFactor",              kaulaFactor,             Config::DEFAULT,  "1e-5", "sigma = kaulaFactor*degree**kaulaPower");
    if(isCreateSchema(config)) return;

    minRegularizationDegree = std::max(minRegularizationDegree, minDegree);
    maxRegularizationDegree = std::min(maxRegularizationDegree, maxDegree);
    Matrix cnm2(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm2(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    // coefficients from gravityfield
    // ------------------------------
    logStatus<<"use accuracies, if not given use signal, if not given use kaulas rule"<<Log::endl;
    if(gravityfield)
    {
      // Use variances
      SphericalHarmonics field = gravityfield->sphericalHarmonics(Time(), maxRegularizationDegree, minRegularizationDegree, GM, R);
      cnm2 = field.sigma2cnm();
      snm2 = field.sigma2snm();
      if(cnm2.size() == 0)
        cnm2 = snm2 = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

      // If no variances, use signal instead
      for(UInt n=minRegularizationDegree; n<=maxRegularizationDegree; n++)
        for(UInt m=0; m<=n; m++)
        {
          if(cnm2(n,m) == 0) cnm2(n,m) = std::pow(field.cnm()(n,m),2);
          if(snm2(n,m) == 0) snm2(n,m) = std::pow(field.snm()(n,m),2);
        }
    }

    // Fill the rest with kaula
    for(UInt n=minRegularizationDegree; n<=maxRegularizationDegree; n++)
    {
      if(cnm2(n,0)==0) cnm2(n,0) = std::pow(kaulaFactor/std::pow(n, kaulaPower), 2);
      for(UInt m=1; m<=n; m++)
      {
        if(cnm2(n,m)==0) cnm2(n,m) = std::pow(kaulaFactor/std::pow(n, kaulaPower), 2);
        if(snm2(n,m)==0) snm2(n,m) = std::pow(kaulaFactor/std::pow(n, kaulaPower), 2);
      }
    }

    // All orders set to the mean accuracy
    if(makeIsotropic)
    {
      logStatus<<"make signals isotrop"<<Log::endl;
      for(UInt n=minRegularizationDegree; n<=maxRegularizationDegree; n++)
      {
        Double kn = 0;
        for(UInt m=0; m<=n; m++)
          kn += (cnm2(n,m) + snm2(n,m))/(2*n+1);
        for(UInt m=0; m<=n; m++)
          cnm2(n,m) = snm2(n,m) = kn;
        snm2(n,0) = 0;
      }
    }

    logStatus<<"create diagonal regularization matrix"<<Log::endl;
    Vector regul(numbering->parameterCount(maxDegree, minDegree));
    logInfo<<"  parameters: "<<regul.rows()<<Log::endl;
    std::vector< std::vector<UInt> > idxC, idxS;
    numbering->numbering(maxDegree, minDegree, idxC, idxS);
    for(UInt n=minRegularizationDegree; n<=maxRegularizationDegree; n++)
    {
      if((idxC[n][0]!=NULLINDEX) && (cnm2(n,0)!=0)) regul(idxC[n][0]) = 1/cnm2(n,0);
      for(UInt m=1; m<=maxRegularizationDegree; m++)
      {
        if((idxC[n][m]!=NULLINDEX) && (cnm2(n,m)!=0)) regul(idxC[n][m]) = 1/cnm2(n,m);
        if((idxS[n][m]!=NULLINDEX) && (snm2(n,m)!=0)) regul(idxS[n][m]) = 1/snm2(n,m);
      }
    }

    // write diagonal matrix
    // ---------------------
    logStatus<<"write diagonal matrix to file"<<Log::endl;
    writeFileMatrix(fileNameDiagonal, regul);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
