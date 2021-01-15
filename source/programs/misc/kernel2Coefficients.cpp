/***********************************************/
/**
* @file kernel2Coefficients.cpp
*
* @brief Compute kernel coefficients.
*
* @author Andreas Kvas
* @date 2019-03-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes and returns the coefficients and inverse coefficients of a \configClass{kernel}{kernelType}
from from \config{minDegree} to \config{maxDegree} at a given \config{height}.

The main purpose is for visualization with \program{PlotGraph}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief  Compute kernel coefficients.
* @ingroup programsGroup */
class Kernel2Coefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Kernel2Coefficients, SINGLEPROCESS, "compute kernel coefficients", Misc)
GROOPS_RENAMED_PROGRAM(KernelComputeCoefficients, Kernel2Coefficients, date2time(2020, 9, 1))

/***********************************************/

void Kernel2Coefficients::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName  fileNameOut;
    KernelPtr kernel;
    Double    R, H;
    UInt      minDegree, maxDegree;

    readConfig(config, "outputfileMatrix", fileNameOut,  Config::MUSTSET, "",  "matrix with columns degree, coefficients and inverse coefficients");
    readConfig(config, "kernel",           kernel,       Config::MUSTSET, "",  "");
    readConfig(config, "minDegree",        minDegree,    Config::DEFAULT, "0", "minimum degree of returned coefficients");
    readConfig(config, "maxDegre",         maxDegree,    Config::MUSTSET, "",  "compute coefficients up to maxDegree");
    readConfig(config, "height",           H,            Config::DEFAULT, "0", "evaluate kernel at R+height [m]");
    readConfig(config, "R",                R,            Config::DEFAULT, STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    Vector kn  = kernel->coefficients(Vector3d(0, 0, R), maxDegree);
    Vector kn2 = kernel->inverseCoefficients(Vector3d(0, 0, R), maxDegree);

    Matrix coeffs(maxDegree+1-minDegree, 3);
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      coeffs(n-minDegree, 0) = n;
      coeffs(n-minDegree, 1) = kn(n) * std::pow((R+H)/R, n+1);
      coeffs(n-minDegree, 2) = kn2(n) * std::pow(R/(R+H), n+1);
    }

    logStatus<<"write kernel coefficients to <"<<fileNameOut<<">"<<Log::endl;
    writeFileMatrix(fileNameOut, coeffs);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
