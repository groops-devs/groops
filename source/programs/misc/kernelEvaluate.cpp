/***********************************************/
/**
* @file kernelEvaluate.cpp
*
* @brief Compute kernel values for distant angles.
*
* @author Torsten Mayer-Guerr
* @date 2017-05-30
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Compute \configClass{kernel}{kernelType} values for distant angles.
The main purpose is for visualization with \program{PlotGraph}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Compute kernel values for distant angles.
* @ingroup programsGroup */
class KernelEvaluate
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(KernelEvaluate, SINGLEPROCESS, "Compute kernel values for distant angles", Misc)

/***********************************************/

void KernelEvaluate::run(Config &config)
{
  try
  {
    FileName  fileNameOut;
    KernelPtr kernel;
    Angle     minv, maxv, sampling;
    Double    R, H;

    readConfig(config, "outputfileMatrix", fileNameOut, Config::MUSTSET,  "", "matrix with first column the angle [degree], second the kernel value");
    readConfig(config, "kernel",           kernel,      Config::MUSTSET,  "", "");
    readConfig(config, "minAngle",         minv,        Config::DEFAULT,  "-180", "[degree]");
    readConfig(config, "maxAngle",         maxv,        Config::DEFAULT,  "180",  "[degree]");
    readConfig(config, "sampling",         sampling,    Config::DEFAULT,  "1",    "[degree]");
    readConfig(config, "height",           H,           Config::DEFAULT,  "0",    "evaluate at R+height [m]");
    readConfig(config, "R",                R,           Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    const UInt count = static_cast<UInt>((maxv-minv)/sampling+1);
    Matrix A(count, 2);
    logTimerStart;
    for(UInt i=0; i<count; i++)
    {
      logTimerLoop(i, count);
      const Double psi = std::min(Double(maxv), minv + i*Double(sampling));
      A(i,0) = psi*RAD2DEG;
      A(i,1) = kernel->kernel(Vector3d((R+H)*sin(psi), 0, (R+H)*cos(psi)), Vector3d(0, 0, R));
    }
    logTimerLoopEnd(count);

    logStatus<<"writing output file <"<<fileNameOut<<">"<<Log::endl;
    writeFileMatrix(fileNameOut, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
