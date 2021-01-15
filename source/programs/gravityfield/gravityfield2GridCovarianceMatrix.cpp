/***********************************************/
/**
* @file gravityfield2GridCovarianceMatrix.cpp
*
* @brief Covariance matrix of values of a gravity field on a grid.
*
* @author Torsten Mayer-Guerr
* @date 2012-05-30
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program propagates the covariance matrix of a \configClass{gravityfield}{gravityfieldType}
evaluated at \config{time} to a \configClass{grid}{gridType}. The full variance-covariance matrix is computed
and written to a \file{matrix file}{matrix}:
\begin{equation}
\mathbf{\Sigma}_\mathbf{y} = \mathbf{F}\mathbf{\Sigma}_\mathbf{x}\mathbf{F}^T
\end{equation}
The \configClass{kernel}{kernelType} determines the quantity of the grid values, for example,
\configClass{kernel:waterHeight}{kernelType:waterHeight}.

See also \program{GravityfieldCovariancesPropagation2GriddedData}, \program{GravityfieldVariancesPropagation2GriddedData}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/grid/grid.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Covariance matrix of values of a gravity field on a grid.
* @ingroup programsGroup */
class Gravityfield2GridCovarianceMatrix
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2GridCovarianceMatrix, SINGLEPROCESS, "Covariance matrix of values of a gravity field on a grid.", Gravityfield)

/***********************************************/

void Gravityfield2GridCovarianceMatrix::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName        outName;
    GridPtr         grid;
    GravityfieldPtr gravityfield;
    KernelPtr       kernel;
    Time            time;

    readConfig(config, "outputfileMatrix",     outName,      Config::MUSTSET,  "", "symmetric grid covariance matrix");
    readConfig(config, "grid",                 grid,         Config::MUSTSET,  "", "");
    readConfig(config, "kernel",               kernel,       Config::MUSTSET,  "", "");
    readConfig(config, "gravityfield",         gravityfield, Config::MUSTSET,  "", "");
    readConfig(config, "time",                 time,         Config::OPTIONAL, "", "at this time the gravity field will be evaluated");
    if(isCreateSchema(config)) return;

    // Fcreate covariance matrix
    // ------------------
    logStatus<<"create covariance matrix"<<Log::endl;
    Matrix C = gravityfield->variance(time, grid->points(), *kernel);

    // write file
    // ----------
    logStatus<<"write covariance matrix to file <"<<outName<<">"<<Log::endl;
    writeFileMatrix(outName, C);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
