/***********************************************/
/**
* @file matrixCalculate.cpp
*
* @brief Matrix manipulation.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program creates a \file{matrix}{matrix} from multiple matrices.
All \configClass{matrices}{matrixGeneratorType} are summed up. The size of the resulting matrix is exandeded to fit all matrices.
The class \configClass{matrixGenerator}{matrixGeneratorType} allows complex matrix operations before.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/matrixGenerator/matrixGenerator.h"

/***********************************************/

/** @brief Matrix manipulation.
* @ingroup programsGroup */
class MatrixCalculate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(MatrixCalculate, SINGLEPROCESS, "Matrix manipulation", Misc, Matrix)

/***********************************************/

void MatrixCalculate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName           fileName;
    MatrixGeneratorPtr matrixGenerator;

    readConfig(config, "outputfileMatrix", fileName, Config::MUSTSET,  "", "");
    readConfig(config, "matrix", matrixGenerator, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"calculate matrix"<<Log::endl;
    Matrix A = matrixGenerator->compute();
    logInfo<<"  dimension of resulting matrix: "<<A.rows()<<" x "<<A.columns()<<Log::endl;

    logStatus<<"write matrix to <"<<fileName<<">"<<Log::endl;
    writeFileMatrix(fileName, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
