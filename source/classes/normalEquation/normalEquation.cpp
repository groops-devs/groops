/***********************************************/
/**
* @file normalEquation.cpp
*
* @brief Systems of normal equations.
* Creation, combination and solving of systems of normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2004-12-10
*
*/
/***********************************************/

#define DOCSTRING_NormalEquation

#include "base/import.h"
#include "config/configRegister.h"
#include "parallel/parallel.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"
#include "misc/varianceComponentEstimation.h"
#include "normalEquationDesign.h"
#include "normalEquationDesignVCE.h"
#include "normalEquationFile.h"
#include "normalEquationRegularization.h"
#include "normalEquationRegularizationGeneralized.h"
#include "normalEquation.h"

/***********************************************/

GROOPS_REGISTER_CLASS(NormalEquation, "normalEquationType",
                      NormalEquationDesign,
                      NormalEquationDesignVCE,
                      NormalEquationFile,
                      NormalEquationRegularization,
                      NormalEquationRegularizationGeneralized)

GROOPS_RENAMED_CLASS(normalequationType, normalEquationType, date2time(2020, 6, 3))

GROOPS_READCONFIG_UNBOUNDED_CLASS(NormalEquation, "normalEquationType")

/***********************************************/

NormalEquation::NormalEquation(Config &config, const std::string &name)
{
  try
  {
    std::string choice;
    while(readConfigChoice(config, name, choice, Config::OPTIONAL, "", "system of normal equations"))
    {
      if(readConfigChoiceElement(config, "design",                    choice, "from observation equations"))
        normalsComponent.push_back(new NormalEquationDesign(config));
      if(readConfigChoiceElement(config, "designVCE",                 choice, "from observation equations with variance component estimation per arc"))
        normalsComponent.push_back(new NormalEquationDesignVCE(config));
      if(readConfigChoiceElement(config, "file",                      choice, "read normal equations from file"))
        normalsComponent.push_back(new NormalEquationFile(config));
      if(readConfigChoiceElement(config, "regularization",            choice, "diagonal matrix"))
        normalsComponent.push_back(new NormalEquationRegularization(config));
      if(readConfigChoiceElement(config, "regularizationGeneralized", choice, "regularization with a sum of partial covariance matrices"))
        normalsComponent.push_back(new NormalEquationRegularizationGeneralized(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    };

    status = UNKNOWN;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalEquation::init(UInt blockSize, Parallel::CommunicatorPtr comm)
{
  try
  {
    status = UNKNOWN;

    // determine dimension
    UInt paraCount = 0;
    UInt rhsCount  = 0;
    for(auto component : normalsComponent)
    {
      paraCount = std::max(paraCount, component->parameterCount());
      rhsCount  = std::max(rhsCount,  component->rightHandSideCount());
    }
    if((paraCount == 0) || (rhsCount == 0))
      throw(Exception("Cannot determine dimension of normals"));

    // init distributed normal matrix
    normals.initEmpty(MatrixDistributed::computeBlockIndex(paraCount, blockSize), comm);
    n        = Matrix(paraCount, rhsCount);
    lPl      = Vector(rhsCount);
    obsCount = 0;
    x        = Matrix(paraCount, rhsCount);
    Wz       = Matrix(paraCount, 100);

    // init normals components
    for(auto component : normalsComponent)
      component->init(normals, rhsCount);

    status = INIT;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

NormalEquation::~NormalEquation()
{
  for(auto component : normalsComponent)
    delete component;
}

/***********************************************/

std::vector<ParameterName> NormalEquation::parameterNames() const
{
  try
  {
    std::vector<ParameterName> names(parameterCount());
    for(auto component : normalsComponent)
      component->parameterNames(names);
    return names;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector NormalEquation::varianceComponentFactors() const
{
  std::vector<Double> s;
  for(auto component : normalsComponent)
  {
    const std::vector<Double> s2 = component->varianceComponentFactors();
    s.insert(s.end(), s2.begin(), s2.end());
  }
  return s;
}

/***********************************************/

void NormalEquation::setApproximateSolution(const const_MatrixSlice &x0)
{
  try
  {
    if(status==UNKNOWN)
      throw(Exception("NormalEquation is not initialized"));

    x = x0;
    Parallel::broadCast(x, 0, normals.communicator());

    if((x.rows() != n.rows()) || (x.columns() != n.columns()))
      throw(Exception("dimension error"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool NormalEquation::build(UInt rightHandSide)
{
  try
  {
    if(status==UNKNOWN)
      throw(Exception("NormalEquation is not initialized"));

    rhsNo = rightHandSide;
    Parallel::broadCast(rhsNo, 0, normals.communicator());

    // accumulate normals
    normals.setNull();
    n.setNull();
    lPl.setNull();
    obsCount = 0;

    Bool ready = TRUE;
    for(auto component : normalsComponent)
      ready = component->addNormalEquation(rhsNo, x, Wz, normals, n, lPl, obsCount) && ready;

    status = NORMAL;
    return ready || (varianceComponentFactors().rows() == 1); // if only one sigma -> iteration is not needed
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalEquation::write(const FileName &fileName)
{
  try
  {
    if(status!=NORMAL)
      build(rhsNo);

    NormalEquationInfo info(parameterNames(), lPl, observationCount());
    writeFileNormalEquation(fileName, info, normals, n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix NormalEquation::solve()
{
  try
  {
    if(status!=NORMAL)
      build(rhsNo);

    // regularize not used parameters
    for(UInt i=0; i<normals.blockCount(); i++)
      if(normals.isMyRank(i,i))
      {
        Matrix &N = normals.N(i,i);
        for(UInt k=0; k<N.rows(); k++)
          if(N(k,k) == 0.)
          {
            N(k,k) += 1.0;
            logWarning<<normals.blockIndex(i)+k<<". parameter has zero diagonal element -> set to one"<<Log::endl;
          }
      }

    x = normals.solve(n, TRUE/*timing*/);
    Parallel::broadCast(x, 0, normals.communicator());

    // N contains now the cholesky decomposition
    Wz = Vce::monteCarlo(x.rows(), 100);
    normals.triangularSolve(Wz);
    Parallel::broadCast(Wz, 0, normals.communicator());

    status = CHOLESKY;
    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double NormalEquation::aposterioriSigma()
{
  try
  {
    Double s;
    if(Parallel::isMaster(normals.communicator()))
      s = std::sqrt((lPl(rhsNo) - inner(x.column(rhsNo), n.column(rhsNo))) / (observationCount() - x.rows()));
    Parallel::broadCast(s, 0, normals.communicator());
    return s;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector NormalEquation::sigmaParameter()
{
  try
  {
    if(status != CHOLESKY)
      solve();

    Double sigma = aposterioriSigma();
    if((sigma <= 0) || std::isnan(sigma))
    {
      logWarningOnce<<"sigma = "<<sigma<<" not applied to covariance matrix"<<Log::endl;
      sigma = 1.;
    }

    normals.choleskyInverse();
    Vector diagonal = Vector(normals.dimension());
    for(UInt i=0; i<normals.blockCount(); i++)
    {
      if(normals.isMyRank(i,i))
      {
        Matrix &N = normals.N(i,i);
        for(UInt z=0; z<N.rows(); z++)
          diagonal(normals.blockIndex(i)+z) += quadsum(N.slice(z,z,1,N.columns()-z));
      }
      for(UInt k=i+1; k<normals.blockCount(); k++)
        if(normals.isMyRank(i,k))
        {
          Matrix &N = normals.N(i,k);
          for(UInt z=0; z<N.rows(); z++)
            diagonal(normals.blockIndex(i)+z) += quadsum(N.row(z));
        }
    }
    Parallel::reduceSum(diagonal, 0, normals.communicator());
    Parallel::broadCast(diagonal, 0, normals.communicator());
    for(UInt i=0; i<diagonal.rows(); i++)
      diagonal(i) = std::sqrt(diagonal(i));

    status = INVERSECHOLESKY;
    return sigma * diagonal;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalEquation::writeCovariance(const FileName &fileName)
{
  try
  {
    if(status != INVERSE)
    {
      if(status != INVERSECHOLESKY)
        sigmaParameter();
      normals.choleskyProduct();
    }

    Double sigma = aposterioriSigma();
    if((sigma <= 0) || std::isnan(sigma))
    {
      logWarningOnce<<"sigma = "<<sigma<<" not applied to covariance matrix"<<Log::endl;
      sigma = 1.;
    }

    for(UInt i=0; i<normals.blockCount(); i++)
      for(UInt k=i; k<normals.blockCount(); k++)
        if(normals.isMyRank(i,k))
          writeFileMatrix(fileName.appendBaseName((normals.blockCount()>1) ? "."+i%"%02i-"s+k%"%02i"s : ""s), sigma*sigma * normals.N(i,k));

    status = INVERSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix NormalEquation::contribution()
{
  try
  {
    if(status != INVERSE)
    {
      if(status != INVERSECHOLESKY)
        sigmaParameter();
      normals.choleskyProduct();
      status = INVERSE;
    }

    Matrix x(normals.dimension(), normalsComponent.size());
    for(UInt i=0; i<normalsComponent.size(); i++)
      copy(normalsComponent.at(i)->contribution(normals), x.column(i));
    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
