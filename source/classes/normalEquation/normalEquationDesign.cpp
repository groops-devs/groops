/***********************************************/
/**
* @file normalEquationDesign.cpp
*
* @brief Accumulate normals from observation equations.
* @f[ N = A^TPA,\quad n=A^TPl @f]
* @see NormalEquation
* @see Observation
*
* @author Torsten Mayer-Guerr
* @date 2004-12-10
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "parallel/parallel.h"
#include "files/fileArcList.h"
#include "classes/observation/observation.h"
#include "classes/normalEquation/normalEquation.h"
#include "classes/normalEquation/normalEquationDesign.h"

/***********************************************/

NormalEquationDesign::NormalEquationDesign(Config &config)
{
  try
  {
    FileName fileNameArcList;

    renameDeprecatedConfig(config, "arcList", "inputfileArcList", date2time(2020, 7, 7));

    readConfig(config, "observation",      observation,     Config::MUSTSET,  "",    "");
    readConfig(config, "aprioriSigma",     sigma2,          Config::DEFAULT,  "1.0", "");
    readConfig(config, "startIndex",       startIndex,      Config::DEFAULT,  "0",   "add this normals at index of total matrix (counting from 0)");
    readConfig(config, "inputfileArcList", fileNameArcList, Config::OPTIONAL, "",    "to accelerate computation");
    if(isCreateSchema(config)) return;


    intervals = {0, observation->arcCount()};
    if(!fileNameArcList.empty())
    {
      std::vector<Time> timesInterval;
      readFileArcList(fileNameArcList, intervals, timesInterval);
    }

    sigma2   *= sigma2;
    sigma2New = sigma2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalEquationDesign::parameterNames(std::vector<ParameterName> &names) const
{
  try
  {
    std::vector<ParameterName> baseNames;
    observation->parameterName(baseNames);

    for(UInt i=0; i<baseNames.size(); i++)
      if(!names.at(i+startIndex).combine(baseNames.at(i)))
        if(Parallel::isMaster())
          logWarning<<"Parameter names do not match at index "<<i+startIndex<<": '"<<names.at(i+startIndex).str()<<"' != '"<<baseNames.at(i).str()<<"'"<< Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalEquationDesign::init(MatrixDistributed &/*normals*/, UInt rhsCount)
{
  try
  {
    if(rhsCount != observation->rightSideCount())
      throw(Exception("number of right hand sides must agree ("+rhsCount%"%i != "s+observation->rightSideCount()%"%i)"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool NormalEquationDesign::addNormalEquation(UInt rhsNo, const const_MatrixSlice &x, const const_MatrixSlice &Wz,
                                             MatrixDistributed &normals, Matrix &n, Vector &lPl, UInt &obsCount)
{
  try
  {
    if(!observation->parameterCount())
      return TRUE;


    const UInt blockStart = normals.index2block(startIndex);
    const UInt blockEnd   = normals.index2block(startIndex+observation->parameterCount()-1);
    for(UInt i=blockStart; i<=blockEnd; i++)
      for(UInt k=i; k<=blockEnd; k++)
      {
        normals.setBlock(i, k);
        normals.N(i,k) = (i==k) ? Matrix(normals.blockSize(i), Matrix::SYMMETRIC) : Matrix(normals.blockSize(i), normals.blockSize(k));
      }

    // compute observation equations
    // -----------------------------
    sigma2   = sigma2New;
    obsCount = 0;
    Double      ePe        = 0;
    Double      redundancy = 0;
    MatrixSlice x0(x.slice(startIndex, rhsNo, observation->parameterCount(), 1));
    MatrixSlice Wz0(Wz.row(startIndex, observation->parameterCount()));

    logStatus<<"accumulate normals from observation equations"<<Log::endl;
    Parallel::forEachInterval(observation->arcCount(), intervals, [&](UInt arcNo)
    {
      // observation equations
      Matrix l, A, B;
      observation->observation(arcNo, l, A, B);
      if(l.rows()==0)
        return;

      // if equations are orthogonal transformed
      // additional residuals appended to l
      Matrix l2;
      if(l.rows()>A.rows())
      {
        l2 = l.row(A.rows(), l.rows()-A.rows());
        l  = l.row(0, A.rows());
      }

      // eliminate arc related parameters
      if(B.size())
        eliminationParameter(B,A,l);

      // right hand side
      // ---------------
      matMult(1/sigma2, A.trans(), l, n);
      for(UInt i=0; i<l.columns(); i++)
        lPl(i) += quadsum(l.column(i)) + quadsum(l2.column(i))/sigma2;
      obsCount   += l.rows() + l2.rows();
      ePe        += (quadsum(l.column(rhsNo) - A*x0) + quadsum(l2.column(rhsNo)))/sigma2;
      redundancy += l.rows() - quadsum(A*Wz0)/sigma2;

      // accumulate normals
      // ------------------
      for(UInt i=blockStart; i<=blockEnd; i++)
      {
        const UInt idxN1 = (normals.blockIndex(i) < startIndex) ? (startIndex-normals.blockIndex(i)) : 0;
        const UInt idxA1 = (normals.blockIndex(i) < startIndex) ? 0 : (normals.blockIndex(i)-startIndex);
        const UInt cols1 = std::min(normals.blockSize(i)-idxN1, A.columns()-idxA1);
        rankKUpdate(1/sigma2, A.column(idxA1, cols1), normals.N(i,i).slice(idxN1, idxN1, cols1, cols1));
        for(UInt k=i+1; k<=blockEnd; k++)
        {
          const UInt idxN2 = (normals.blockIndex(k) < startIndex) ? (startIndex-normals.blockIndex(k)) : 0;
          const UInt idxA2 = (normals.blockIndex(k) < startIndex) ? 0 : (normals.blockIndex(k)-startIndex);
          const UInt cols2 = std::min(normals.blockSize(k)-idxN1, A.columns()-idxA1);
          matMult(1/sigma2, A.column(idxA1, cols1).trans(), A.column(idxA2, cols2), normals.N(i,k).slice(idxN1, idxN2, cols1, cols2));
        }
      }
    });

    normals.reduceSum(FALSE);
    Parallel::reduceSum(n);
    Parallel::reduceSum(lPl);
    Parallel::reduceSum(obsCount);

    UInt ready = 0;
    if(quadsum(x0) > 0)
    {
      Parallel::reduceSum(ePe);
      Parallel::reduceSum(redundancy);
      sigma2New = ePe/redundancy;
      ready = (std::fabs(sqrt(sigma2New)-std::sqrt(sigma2))/std::sqrt(sigma2New) < 0.01);
      Parallel::broadCast(sigma2New);
      Parallel::broadCast(ready);
    }

    return ready;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector NormalEquationDesign::contribution(MatrixDistributed &Cov)
{
  try
  {
    if(Parallel::isMaster())
      logWarning<<"In NormalEquationDesign: contribution is not implemented"<<Log::endl;
    return Vector(Cov.dimension());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
