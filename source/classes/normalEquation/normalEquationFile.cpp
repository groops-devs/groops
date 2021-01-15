/***********************************************/
/**
* @file normalEquationFile.cpp
*
* @brief Read normals from file.
* @see NormalEquation
*
* @author Torsten Mayer-Guerr
* @date 2004-12-10
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"
#include "parallel/parallel.h"
#include "parallel/matrixDistributed.h"
#include "classes/normalEquation/normalEquation.h"
#include "classes/normalEquation/normalEquationFile.h"

/***********************************************/

NormalEquationFile::NormalEquationFile(Config &config)
{
  try
  {
    Double sigma;

    renameDeprecatedConfig(config, "inputfileNormalequation", "inputfileNormalEquation", date2time(2020, 6, 3));

    readConfig(config, "inputfileNormalEquation", fileName,   Config::MUSTSET,  "",    "");
    readConfig(config, "aprioriSigma",            sigma,      Config::DEFAULT,  "1.0", "");
    readConfig(config, "startIndex",              startIndex, Config::DEFAULT,  "0",   "add this normals at index of total matrix (counting from 0)");
    if(isCreateSchema(config)) return;

    NormalEquationInfo info;
    Matrix             n;
    readFileNormalEquation(fileName, info, n);
    lPl       = info.lPl;
    rhsCount  = n.columns();
    paraCount = n.rows();
    obsCount  = info.observationCount;
    names     = info.parameterName;
    sigma2    = sigma * sigma;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalEquationFile::parameterNames(std::vector<ParameterName> &names) const
{
  try
  {
    for(UInt i=0; i<this->names.size(); i++)
      if(!names.at(i+startIndex).combine(this->names.at(i)))
        logWarningOnce<<"Parameter names do not match at index "<<i+startIndex<<": '"<<names.at(i+startIndex).str()<<"' != '"<<this->names.at(i).str()<<"'"<< Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalEquationFile::init(MatrixDistributed &normalsTotal, UInt rhsCount)
{
  try
  {
    NormalEquationInfo info;
    readFileNormalEquation(fileName, info, normals, n, normalsTotal.communicator());

    std::vector<UInt> index(normalsTotal.parameterCount(), NULLINDEX);
    std::iota(index.begin()+startIndex, index.begin()+startIndex+normals.parameterCount(), 0);
    normals.reorder(index, normalsTotal.blockIndex(), normals.getCalculateRank());

    // init right hand sides
    // ---------------------
    if(Parallel::isMaster(normalsTotal.communicator()))
    {
      if(rhsCount != n.columns())
        throw(Exception("number of right hand sides must agree ("+rhsCount%"%i != "s+n.columns()%"%i)"s));
      if((startIndex > 0) || (normalsTotal.parameterCount() > n.rows()))
      {
        Matrix tmp(normalsTotal.parameterCount(), rhsCount);
        copy(n, tmp.row(startIndex, n.rows()));
        std::swap(tmp, n);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool NormalEquationFile::addNormalEquation(UInt rhsNo, const const_MatrixSlice &x, const const_MatrixSlice &Wz,
                                           MatrixDistributed &normalsTotal, Matrix &nTotal, Vector &lPlTotal, UInt &obsCountTotal)
{
  try
  {
    UInt ready = 0;
    const Vector vx = x.column(rhsNo);
    if(!isStrictlyZero(vx))
    {
      Matrix Nvx(vx.rows(), vx.columns());
      Matrix NWz(Wz.rows(), Wz.columns());
      for(UInt i=0; i<normals.blockCount(); i++)
        for(UInt k=i; k<normals.blockCount(); k++)
          if(normals.isMyRank(i, k))
          {
            matMult(1., normals.N(i,k), vx.row(normals.blockIndex(k), normals.blockSize(k)), Nvx.row(normals.blockIndex(i), normals.blockSize(i)));
            matMult(1., normals.N(i,k), Wz.row(normals.blockIndex(k), normals.blockSize(k)), NWz.row(normals.blockIndex(i), normals.blockSize(i)));
            if(i == k)
              continue;
            matMult(1., normals.N(i,k).trans(), vx.row(normals.blockIndex(i), normals.blockSize(i)), Nvx.row(normals.blockIndex(k), normals.blockSize(k)));
            matMult(1., normals.N(i,k).trans(), Wz.row(normals.blockIndex(i), normals.blockSize(i)), NWz.row(normals.blockIndex(k), normals.blockSize(k)));
          }
      Parallel::reduceSum(Nvx, 0, normalsTotal.communicator());
      Parallel::reduceSum(NWz, 0, normalsTotal.communicator());

      if(Parallel::isMaster(normalsTotal.communicator()) )
      {
        // aposteriori sigma
        const Double sigma2Old = sigma2;
        const Double ePe = inner(vx, Nvx) - 2*inner(vx, n.column(rhsNo)) + lPl(rhsNo);
        const Double r   = obsCount - inner(Wz, NWz)/sigma2;
        sigma2 = ePe/r;
        ready  = (std::fabs(std::sqrt(sigma2)-std::sqrt(sigma2Old))/std::sqrt(sigma2) < 0.01);
      }
      Parallel::broadCast(sigma2, 0, normalsTotal.communicator());
      Parallel::broadCast(ready,  0, normalsTotal.communicator());
    }

    // accumulate normals
    // ------------------
    for(UInt i=0; i<normals.blockCount(); i++)
      for(UInt k=i; k<normals.blockCount(); k++)
        if(normals.isBlockUsed(i, k))
        {
          normalsTotal.setBlock(i, k);
          if(normals.isMyRank(i, k))
            axpy(1./sigma2, normals.N(i,k), normalsTotal.N(i,k));
        }

    // accumulate right hand sides
    // ---------------------------
    if(Parallel::isMaster(normalsTotal.communicator()))
    {
      axpy(1./sigma2, n,   nTotal);
      axpy(1./sigma2, lPl, lPlTotal);
      obsCountTotal += obsCount;
    }

    return ready;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector NormalEquationFile::contribution(MatrixDistributed &Cov)
{
  try
  {
    Vector contrib(Cov.dimension());
    for(UInt i=0; i<normals.blockCount(); i++)
    {
      // diagonal
      if(normals.isMyRank(i, i))
      {
        fillSymmetric(normals.N(i,i));
        normals.N(i,i).setType(Matrix::GENERAL);
        Matrix NN = Cov.N(i,i) * normals.N(i,i);
        for(UInt z=0; z<NN.rows(); z++)
          contrib(z+Cov.blockIndex(i)) += NN(z,z);
        normals.N(i,i).setType(Matrix::SYMMETRIC);
      }

      for(UInt k=i+1; k<normals.blockCount(); k++)
        if(normals.isMyRank(i, k))
        {
          Matrix NN = Cov.N(i,k) * normals.N(i,k).trans();
          for(UInt z=0; z<NN.rows(); z++)
            contrib(z+Cov.blockIndex(i)) += NN(z,z);
          NN = Cov.N(i,k).trans() * normals.N(i,k);
          for(UInt z=0; z<NN.rows(); z++)
            contrib(z+Cov.blockIndex(k)) += NN(z,z);
        }
    }

    Parallel::reduceSum(contrib, 0, Cov.communicator());
    Parallel::broadCast(contrib, 0, Cov.communicator());
    return contrib;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
