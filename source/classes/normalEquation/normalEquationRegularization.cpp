/***********************************************/
/**
* @file normalEquationRegularization.cpp
*
* @brief Regularization with diagonal matrix.
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
#include "parallel/parallel.h"
#include "classes/normalEquation/normalEquation.h"
#include "classes/normalEquation/normalEquationRegularization.h"

/***********************************************/

NormalEquationRegularization::NormalEquationRegularization(Config &config)
{
  try
  {
    FileName diagName, biasName;
    Double   sigma;

    readConfig(config, "inputfileDiagonalMatrix", diagName,   Config::OPTIONAL, "",    "Vector with the diagonal elements of the weight matrix");
    readConfig(config, "inputfileBias",           biasName,   Config::OPTIONAL, "",    "Matrix with right hand sides");
    readConfig(config, "aprioriSigma",            sigma,      Config::DEFAULT,  "1.0", "");
    readConfig(config, "startIndex",              startIndex, Config::DEFAULT,  "0",   "regularization of parameters starts at this index (counting from 0)");
    if(isCreateSchema(config)) return;

    paraCount = 0;
    obsCount  = 0;
    rhsCount  = 1;
    sigma2    = sigma*sigma;

    if(!diagName.empty())
    {
      readFileMatrix(diagName, K);
      paraCount = K.rows();
    }

    if(!biasName.empty())
    {
      readFileMatrix(biasName, bias);
      paraCount = bias.rows();
      rhsCount  = bias.columns();
    }

    if(K.rows() && bias.rows() && (K.rows() != bias.rows()))
      throw(Exception("dimension of diagonal matrix ("+K.rows()%"%i) and bias ("s+bias.rows()%"%i x "s+bias.columns()%"%i) must agree!"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalEquationRegularization::init(MatrixDistributed &normals, UInt rhsCount)
{
  try
  {
    if(paraCount == 0)
      paraCount = normals.parameterCount()-startIndex;
    this->rhsCount = rhsCount;

    // adjust right hand side
    // ----------------------
    if(!bias.size())
      bias = Matrix(paraCount, rhsCount);
    if(bias.columns() != rhsCount)
      throw(Exception("Dimension error (right hand side count): "+bias.columns()%"%i != "s+rhsCount%"%i"s));

    // adjust diagonal matrix
    // ----------------------
    if(!K.size())
      K = Vector(paraCount, 1.);

    // compute lPl
    // -----------
    lPl = Vector(rhsCount);
    for(UInt rhsNo=0; rhsNo<rhsCount; rhsNo++)
      for(UInt i=0; i<paraCount; i++)
        lPl(rhsNo) += bias(i, rhsNo) * K(i) * bias(i, rhsNo);

    // compute obsCount
    // ----------------
    obsCount = 0;
    for(UInt i=0; i<paraCount; i++)
      if(K(i) != 0.)
        obsCount++;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Bool NormalEquationRegularization::addNormalEquation(UInt rhsNo, const const_MatrixSlice &x, const const_MatrixSlice &Wz,
                                                     MatrixDistributed &normals, Matrix &n, Vector &lPl, UInt &obsCount)
{
  try
  {
    UInt ready = 0;
    if(Parallel::isMaster())
    {
      Vector vx = x.slice(startIndex, rhsNo, bias.rows(), 1);
      if(quadsum(vx) > 0.)
      {
        // square sum of residuals
        Double ePe = 0.0;
        for(UInt i=0; i<paraCount; i++)
          ePe += (vx(i)-bias(i, rhsNo)) * K(i) * (vx(i)-bias(i, rhsNo)); // = (x-n)'K(x-n)
        // redundancy (monte carlo estimation)
        Double r = static_cast<Double>(this->obsCount);
        for(UInt i=0; i<paraCount; i++)
          r -= 1/sigma2 * K(i) * quadsum(Wz.row(i+startIndex)); // = z'W'KWz
        // aposteriori sigma
        const Double sigma2Old = sigma2;
        sigma2 = ePe/r;
        ready  = (std::fabs(std::sqrt(sigma2)-std::sqrt(sigma2Old))/std::sqrt(sigma2) < 0.01);
      }
    }
    Parallel::broadCast(sigma2);
    Parallel::broadCast(ready);

    // add normals
    const Double factor = 1/sigma2;
    for(UInt i=0; i<normals.blockCount(); i++)
      if((normals.blockIndex(i+1) > startIndex) && (normals.blockIndex(i) < startIndex+paraCount))
      {
        normals.setBlock(i,i);
        if(normals.isMyRank(i,i))
          for(UInt z=0; z<normals.blockSize(i); z++)
            if((z+normals.blockIndex(i) >= startIndex) && (z+normals.blockIndex(i)-startIndex < paraCount))
              normals.N(i,i)(z,z) += factor * K(z+normals.blockIndex(i)-startIndex);
      }

    // accumulate right hand sides
    if(Parallel::isMaster())
    {
      for(UInt z=0; z<paraCount; z++)
        for(UInt s=0; s<rhsCount; s++)
          n(z+startIndex, s) += factor * K(z) * bias(z,s);
      axpy(factor, this->lPl, lPl);
      obsCount += this->obsCount;
    }

    return ready;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector NormalEquationRegularization::contribution(MatrixDistributed &Cov)
{
  try
  {
    Vector contrib(Cov.dimension());

    for(UInt i=0; i<Cov.blockCount(); i++)
      if(Cov.isMyRank(i,i))
        for(UInt z=0; z<Cov.blockSize(i); z++)
          if((z+Cov.blockIndex(i) >= startIndex) && (z+Cov.blockIndex(i) < startIndex+paraCount))
            contrib(z+Cov.blockIndex(i)) += Cov.N(i,i)(z,z) * K(z+Cov.blockIndex(i)-startIndex)/sigma2;

    Parallel::reduceSum(contrib);
    Parallel::broadCast(contrib);
    return contrib;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
