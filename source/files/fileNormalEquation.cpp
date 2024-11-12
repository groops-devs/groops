/***********************************************/
/**
* @file fileNormalEquation.cpp
*
* @brief Read/write a system of normal equations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-05-27
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_NormalEquation

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "parallel/matrixDistributed.h"
#include "files/fileFormatRegister.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"
#include "files/fileNormalEquation.h"

GROOPS_REGISTER_FILEFORMAT(NormalEquation, FILE_NORMALEQUATION_TYPE)

/***********************************************/

static void writeInfoFile(const FileName &name, const NormalEquationInfo &info, const Matrix &n)
{
  try
  {
    if((info.blockIndex.back() != n.rows()) || (info.lPl.rows() != n.columns()) || (info.parameterName.size() != n.rows()))
      throw(Exception("<"+name.str()+"> dimension error"));
    if(n.rows() && !info.usedBlocks.size())
      throw(Exception("no used blocks defined"));

    FileName nameInfo(name.replaceFullExtension(".info.xml"));
    OutFileArchive file(nameInfo, FILE_NORMALEQUATION_TYPE, FILE_NORMALEQUATION_VERSION);
    file<<nameValue("observationCount",   info.observationCount);
    file<<nameValue("parameterCount",     info.blockIndex.back());
    file<<nameValue("rightHandSideCount", info.lPl.rows());
    for(UInt i=0; i<info.lPl.rows(); i++)
      file<<nameValue("lPl", info.lPl(i));
    file<<nameValue("blockCount", info.blockIndex.size()-1);
    for(UInt i=0; i<info.blockIndex.size()-1; i++)
      file<<nameValue("blockSize", info.blockIndex.at(i+1)-info.blockIndex.at(i));
//     file<<nameValue("extension",  name.extension().str());
    file<<nameValue("usedBlocks", info.usedBlocks);

    // parameter name file
    writeFileParameterName(name.replaceFullExtension(".parameterNames.txt"), info.parameterName);

    // right hand side
    writeFileMatrix(name.appendBaseName(".rightHandSide"), n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileNormalEquation(const FileName &name, NormalEquationInfo info, const Matrix &N, const Matrix &n)
{
  try
  {
    Bool writeBlock = !isStrictlyZero(N);
    // info file
    info.blockIndex = {0, N.rows()};
    info.usedBlocks = identityMatrix(writeBlock ? 1 : 0);
    writeInfoFile(name, info, n);

    // normal matrix
    if(writeBlock)
      writeFileMatrix(name, N);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileNormalEquation(const FileName &name, NormalEquationInfo info, const std::vector<std::vector<Matrix>> &N, const Matrix &n)
{
  try
  {
    const UInt blockCount = N.size();
    info.blockIndex = std::vector<UInt>(blockCount+1, 0);
    for(UInt i=0; i<blockCount; i++)
      info.blockIndex.at(i+1) = info.blockIndex.at(i) + N.at(i).at(i).rows();

    info.usedBlocks = Matrix(blockCount, Matrix::SYMMETRIC);
    for(UInt i=0; i<blockCount; i++)
      for(UInt k=i; k<blockCount; k++)
        if(N.at(i).at(k).size() && !isStrictlyZero(N.at(i).at(k)))
          info.usedBlocks(i,k) = 1;

    // info file
    writeInfoFile(name, info, n);

    // normal matrix
    for(UInt i=0; i<blockCount; i++)
      for(UInt k=i; k<blockCount; k++)
        if(info.usedBlocks(i,k) > 0)
          writeFileMatrix(name.appendBaseName((blockCount>1) ? "."+i%"%02i-"s+k%"%02i"s : ""s), N.at(i).at(k));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileNormalEquation(const FileName &name, NormalEquationInfo info, const MatrixDistributed &normal, const Matrix &n)
{
  try
  {
    info.blockIndex = normal.blockIndex();
    info.usedBlocks = Matrix(normal.blockCount(), Matrix::SYMMETRIC);
    for(UInt i=0; i<normal.blockCount(); i++)
      for(UInt k=i; k<normal.blockCount(); k++)
        if(normal.isMyRank(i,k))
          info.usedBlocks(i,k) = isStrictlyZero(normal.N(i, k)) ? 0 : 1;
    Parallel::reduceSum(info.usedBlocks, 0, normal.communicator());

    // info file
    if(Parallel::isMaster(normal.communicator()))
      writeInfoFile(name, info, n);

    // normal matrix
    for(UInt i=0; i<normal.blockCount(); i++)
      for(UInt k=i; k<normal.blockCount(); k++)
        if(normal.isMyRank(i,k) && !isStrictlyZero(normal.N(i, k)))
          writeFileMatrix(name.appendBaseName((normal.blockCount()>1) ? "."+i%"%02i-"s+k%"%02i"s : ""s), normal.N(i, k));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileNormalEquation(const FileName &name, NormalEquationInfo info, const Matrix &n)
{
  try
  {
    writeInfoFile(name, info, n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static void readInfoFile(const FileName &name, NormalEquationInfo &info, Matrix &n)
{
  try
  {
    info.blockIndex.clear();
    info.blockIndex.push_back(0);

    FileName nameInfo(name.replaceFullExtension(".info.xml"));
    InFileArchive file(nameInfo, FILE_NORMALEQUATION_TYPE, FILE_NORMALEQUATION_VERSION);
    UInt paramCount, rhsCount, blockCount;
    file>>nameValue("observationCount",   info.observationCount);
    file>>nameValue("parameterCount",     paramCount);
    file>>nameValue("rightHandSideCount", rhsCount);
    // old version
    if(file.version() == 1)
    {
      logWarning<<"file <"<<name<<"> is saved as an old version, please save as new version"<<Log::endl;
      file>>nameValue("lPl", info.lPl);
      info.blockIndex.push_back(info.blockIndex.back() + paramCount);
      blockCount = 1;
    }
    else // new version
    {
      info.lPl = Vector(rhsCount);
      for(UInt i=0; i<rhsCount; i++)
        file>>nameValue("lPl", info.lPl(i));
      file>>nameValue("blockCount", blockCount);
      UInt size;
      for(UInt i=0; i<blockCount; i++)
      {
        file>>nameValue("blockSize", size);
        info.blockIndex.push_back(info.blockIndex.back() + size);
      }
//     file>>nameValue("extension", name.stripExtension());
    }

    try
    {
      file>>nameValue("usedBlocks", info.usedBlocks);
    }
    catch(std::exception &/*e*/)
    {
      info.usedBlocks = Matrix(blockCount, Matrix::SYMMETRIC);
      for(UInt i=0; i<blockCount; i++)
        for(UInt k=i; k<blockCount; k++)
          info.usedBlocks(i,k) = 1;
    }

    // read parameter name file
    // ------------------------
    info.parameterName.clear();
    FileName fileNameParameterNames(name.replaceFullExtension(".parameterNames.txt"));
    try
    {
      readFileParameterName(fileNameParameterNames, info.parameterName);
    }
    catch(std::exception &/*e*/)
    {
      info.parameterName.clear();
      info.parameterName.resize(paramCount);
    }

    // read ride hand side
    // ------------------------
    readFileMatrix(name.appendBaseName(".rightHandSide"), n);
    if((n.columns() != info.lPl.rows()) || (n.rows() != info.blockIndex.back()))
      throw(Exception("dimension error"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileNormalEquation(const FileName &name, NormalEquationInfo &info, Matrix &N, Matrix &n)
{
  try
  {
    readInfoFile(name, info, n);
    const std::vector<UInt> *blockIndex = &info.blockIndex; // just for readability

    // read normal equation
    if(blockIndex->size() == 2)
    {
      readFileMatrix(name, N);
      if(N.rows() != blockIndex->back())
        throw(Exception("<"+name.str()+"> dimension error"));
    }
    else
    {
      N = Matrix(blockIndex->back(), Matrix::SYMMETRIC);
      for(UInt i=0; i<blockIndex->size()-1; i++)
        for(UInt k=i; k<blockIndex->size()-1; k++)
          if(!info.usedBlocks.size() || (info.usedBlocks.size() && info.usedBlocks(i,k) > 0))
          {
            Matrix M;
            readFileMatrix(name.appendBaseName("."+i%"%02i-"s+k%"%02i"s), M);
            copy(M, N.slice(blockIndex->at(i), blockIndex->at(k), blockIndex->at(i+1)-blockIndex->at(i), blockIndex->at(k+1)-blockIndex->at(k)));
          }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileNormalEquation(const FileName &name, NormalEquationInfo &info, MatrixDistributed &normal, Matrix &n, Parallel::CommunicatorPtr comm)
{
  try
  {
    readInfoFile(name, info, n);

    // read normal equation
    normal.initEmpty(info.blockIndex, comm);
    for(UInt i=0; i<normal.blockCount(); i++)
      for(UInt k=i; k<normal.blockCount(); k++)
      {
        if((!info.usedBlocks.size() || (info.usedBlocks.size() && info.usedBlocks(i,k) > 0)))
        {
          normal.setBlock(i, k);
          if(normal.isMyRank(i,k))
            readFileMatrix(name.appendBaseName((normal.blockCount()>1) ? "."+i%"%02i-"s+k%"%02i"s : ""s), normal.N(i,k));
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileNormalEquation(const FileName &name, NormalEquationInfo &info, Matrix &n)
{
  try
  {
    readInfoFile(name, info, n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
