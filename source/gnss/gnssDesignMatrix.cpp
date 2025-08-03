/***********************************************/
/**
* @file gnssDesignMatrix.cpp
*
* @brief Management of sparse design matrix.
*
* @author Torsten Mayer-Guerr
* @date 2018-04-01
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "gnss/gnssNormalEquationInfo.h"
#include "gnssDesignMatrix.h"

/***********************************************/

GnssDesignMatrix::GnssDesignMatrix(const GnssNormalEquationInfo &normalEquationInfo_, UInt rows_) :
      normalEquationInfo(normalEquationInfo_),
      blockIndices(normalEquationInfo.blockIndices()),
      indexUsedParameter(normalEquationInfo.blockCount()),
      countUsedParameter(normalEquationInfo.blockCount()),
      row(0),
      rows(rows_),
      A(rows_, normalEquationInfo.parameterCount())
{
}

/***********************************************/

void GnssDesignMatrix::init(UInt rows_)
{
  try
  {
    if(A.rows() < rows_)
      A = Matrix(rows_, blockIndices.back());
    row  = 0;
    rows = rows_;

    for(UInt i : indexUsedBlock)
    {
      for(UInt k=0; k<indexUsedParameter[i].size(); k++)
        A.column(blockIndices[i]+indexUsedParameter[i][k], countUsedParameter[i][k]).setNull();
      indexUsedParameter[i].clear();
      countUsedParameter[i].clear();
    }

    indexUsedBlock.clear();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssDesignMatrix &GnssDesignMatrix::selectRows(UInt row_, UInt rows_)
{
  try
  {
    row  = row_;
    rows = rows_;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

MatrixSlice GnssDesignMatrix::column(const GnssParameterIndex &index)
{
  try
  {
    const UInt block = normalEquationInfo.block(index);
    const UInt col   = normalEquationInfo.index(index) - normalEquationInfo.blockIndex(block);
    const UInt cols  = normalEquationInfo.count(index);
    return column(block, col, cols);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

MatrixSlice GnssDesignMatrix::column(UInt block, UInt col, UInt cols)
{
  try
  {
    const auto iter = std::lower_bound(indexUsedBlock.begin(), indexUsedBlock.end(), block);
    if(iter == indexUsedBlock.end() || block < *iter)
      indexUsedBlock.insert(iter, block);

    constexpr UInt gap = 32;
    UInt idx = 0;
    while((idx < indexUsedParameter[block].size()) && (indexUsedParameter[block][idx]+countUsedParameter[block][idx]+gap < col))
      idx++;
    if((idx >= indexUsedParameter[block].size()) || (col+cols+gap < indexUsedParameter[block][idx]))
    {
      indexUsedParameter[block].insert(indexUsedParameter[block].begin()+idx, col);
      countUsedParameter[block].insert(countUsedParameter[block].begin()+idx, cols);
    }
    else // merge
    {
      const UInt end = std::max(col+cols, indexUsedParameter[block][idx]+countUsedParameter[block][idx]);
      indexUsedParameter[block][idx] = std::min(col, indexUsedParameter[block][idx]);
      countUsedParameter[block][idx] = end - indexUsedParameter[block][idx];
      // merge with following parameter group?
      while((idx+1 < indexUsedParameter[block].size()) && (indexUsedParameter[block][idx]+countUsedParameter[block][idx]+gap >= indexUsedParameter[block][idx+1]))
      {
        countUsedParameter[block][idx] = std::max(countUsedParameter[block][idx], indexUsedParameter[block][idx+1]+countUsedParameter[block][idx+1]-indexUsedParameter[block][idx]);
        indexUsedParameter[block].erase(indexUsedParameter[block].begin()+idx+1);
        countUsedParameter[block].erase(countUsedParameter[block].begin()+idx+1);
      }
    }

    return A.slice(row, blockIndices[block]+col, rows, cols);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix GnssDesignMatrix::mult(const_MatrixSliceRef x)
{
  try
  {
    Matrix y(rows, x.columns());

    for(UInt block : indexUsedBlock)
      for(UInt k=0; k<indexUsedParameter[block].size(); k++)
      {
        const UInt index = blockIndices[block]+indexUsedParameter[block][k];
        const UInt count = countUsedParameter[block][k];
        matMult(1., A.slice(row, index, rows, count), x.row(index, count), y);
      }

    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix GnssDesignMatrix::mult(const std::vector<Matrix> &x, UInt startBlock, UInt countBlock)
{
  try
  {
    Matrix y(rows, x.at(startBlock).columns());
    for(UInt block : indexUsedBlock)
      if((startBlock <= block) && (block < startBlock+countBlock))
        for(UInt k=0; k<indexUsedParameter[block].size(); k++)
        {
          const UInt index = blockIndices[block]+indexUsedParameter[block][k];
          const UInt count = countUsedParameter[block][k];
          matMult(1., A.slice(row, index, rows, count), x.at(block).row(indexUsedParameter[block][k], count), y);
        }
    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssDesignMatrix::transMult(const_MatrixSliceRef l, std::vector<Matrix> &x, UInt startBlock, UInt countBlock)
{
  try
  {
    for(UInt block : indexUsedBlock)
      if((startBlock <= block) && (block < startBlock+countBlock))
        for(UInt k=0; k<indexUsedParameter[block].size(); k++)
        {
          const UInt index = blockIndices[block]+indexUsedParameter[block][k];
          const UInt count = countUsedParameter[block][k];
          matMult(1., A.slice(row, index, rows, count).trans(), l, x.at(block).row(indexUsedParameter[block][k], count));
        }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssDesignMatrix::accumulateNormals(const GnssDesignMatrix &A, const_MatrixSliceRef l, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount)
{
  try
  {
    for(UInt ii=0; ii<A.indexUsedBlock.size(); ii++)
    {
      const UInt blocki = A.indexUsedBlock[ii];
      normals.setBlock(blocki, blocki);
      if(!normals.N(blocki, blocki).size())
        normals.N(blocki, blocki) = Matrix(normals.blockSize(blocki), Matrix::SYMMETRIC);

      for(UInt kk=ii+1; kk<A.indexUsedBlock.size(); kk++)
      {
        const UInt blockk = A.indexUsedBlock[kk];
        normals.setBlock(blocki, blockk);
        if(!normals.N(blocki, blockk).size())
          normals.N(blocki, blockk) = Matrix(normals.blockSize(blocki), normals.blockSize(blockk));
      }

      for(UInt i=0; i<A.indexUsedParameter[blocki].size(); i++)
      {
        const UInt index = A.indexUsedParameter[blocki][i];
        const UInt count = A.countUsedParameter[blocki][i];
        const const_MatrixSlice Ai(A.A.slice(A.row, A.blockIndices[blocki]+index, A.rows, count).trans());

        // right hand side
        matMult(1., Ai, l, n.at(blocki).row(index, count));

        // diagonal block
        rankKUpdate(1., Ai.trans(), normals.N(blocki, blocki).slice(index, index, count, count));
        for(UInt k=i+1; k<A.indexUsedParameter[blocki].size(); k++)
          matMult(1., Ai, A.A.slice(A.row, A.blockIndices[blocki]+A.indexUsedParameter[blocki][k], A.rows, A.countUsedParameter[blocki][k]),
                  normals.N(blocki, blocki).slice(index, A.indexUsedParameter[blocki][k], count, A.countUsedParameter[blocki][k]));

        // other blocks
        for(UInt kk=ii+1; kk<A.indexUsedBlock.size(); kk++)
        {
          const UInt blockk = A.indexUsedBlock[kk];
          for(UInt k=0; k<A.indexUsedParameter[blockk].size(); k++)
            matMult(1., Ai, A.A.slice(A.row, A.blockIndices[blockk]+A.indexUsedParameter[blockk][k], A.rows, A.countUsedParameter[blockk][k]),
                    normals.N(blocki, blockk).slice(index, A.indexUsedParameter[blockk][k], count, A.countUsedParameter[blockk][k]));
        }
      }
    }

    // accumulate right hand side
    obsCount += l.rows();
    lPl      += quadsum(l);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssDesignMatrix::axpy(const std::vector<UInt> &rowInA, const std::vector<Double> &factors, const GnssDesignMatrix &B, GnssDesignMatrix &A)
{
  try
  {
    UInt idxUsedBlock = 0;
    for(UInt block : B.indexUsedBlock)
    {
      while((idxUsedBlock < A.indexUsedBlock.size()) && (A.indexUsedBlock[idxUsedBlock] < block))
        idxUsedBlock++;
      if((idxUsedBlock >= A.indexUsedBlock.size()) || block < A.indexUsedBlock[idxUsedBlock])
        A.indexUsedBlock.insert(A.indexUsedBlock.begin()+idxUsedBlock, block);

      UInt idx = 0;
      for(UInt k=0; k<B.indexUsedParameter[block].size(); k++)
      {
        const UInt col  = B.indexUsedParameter[block][k];
        const UInt cols = B.countUsedParameter[block][k];
        {
          while((idx < A.indexUsedParameter[block].size()) && (A.indexUsedParameter[block][idx]+A.countUsedParameter[block][idx] < col))
            idx++;
          if((idx >= A.indexUsedParameter[block].size()) || (col+cols < A.indexUsedParameter[block][idx])) // new block
          {
            A.indexUsedParameter[block].insert(A.indexUsedParameter[block].begin()+idx, col);
            A.countUsedParameter[block].insert(A.countUsedParameter[block].begin()+idx, cols);
          }
          else if((A.indexUsedParameter[block][idx] != col) || (A.countUsedParameter[block][idx] != cols)) // merge
          {
            const UInt end = std::max(col+cols, A.indexUsedParameter[block][idx]+A.countUsedParameter[block][idx]);
            A.indexUsedParameter[block][idx] = std::min(col, A.indexUsedParameter[block][idx]);
            A.countUsedParameter[block][idx] = end - A.indexUsedParameter[block][idx];
            // merge with following parameter group?
            while((idx+1 < A.indexUsedParameter[block].size()) && (A.indexUsedParameter[block][idx]+A.countUsedParameter[block][idx] >= A.indexUsedParameter[block][idx+1]))
            {
              A.countUsedParameter[block][idx] = std::max(A.countUsedParameter[block][idx], A.indexUsedParameter[block][idx+1]+A.countUsedParameter[block][idx+1]-A.indexUsedParameter[block][idx]);
              A.indexUsedParameter[block].erase(A.indexUsedParameter[block].begin()+idx+1);
              A.countUsedParameter[block].erase(A.countUsedParameter[block].begin()+idx+1);
            }
          }
        }
        MatrixSlice       As(A.A.slice(A.row, A.blockIndices[block]+col, A.rows, cols));
        const_MatrixSlice Bs(B.A.slice(B.row, B.blockIndices[block]+col, B.rows, cols));
        for(UInt i=0; i<rowInA.size(); i++)
          if(rowInA.at(i) != NULLINDEX)
            ::axpy(factors.at(i), Bs.row(i), As.row(rowInA.at(i)));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
