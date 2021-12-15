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

GnssDesignMatrix::GnssDesignMatrix(const GnssNormalEquationInfo &normalEquationInfo_, const_MatrixSliceRef l_) :
      normalEquationInfo(normalEquationInfo_),
      blockIndices(normalEquationInfo.blockIndices()),
      indexUsedParameter(normalEquationInfo.blockCount()),
      countUsedParameter(normalEquationInfo.blockCount()),
      row(0),
      rows(l_.rows()),
      A(l_.rows(), normalEquationInfo.parameterCount()),
      l(l_)
{
}

/***********************************************/

void GnssDesignMatrix::init(const_MatrixSliceRef l)
{
  try
  {
    this->l = l;
    if(A.rows() < l.rows())
      A = Matrix(l.rows(), blockIndices.back());
    row  = 0;
    rows = l.rows();

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
    if(!rows)
      rows = l.rows()-row_;
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

void GnssDesignMatrix::accumulateNormals(MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount)
{
  try
  {
    for(UInt ii=0; ii<indexUsedBlock.size(); ii++)
    {
      const UInt blocki = indexUsedBlock[ii];
      normals.setBlock(blocki, blocki);
      if(!normals.N(blocki, blocki).size())
        normals.N(blocki, blocki) = Matrix(normals.blockSize(blocki), Matrix::SYMMETRIC);

      for(UInt kk=ii+1; kk<indexUsedBlock.size(); kk++)
      {
        const UInt blockk = indexUsedBlock[kk];
        normals.setBlock(blocki, blockk);
        if(!normals.N(blocki, blockk).size())
          normals.N(blocki, blockk) = Matrix(normals.blockSize(blocki), normals.blockSize(blockk));
      }

      for(UInt i=0; i<indexUsedParameter[blocki].size(); i++)
      {
        const UInt index = indexUsedParameter[blocki][i];
        const UInt count = countUsedParameter[blocki][i];
        const const_MatrixSlice Ai(A.slice(row, blockIndices[blocki]+index, rows, count).trans());

        // right hand side
        matMult(1., Ai, l.row(row, rows), n.at(blocki).row(index, count));

        // diagonal block
        rankKUpdate(1., Ai.trans(), normals.N(blocki, blocki).slice(index, index, count, count));
        for(UInt k=i+1; k<indexUsedParameter[blocki].size(); k++)
          matMult(1., Ai, A.slice(row, blockIndices[blocki]+indexUsedParameter[blocki][k], rows, countUsedParameter[blocki][k]),
                  normals.N(blocki, blocki).slice(index, indexUsedParameter[blocki][k], count, countUsedParameter[blocki][k]));

        // other blocks
        for(UInt kk=ii+1; kk<indexUsedBlock.size(); kk++)
        {
          const UInt blockk = indexUsedBlock[kk];
          for(UInt k=0; k<indexUsedParameter[blockk].size(); k++)
            matMult(1., Ai, A.slice(row, blockIndices[blockk]+indexUsedParameter[blockk][k], rows, countUsedParameter[blockk][k]),
                    normals.N(blocki, blockk).slice(index, indexUsedParameter[blockk][k], count, countUsedParameter[blockk][k]));
        }
      }
    }

    // accumulate right hand side
    obsCount += rows;
    lPl      += quadsum(l);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
