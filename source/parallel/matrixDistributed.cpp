/***********************************************/
/**
* @file matrixDistributed.cpp
*
* @brief Positve definte matrix in distributed memory.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2011-01-30
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "matrixDistributed.h"

/***********************************************/

MatrixDistributed::MatrixDistributed() : _blockIndex{0}
{
  setCalculateRank(nullptr);
}

/***********************************************/

MatrixDistributed::MatrixDistributed(const std::vector<UInt> &blockIndex, Parallel::CommunicatorPtr comm, std::function<UInt(UInt, UInt, UInt)> calcRank)
{
  init(blockIndex, comm, calcRank);
}

/***********************************************/

void MatrixDistributed::init(const std::vector<UInt> &blockIndex, Parallel::CommunicatorPtr comm, std::function<UInt(UInt, UInt, UInt)> calcRank)
{
  try
  {
    initEmpty(blockIndex, comm, calcRank);

    for(UInt i=0; i<blockCount(); i++)
      if(blockSize(i))
        for(UInt k=i; k<blockCount(); k++)
          if(blockSize(k))
            setBlock(i, k);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::initEmpty(const std::vector<UInt> &blockIndex, Parallel::CommunicatorPtr comm, std::function<UInt(UInt, UInt, UInt)> calcRank)
{
  try
  {
    if(!std::is_sorted(blockIndex.begin(), blockIndex.end()))
      throw(Exception("Block indices must be given in ascending order."));

    if(blockIndex.front() != static_cast<UInt>(0))
      throw(Exception("Block index must start with zero."));

    this->comm = comm;
    _blockIndex = blockIndex;
    _row.clear();    _row.resize(blockCount());
    _column.clear(); _column.resize(blockCount());
    _rank.clear();
    _N.clear();
    setCalculateRank(calcRank);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

UInt MatrixDistributed::index(UInt row, UInt col) const
{
  if((row > blockCount()) || (col > blockCount()))
    throw(Exception("Access to block ("+row%"%i"s+" x "+col%"%i"s+") out of range, matrix size ("+blockCount()%"%i"s+" x "+blockCount()%"%i"s+")."));

  const auto iter = std::lower_bound(_column[row].begin(), _column[row].end(), col, [](const auto &x, UInt col) {return x.first < col;});
  if((iter != _column[row].end()) && (iter->first == col))
    return iter->second;
  return NULLINDEX;
}

/***********************************************/

void MatrixDistributed::loopBlockColumn(const std::array<UInt,2> &rows, UInt col, std::function<void(UInt, UInt)> block) const
{
  try
  {
    const auto &_row = this->_row.at(col);
    for(UInt idx=0; idx<_row.size(); idx++)
      if((rows[0] <= _row[idx].first) && (_row[idx].first < rows[1]))
        block(_row[idx].first, _row[idx].second);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::loopBlockRow(UInt row, const std::array<UInt,2> &cols, std::function<void(UInt, UInt)> block) const
{
  try
  {
    const auto &_column = this->_column.at(row);
    for(UInt idx=0; idx<_column.size(); idx++)
      if((cols[0] <= _column[idx].first) && (_column[idx].first < cols[1]))
        block(_column[idx].first, _column[idx].second);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::vector<Bool> MatrixDistributed::usedRanksInColumn(const std::array<UInt,2> &rows, UInt col) const
{
  try
  {
    std::vector<Bool> usedRank(Parallel::size(comm), FALSE);
    loopBlockColumn(rows, col, [&](UInt /*row*/, UInt idx) {usedRank.at(_rank[idx]) = TRUE;});
    return usedRank;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Bool> MatrixDistributed::usedRanksInRow(UInt row, const std::array<UInt,2> &cols) const
{
  try
  {
    std::vector<Bool> usedRank(Parallel::size(comm), FALSE);
    loopBlockRow(row, cols, [&](UInt /*col*/, UInt idx) {usedRank.at(_rank[idx]) = TRUE;});
    return usedRank;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::broadCast(Matrix &x, UInt idx, const std::vector<Bool> &usedRank)
{
  try
  {
    std::vector<UInt> ranks = {_rank[idx]};
    for(UInt idProcess=0; idProcess<usedRank.size(); idProcess++)
      if(usedRank.at(idProcess) && (idProcess != _rank[idx]))
        ranks.push_back(idProcess);
    Parallel::CommunicatorPtr commNew = Parallel::createCommunicator(ranks, comm);
    if(commNew)
      Parallel::broadCast(x, 0, commNew);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::reduceSum(Matrix &x, UInt idx, const std::vector<Bool> &usedRank, Bool free)
{
  try
  {
    std::vector<UInt> ranks = {_rank[idx]};
    for(UInt idProcess=0; idProcess<usedRank.size(); idProcess++)
      if(usedRank.at(idProcess) && (idProcess != _rank[idx]))
        ranks.push_back(idProcess);
    Parallel::CommunicatorPtr commNew = Parallel::createCommunicator(ranks, comm);
    if(!commNew)
      return;
    Parallel::reduceSum(x, 0, commNew);
    if(!isMyRank(idx))
    {
      if(free)
        x = Matrix();
      else
        x.setNull();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

UInt MatrixDistributed::rank(UInt i, UInt k) const
{
  const UInt ik = index(i, k);
  if(ik == NULLINDEX)
    return NULLINDEX;
  return _rank[ik];
}

/***********************************************/

UInt MatrixDistributed::index2block(UInt i) const
{
  return std::distance(_blockIndex.begin(), std::upper_bound(_blockIndex.begin(), _blockIndex.end(), i))-1;
}

/***********************************************/

Matrix &MatrixDistributed::N(UInt i, UInt k)
{
  const UInt ik = index(i, k);
  if(ik == NULLINDEX)
    throw(Exception("In MatrixDistributed::N("+i%"%i, "s+k%"%i): block not exist"s));
  return _N[ik];
}

/***********************************************/

const Matrix &MatrixDistributed::N(UInt i, UInt k) const
{
  const UInt ik = index(i, k);
  if(ik == NULLINDEX)
    throw(Exception("In MatrixDistributed::N("+i%"%i, "s+k%"%i): block not exist"s));
  return _N[ik];
}

/***********************************************/
/***********************************************/

void MatrixDistributed::setCalculateRank(std::function<UInt(UInt, UInt, UInt)> calcRank_)
{
  try
  {
    calcRank = calcRank_;
    if(calcRank == nullptr)
      calcRank = calculateRankBlockCyclic;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt MatrixDistributed::calculateRankBlockCyclic(UInt i, UInt k, UInt commSize)
{
  try
  {
    // find optimal process grid (nearly quadratic)
    UInt pRows = static_cast<UInt>(std::floor(std::sqrt(commSize)));
    while(commSize % pRows)
      pRows++;
    const UInt pCols = commSize/pRows;

    return (i%pRows)*pCols+(k%pCols);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt MatrixDistributed::setBlock(UInt i, UInt k, UInt rank)
{
  try
  {
    if(rank == NULLINDEX)
      rank = this->rank(i, k); // rank of used block
    if(rank == NULLINDEX)
      rank = calcRank(i, k, Parallel::size(comm)); // new block? => calculate new rank

    UInt ik = index(i, k);
    if(ik == NULLINDEX)
    {
      ik = _N.size();
      _N.push_back(Matrix());
      _rank.push_back(rank);

      _column[i].insert(std::lower_bound(_column[i].begin(), _column[i].end(), k, [](const auto &x, UInt k) {return x.first < k;}), std::pair<UInt, UInt>(k, ik));
      _row[k].insert   (std::lower_bound(_row[k].begin(),   _row[k].end(),     i, [](const auto &x, UInt i) {return x.first < i;}), std::pair<UInt, UInt>(i, ik));
    }

    _rank[ik] = rank;
    if(isMyRank(ik) && (_N[ik].size() == 0))
      _N[ik] = ((i==k) ? Matrix(blockSize(i), Matrix::SYMMETRIC) : Matrix(blockSize(i), blockSize(k)));
    return ik;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::eraseBlocks(UInt start, UInt count)
{
  try
  {
    if(start+count > blockCount())
      throw(Exception("Range of block rows/columns to remove ("+start%"%i"s+" to "+(start+count-1)%"%i"s+") exceeds the number of block rows/columns in the matrix ("+blockCount()%"%i"s+")."));

    std::vector<std::vector<std::pair<UInt, UInt>>> _rowNew(blockCount()-count);      // each column, used row -> idx to _N and _rank
    std::vector<std::vector<std::pair<UInt, UInt>>> _columnNew(blockCount()-count);   // each row, used column -> idx to _N and _rank
    std::vector<Matrix>               _NNew;                            // unorderd list of used blocks
    std::vector<UInt>                 _rankNew;                         // unorderd list of rank of used blocks
    for(UInt i=0; i<blockCount(); i++)
      if((i < start) || (i >= start+count))
        loopBlockRow(i, {0, blockCount()}, [&](UInt k, UInt ik)
        {
          if((k < start) || (k >= start+count))
          {
            _columnNew[(i < start) ? i : i-count].push_back(std::pair<UInt, UInt>((k < start) ? k : k-count, _NNew.size()));
            _rowNew   [(k < start) ? k : k-count].push_back(std::pair<UInt, UInt>((i < start) ? i : i-count, _NNew.size()));
            _NNew.push_back(_N[ik]);
            _rankNew.push_back(_rank[ik]);
          }
        });

    _row    = _rowNew;
    _column = _columnNew;
    _N      = _NNew;
    _rank   = _rankNew;

    // adjust block indices
    const UInt erasedParameterCount = blockIndex(start+count) - blockIndex(start);
    for(UInt i=start; i<blockCount(); i++)
      _blockIndex[i+1] -= erasedParameterCount;
    _blockIndex.erase(_blockIndex.begin() + start, _blockIndex.begin() + start+count);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void MatrixDistributed::setNull()
{
  try
  {
    for(UInt i=0; i<blockCount(); i++)
      loopBlockRow(i, {i, blockCount()}, [&](UInt k, UInt ik)
      {
        if(isMyRank(ik))
        {
          if(i == k)
            _N[ik].setType(Matrix::SYMMETRIC, Matrix::UPPER);
          _N[ik].setNull();
        }
      });
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::reduceSum(UInt i, UInt k)
{
  try
  {
    if(Parallel::size(comm)<=1)
      return;

    const UInt ik = index(i, k);
    if(ik == NULLINDEX)
      throw(Exception("N("+i%"%i, "s+k%"%i): block not exist"s));
    const UInt color = _N[ik].size()  ? k : NULLINDEX;
    const UInt key   = isMyRank(ik) ? 0 : Parallel::myRank(comm)+1;
    Parallel::CommunicatorPtr commNew = Parallel::splitCommunicator(color, key, comm);
    if(commNew && (Parallel::size(commNew)>1))
      Parallel::reduceSum(_N[ik], 0, commNew);
    if(!isMyRank(ik))
      _N[ik] = Matrix();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::reduceSum(Bool timing)
{
  try
  {
    if(Parallel::size(comm)<=1)
      return;

    if(timing) logTimerStart;
    UInt idxBlock = 0;
    for(UInt i=0; i<blockCount(); i++)
      loopBlockRow(i, {i, blockCount()}, [&](UInt /*k*/, UInt ik)
      {
        if(timing) logTimerLoop(idxBlock++, _N.size());
        UInt color = _N[ik].size()  ? idxBlock : NULLINDEX;
        UInt key   = isMyRank(ik) ? 0 : Parallel::myRank(comm)+1;
        Parallel::CommunicatorPtr commNew = Parallel::splitCommunicator(color, key, comm);
        if(commNew && (Parallel::size(commNew)>1))
          Parallel::reduceSum(_N[ik], 0, commNew);
        if(!isMyRank(ik))
          _N[ik] = Matrix();
      });
    Parallel::barrier(comm);
    if(timing) logTimerLoopEnd(_N.size());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

void MatrixDistributed::cholesky(Bool timing, UInt startBlock, UInt countBlock, Bool collect)
{
  try
  {
    if(timing) logTimerStart;
    for(UInt i=startBlock; i<blockCount(); i++)
      if(blockSize(i))
      {
        if(timing) logTimerLoop(i-startBlock, blockCount()-startBlock);
        UInt ii = index(i,i);
        if((ii == NULLINDEX) && (i < startBlock+countBlock))
          throw(Exception("Diagonal block ("+i%"%i, "s+i%"%i) is not set."s));

        loopBlockColumn({startBlock, std::min(i, startBlock+countBlock)}, i, [&](UInt z, UInt zi)
        {
          if(ii == NULLINDEX)
            ii = setBlock(i,i);

          // column rank k update
          if(isMyRank(zi))
          {
            if(_N[ii].size() == 0)
              _N[ii] = Matrix(blockSize(i), Matrix::SYMMETRIC, Matrix::UPPER);
            rankKUpdate(-1., _N[zi], _N[ii]);
          }

          // distribute top column to right hand side blocks
          if(Parallel::size(comm) > 1)
            broadCast(_N[zi], zi, usedRanksInRow(z, {i+1, blockCount()}));

          // dgemm
          loopBlockRow(z, {i+1, blockCount()}, [&](UInt s, UInt zs)
          {
            const UInt is = setBlock(i, s);
            if(isMyRank(zs))
            {
              if(_N[is].size() == 0)
                _N[is] = Matrix(blockSize(i), blockSize(s));
              matMult(-1., _N[zi].trans(), _N[zs], _N[is]);
            }
          });

          // free column
          if(!isMyRank(zi) && _N[zi].size())
            _N[zi] = Matrix();
        }); // for(row z)

        // collect right row elements from top block
        if((i>0) && (Parallel::size(comm) > 1) && ((i < startBlock+countBlock) || collect))
          loopBlockRow(i, {i, blockCount()}, [&](UInt s, UInt is)
          {
            std::vector<Bool> usedRank(Parallel::size(comm), FALSE);
            loopBlockColumn({startBlock, std::min(i, startBlock+countBlock)}, i, [&](UInt z, UInt /*zi*/)
            {
              const UInt zs = index(z,s);
              if(zs != NULLINDEX)
                usedRank.at(_rank[zs]) = TRUE;
            });
            reduceSum(_N[is], is, usedRank);
          });

        if(i < startBlock+countBlock)
        {
          // cholesky
          if(isMyRank(ii))
            ::cholesky(_N[ii]);

          // distribute diagonal element to row
          if(Parallel::size(comm) > 1)
            broadCast(_N[ii], ii, usedRanksInRow(i, {i+1, blockCount()}));

          // triangularSolve to row
          loopBlockRow(i, {i+1, blockCount()}, [&](UInt /*s*/, UInt is)
          {
            if(isMyRank(is))
              ::triangularSolve(1., _N[ii].trans(), _N[is]);
          });

          // free diagonal
          if(ii != NULLINDEX && !isMyRank(ii) && _N[ii].size())
            _N[ii] = Matrix();
        }
      }
    Parallel::barrier(comm);
    if(timing) logTimerLoopEnd(blockCount()-startBlock);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix MatrixDistributed::solve(const_MatrixSliceRef n, Bool timing)
{
  try
  {
    cholesky(timing);

    UInt rhsCount = n.columns();
    Parallel::broadCast(rhsCount, 0, comm);

    std::vector<Matrix> x(blockCount());
    if(Parallel::isMaster(comm))
      for(UInt i=0; i<blockCount(); i++)
        x.at(i) = n.row(blockIndex(i), blockSize(i));
    else
      for(UInt i=0; i<blockCount(); i++)
        x.at(i) = Matrix(blockSize(i), rhsCount);

    triangularTransSolve(x);
    triangularSolve(x);

    Matrix x2(n.rows(), n.columns());
    if(Parallel::isMaster(comm))
      for(UInt i=0; i<blockCount(); i++)
        copy(x.at(i), x2.row(blockIndex(i), blockSize(i)));

    return x2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::triangularSolve(MatrixSliceRef x2)
{
  try
  {
    UInt rhsCount = x2.columns();
    Parallel::broadCast(rhsCount, 0, comm);

    std::vector<Matrix> x(blockCount());
    if(Parallel::isMaster(comm))
      for(UInt i=0; i<blockCount(); i++)
        x.at(i) = x2.row(blockIndex(i), blockSize(i));
    else
      for(UInt i=0; i<blockCount(); i++)
        x.at(i) = Matrix(blockSize(i), rhsCount);

    triangularSolve(x);

    if(Parallel::isMaster(comm))
      for(UInt i=0; i<blockCount(); i++)
        copy(x.at(i), x2.row(blockIndex(i), blockSize(i)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::triangularSolve(std::vector<Matrix> &x, UInt startBlock, UInt countBlock)
{
  try
  {
    for(UInt i=startBlock+countBlock; i-->startBlock;)
      if(blockSize(i))
      {
        const UInt ii = index(i,i);

        // collect
        if(Parallel::size(comm) > 1)
        {
          std::vector<Bool> usedRank = usedRanksInRow(i, {i, startBlock+countBlock});
          usedRank.at(0) = TRUE; // master
          reduceSum(x.at(i), ii, usedRank);
        }

        // solve
        if(isMyRank(ii))
          ::triangularSolve(1., _N[ii], x.at(i));

        // distribute to top column
        if(Parallel::size(comm) > 1)
        {
          std::vector<Bool> usedRank = usedRanksInColumn({startBlock, i}, i);
          usedRank.at(0) = TRUE; // master
          broadCast(x.at(i), ii, usedRank);
        }

        // reduce
        loopBlockColumn({startBlock, i}, i, [&](UInt z, UInt zi)
        {
          if(isMyRank(zi))
            matMult(-1., _N[zi], x.at(i), x.at(z));
        });

        // free
        if(!Parallel::isMaster(comm))
          x.at(i).setNull();
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::triangularTransSolve(MatrixSliceRef x2, UInt startBlock, UInt countBlock)
{
  try
  {
    UInt rhsCount = x2.columns();
    Parallel::broadCast(rhsCount, 0, comm);

    std::vector<Matrix> x(blockCount());
    if(Parallel::isMaster(comm))
      for(UInt i=0; i<blockCount(); i++)
        x.at(i) = x2.row(blockIndex(i), blockSize(i));
    else
      for(UInt i=0; i<blockCount(); i++)
        x.at(i) = Matrix(blockSize(i), rhsCount);

    triangularTransSolve(x, startBlock, countBlock, TRUE);

    if(Parallel::isMaster(comm))
      for(UInt i=0; i<blockCount(); i++)
        copy(x.at(i), x2.row(blockIndex(i), blockSize(i)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::triangularTransSolve(std::vector<Matrix> &x, UInt startBlock, UInt countBlock, Bool collect)
{
  try
  {
    for(UInt i=startBlock; i<startBlock+countBlock; i++)
      if(blockSize(i))
      {
        const UInt ii = index(i,i);

        // collect
        if(Parallel::size(comm) > 1)
        {
          std::vector<Bool> usedRank = usedRanksInColumn({startBlock, i}, i);
          usedRank.at(0) = TRUE; // master
          reduceSum(x.at(i), ii, usedRank, collect);
        }

        // solve
        if(isMyRank(ii))
          ::triangularSolve(1., _N[ii].trans(), x.at(i));

        // distribute to row
        if(Parallel::size(comm) > 1)
        {
          std::vector<Bool> usedRank = usedRanksInRow(i, {i, blockCount()});
          usedRank.at(0) = TRUE; // master
          broadCast(x.at(i), ii, usedRank);
        }

        // reduce
        loopBlockRow(i, {i+1, blockCount()}, [&](UInt s, UInt is)
        {
          if(isMyRank(is) && _N[is].size())
            matMult(-1., _N[is].trans(), x.at(i), x.at(s));
        });

        // free
        if(!Parallel::isMaster(comm))
          x.at(i).setNull();
      }

    // reduce special block
    // --------------------
    if(collect && (Parallel::size(comm) > 1))
      for(UInt i=startBlock+countBlock; i<blockCount(); i++)
        if(blockSize(i))
        {
          std::vector<UInt> ranks = {0};
          std::vector<Bool> usedRank(Parallel::size(comm), FALSE);
          usedRank.at(0) = TRUE; // master
          loopBlockColumn({startBlock, startBlock+countBlock}, i, [&](UInt /*z*/, UInt zi)
          {
            if(!usedRank.at(_rank[zi]))
            {
              usedRank.at(_rank[zi]) = TRUE;
              ranks.push_back(_rank[zi]);
            }
          });
          Parallel::CommunicatorPtr commNew = Parallel::createCommunicator(ranks, comm);
          if(commNew)
            Parallel::reduceSum(x.at(i), 0, commNew);
          // free
          if(!Parallel::isMaster(comm))
            x.at(i).setNull();
        }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::choleskyInverse(Bool timing, UInt startBlock, UInt countBlock)
{
  try
  {
    if(timing) logTimerStart;
    for(UInt i=startBlock; i<startBlock+countBlock; i++)
      if(blockSize(i))
      {
        if(timing) logTimerLoop(i-startBlock, countBlock);
        const UInt ii = index(i,i);

        // distribute top column elements to left triangular
        if(Parallel::size(comm) > 1)
          loopBlockColumn({startBlock, i}, i, [&](UInt z, UInt zi)
          {
            broadCast(_N[zi], zi, usedRanksInColumn({0, z+1}, z));
          });

        // compute triangularMult
        for(UInt z=startBlock; z<i; z++)
        {
          const UInt zz = index(z,z);
          const UInt zi = index(z,i);
          if(zi != NULLINDEX)
          {
            if((zz != NULLINDEX) && isMyRank(zz))
              triangularMult(1., _N[zz], _N[zi]);
            else if(isMyRank(zi))
              _N[zi].setNull();
            else if(_N[zi].size())
              _N[zi] = Matrix();
          }

          loopBlockRow(z, {z+1, i}, [&](UInt s, UInt zs)
          {
            if(isMyRank(zs))
            {
              const UInt si = index(s,i);
              if(si != NULLINDEX)
              {
                const UInt zi = setBlock(z, i);
                if(_N[zi].size() == 0)
                  _N[zi] = Matrix(blockSize(z), blockSize(i));
                matMult(1.,  _N[zs],  _N[si], _N[zi]);
              }
            }
          });
        }

        // reduceSum top column elements from left triangular
        if(Parallel::size(comm) > 1)
          loopBlockColumn({startBlock, i}, i, [&](UInt z, UInt zi)
          {
            std::vector<Bool> usedRank(Parallel::size(comm), FALSE);
            loopBlockRow(z, {z, i}, [&](UInt s, UInt zs)
            {
              if(isBlockUsed(s, i))
                usedRank.at(_rank[zs]) = TRUE;
            });
            if(usedRank.at(Parallel::myRank(comm)) && !_N[zi].size())
              _N[zi] = Matrix(blockSize(z), blockSize(i));
            reduceSum(_N[zi], zi, usedRank);
          });

        // distribute diagonal element to top column
        if((i>startBlock) && (Parallel::size(comm) > 1))
          broadCast(_N[ii], ii, usedRanksInColumn({startBlock, i+1}, i));

        // triangularSolve to column
        loopBlockColumn({startBlock, i}, i, [&](UInt /*z*/, UInt zi)
        {
          if(isMyRank(zi))
            ::triangularSolve(-1., _N[ii].trans(), _N[zi].trans());
        });

        // free diagonal
        if((!isMyRank(ii)) && _N[ii].size())
          _N[ii] = Matrix();

        // inverte triangular element
        if(isMyRank(ii))
          ::inverse(_N[ii]);
      }
    Parallel::barrier(comm);
    if(timing) logTimerLoopEnd(countBlock);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void MatrixDistributed::choleskyProduct(Bool timing)
{
  try
  {
    if(timing) logTimerStart;
    for(UInt i=0; i<blockCount(); i++)
      if(blockSize(i))
      {
        if(timing) logTimerLoop(i, blockCount());
        const UInt ii = index(i,i);

        // distribute diagonal element to top column
        if((i>0) && (Parallel::size(comm) > 1))
          broadCast(_N[ii], ii, usedRanksInColumn({0, i+1}, i));

        // compute triangularMult
        loopBlockColumn({0, i}, i, [&](UInt /*z*/, UInt zi)
        {
          if(isMyRank(zi))
            triangularMult(1., _N[ii], _N[zi].trans());
        });

        // W'W
        if(isMyRank(ii))
          ::choleskyProduct(_N[ii]);
        else if(_N[ii].size())
          _N[ii] = Matrix(); // free diagonal

        // distribute right row elements to top block
        if((i>0) && (Parallel::size(comm) > 1))
          loopBlockRow(i, {i+1, blockCount()}, [&](UInt s, UInt is)
          {
            broadCast(_N[is], is, usedRanksInColumn({0, i+1}, s));
          });

        // dgemm
        for(UInt z=0; z<i; z++)
          loopBlockRow(i, {i+1, blockCount()}, [&](UInt s, UInt is)
          {
            const UInt zs = index(z,s);
            if((zs != NULLINDEX) && isMyRank(zs))
            {
              const UInt zi = setBlock(z, i);
              if(!_N[zi].size())
                _N[zi] = Matrix(blockSize(z), blockSize(i));
              matMult(1., _N[zs], _N[is].trans(), _N[zi]);
            }
          });

        // free row elements
        loopBlockRow(i, {i+1, blockCount()}, [&](UInt /*s*/, UInt is)
        {
          if(!isMyRank(is) && _N[is].size())
            _N[is] = Matrix();
        });

        // reduceSum column elements from right hand side block
        if(Parallel::size(comm) > 1)
          loopBlockColumn({0, i}, i, [&](UInt z, UInt zi)
          {
            std::vector<Bool> usedRank(Parallel::size(comm), FALSE);
            loopBlockRow(z, {i+1, blockCount()}, [&](UInt s, UInt zs)
            {
              if(isBlockUsed(i, s))
                usedRank.at(_rank[zs]) = TRUE;
            });
            reduceSum(_N[zi], zi, usedRank);
          });

        // row rank k update
        loopBlockRow(i, {i+1, blockCount()}, [&](UInt /*s*/, UInt is)
        {
          if(isMyRank(is))
          {
            if(_N[ii].size() == 0)
              _N[ii] = Matrix(blockSize(i), Matrix::SYMMETRIC, Matrix::UPPER);
            rankKUpdate(1., _N[is].trans(), _N[ii]);
          }
        });

        // reduceSum diagonal element from row
        if(Parallel::size(comm) > 1)
          reduceSum(_N[ii], ii, usedRanksInRow(i, {i+1, blockCount()}));

        // free row elements
        loopBlockRow(i, {i+1, blockCount()}, [&](UInt /*s*/, UInt is)
        {
          if((!isMyRank(is)) && _N[is].size())
            _N[is] = Matrix();
        });
      }
    Parallel::barrier(comm);
    if(timing) logTimerLoopEnd(blockCount());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void MatrixDistributed::cholesky2SparseInverse(Bool timing)
{
  try
  {
    if(timing) logTimerStart;
    for(UInt i=blockCount(); i-->0;)
      if(blockSize(i))
      {
        if(timing) logTimerLoop(blockCount()-1-i, blockCount());

        const UInt ii = index(i, i);

        std::vector<Bool> rowRanks = usedRanksInRow(i, {i, blockCount()});
        broadCast(_N[ii], ii, rowRanks); // broadcast diagonal(i,i) to whole row;

        // update off diagonal blocks in row with inverse diagonal
        // S_{i,j} = -W_{i,i}^{-1} * W_{i, i+1:n} * S_{i+1:n, j}
        loopBlockRow(i, {i+1, blockCount()}, [&](UInt /*k*/, UInt ik)
        {
          if(isMyRank(ik))
            ::triangularSolve(-1.0, _N[ii], _N[ik]);
        });

        if(isMyRank(ii))
          cholesky2Inverse(_N[ii]);
        else if(_N[ii].size())
          _N[ii].setNull();

        // distribute row to lower symmetric
        if(Parallel::size(comm) > 1)
          loopBlockRow(i, {i+1, blockCount()}, [&](UInt j, UInt ij)
          {
            std::vector<Bool> usedRank(Parallel::size(comm), FALSE);
            loopBlockRow(i, {i+1, blockCount()}, [&](UInt k, UInt /*ik*/)
            {
              const UInt jk = index(std::min(j,k), std::max(j,k)); // upper triangle of symm.
              if(jk != NULLINDEX)
                usedRank.at(_rank[jk]) = TRUE;
            });
            broadCast(_N[ij], ij, usedRank);
          });

        // symm. matMult
        std::vector<Matrix> Uij(blockCount());
        loopBlockRow(i, {i+1, blockCount()}, [&](UInt j, UInt ij) // loop over columns
        {
          std::vector<Bool> usedRank(Parallel::size(comm), FALSE);
          loopBlockRow(i, {i+1, blockCount()}, [&](UInt k, UInt ik) // loop over rows
          {
            const UInt jk = index(std::min(j,k), std::max(j,k)); // upper triangle of symm.
            if(jk != NULLINDEX)
            {
              usedRank.at(_rank[jk]) = TRUE;
              if(isMyRank(jk))
              {
                if(!Uij.at(j).size())
                  Uij.at(j) = Matrix(blockSize(i), blockSize(j));
                matMult(1., _N[ik], ((k>j) ? _N[jk].trans() : _N[jk]), Uij.at(j));
              }
            }
          });
          if(isMyRank(ij) && !Uij.at(j).size())
            Uij.at(j) = Matrix(blockSize(i), blockSize(j));
          reduceSum(Uij.at(j), ij, usedRank);
        });

        // free row elements
        if(Parallel::size(comm) > 1)
          loopBlockRow(i, {i+1, blockCount()}, [&](UInt /*j*/, UInt ij)
          {
            if((!isMyRank(ij)) && _N[ij].size())
              _N[ij] = Matrix();
          });

        // update diagonal from row
        // S_{i,i} = W_{i,i}^-1 W_{i,i}^-T - S_{i, i+1:n}*[ W_{i,i}^-1 W_{i,i+1:n}]
        loopBlockRow(i, {i+1, blockCount()}, [&](UInt j, UInt ij) // loop over columns
        {
          if(isMyRank(ij))
          {
            _N[ii].setType(Matrix::GENERAL);
            matMult(1.0, Uij.at(j), _N[ij].trans(), _N[ii]);
            _N[ij] = Uij.at(j);
          }
        });

        reduceSum(_N[ii], ii, rowRanks);
        _N[ii].setType(Matrix::SYMMETRIC);
      }
    if(timing) logTimerLoopEnd(blockCount());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void MatrixDistributed::reorder(const std::vector<UInt> &index, const std::vector<UInt> &blockIndexNew, std::function<UInt(UInt, UInt, UInt)> calcRank)
{
  try
  {
    if(index.size() != blockIndexNew.back())
      throw(Exception("index and blockIndex do not match."));

    MatrixDistributed matrixNew;
    matrixNew.initEmpty(blockIndexNew, comm, calcRank);

    // create new sort index {blockOld, indexOld, blockNew, indexNew}
    std::vector<std::array<UInt,4>> idx;
    idx.reserve(index.size());
    for(UInt i=0; i<index.size(); i++)
      if(index.at(i) != NULLINDEX)
      {
        const UInt blockOld = index2block(index.at(i));
        const UInt blockNew = matrixNew.index2block(i);
        idx.push_back({blockOld, index.at(i)-blockIndex(blockOld), blockNew, i-matrixNew.blockIndex(blockNew)});
      }
    std::sort(idx.begin(), idx.end(), [](const std::array<UInt,4> &a, const std::array<UInt,4> &b) {return (a[0] == b[0]) ? (a[1] < b[1]) : (a[0] < b[0]);});

    // ----------------------------
    auto continuousRange = [&](UInt start)
    {
      UInt i = start;
      while((i < idx.size()) &&
            (idx[i][0] == idx[start][0]) &&
            (idx[i][1] == idx[start][1]+i-start) &&
            (idx[i][2] == idx[start][2]) &&
            (idx[i][3] == idx[start][3]+i-start))
        i++;
      return i-start;
    };
    // ----------------------------

    UInt idxRowStart = 0;
    UInt idxRowEnd   = 0;
    for(UInt i=0; i<blockCount(); i++)   // loop over all block rows
    {
      idxRowStart = std::distance(idx.begin(), std::lower_bound(idx.begin()+idxRowEnd,   idx.end(), i, [](const std::array<UInt,4> &a, UInt i){return a[0] < i;}));
      idxRowEnd   = std::distance(idx.begin(), std::upper_bound(idx.begin()+idxRowStart, idx.end(), i, [](UInt i, const std::array<UInt,4> &a){return i < a[0];}));

      loopBlockRow(i, {i, blockCount()}, [&](UInt k, UInt ik)
      {
        const UInt idxColStart = std::distance(idx.begin(), std::lower_bound(idx.begin()+idxRowStart, idx.end(), k, [](const std::array<UInt,4> &a, UInt k){return a[0] < k;}));
        const UInt idxColEnd   = std::distance(idx.begin(), std::upper_bound(idx.begin()+idxColStart, idx.end(), k, [](UInt k, const std::array<UInt,4> &a){return k < a[0];}));

        // distribute block
        std::vector<Bool> usedRank(Parallel::size(comm), FALSE);
        for(UInt z=idxRowStart; z<idxRowEnd; z += continuousRange(z))
          for(UInt s=((i==k) ? z : idxColStart); s<idxColEnd; s += continuousRange(s))
          {
            UInt iNew = idx.at(z)[2];
            UInt kNew = idx.at(s)[2];
            if(kNew < iNew)
              std::swap(iNew, kNew);
            const UInt ikNew = matrixNew.setBlock(iNew, kNew);
            usedRank.at(matrixNew._rank[ikNew]) = TRUE;
          }
        broadCast(_N[ik], ik, usedRank);

        // copy elements
        UInt rows, cols;
        if(usedRank.at(Parallel::myRank(comm)))
          for(UInt z=idxRowStart; z<idxRowEnd; z+=rows)
          {
            UInt iNew = idx.at(z)[2];
            UInt row  = idx.at(z)[3];
            rows = continuousRange(z);
            for(UInt s=((i==k) ? z : idxColStart); s<idxColEnd; s += cols)
            {
              UInt kNew = idx.at(s)[2];
              UInt col  = idx.at(s)[3];
              cols = continuousRange(s);
              if((iNew > kNew) || ((iNew == kNew) && (row > col))) // transpose to access upper triangle?
              {
                const UInt kiNew = matrixNew.index(kNew, iNew);
                if(matrixNew.isMyRank(kiNew))
                  copy(_N[ik].slice(idx.at(z)[1], idx.at(s)[1], rows, cols), matrixNew._N[kiNew].trans().slice(row, col, rows, cols));
              }
              else
              {
                const UInt ikNew = matrixNew.index(iNew, kNew);
                if(matrixNew.isMyRank(ikNew))
                  copy(_N[ik].slice(idx.at(z)[1], idx.at(s)[1], rows, cols), matrixNew._N[ikNew].slice(row, col, rows, cols));
              }
            }
          }

        _N[ik] = Matrix();
      }); // for(block[ik])
    } // for(block row i)

    *this = matrixNew;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<UInt> MatrixDistributed::computeBlockIndex(UInt parameterCount, UInt blockSize)
{
  if(parameterCount==0)
    return {0};

  if(blockSize==0)
    return {0, parameterCount};

  std::vector<UInt> blockIndex(1, 0);
  while(blockIndex.back()<parameterCount)
    blockIndex.push_back(std::min(blockIndex.back()+blockSize, parameterCount));

  return blockIndex;
}

/***********************************************/
