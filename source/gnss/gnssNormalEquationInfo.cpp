/***********************************************/
/**
* @file gnssNormalEquationInfo.cpp
*
* @brief GNSS normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2010-08-03
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "parallel/matrixDistributed.h"
#include "base/parameterName.h"
#include "gnssNormalEquationInfo.h"

/***********************************************/

GnssNormalEquationInfo::GnssNormalEquationInfo(UInt countEpoch, UInt countReceiver, UInt countTransmitter, Parallel::CommunicatorPtr comm_) :
    comm(comm_),
    isEachReceiverSeparately(FALSE),
    estimateReceiver(countReceiver, TRUE),
    idEpochs(countEpoch),
    defaultBlockSizeEpoch(0),
    defaultBlockSizeInterval(64),
    defaultBlockSizeAmbiguity(64),
    defaultBlockReceiverCount(0),
    defaultBlockCountReduction(32),
    keepEpochNormalsInMemory(TRUE),
    accumulateEpochObservations(FALSE),
    blockCountEpoch_(countEpoch, 0),
    countTransmitter_(countTransmitter)
{
  std::iota(idEpochs.begin(), idEpochs.end(), 0);
}

/***********************************************/

void GnssNormalEquationInfo::initNewParameterNames()
{
  try
  {
    parameters.clear();
    block_.clear();
    index_.clear();
    count_.clear();
    parameterNames_.clear();
    blockIndices_.clear();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssParameterIndex GnssNormalEquationInfo::addParameters(UInt idEpoch, UInt idRecv, UInt idTrans, const std::vector<ParameterName> &parameterNames)
{
  if(!parameterNames.size())
    return GnssParameterIndex(NULLINDEX);
  const UInt idx = parameters.size();
  parameters.push_back(Parameter{idEpoch, idRecv, idTrans, NULLINDEX, idx, parameterNames});
  return GnssParameterIndex(parameters.back().idx);
}

/***********************************************/

void GnssNormalEquationInfo::calculateIndex(const Vector &recvProcess)
{
  try
  {
    // ------------------
    auto newBlock = [&]()
    {
      const UInt parameterCount = this->parameterCount();
      if(blockSize(blockCount()-1))
      {
        blockIndices_.push_back(parameterCount);
        blockRank_.push_back(NULLINDEX);
      }
    };
    // ------------------

    // ------------------
    auto insert = [&](auto iter, UInt defaultBlockSize)
    {
      if((defaultBlockSize > 0) && (blockSize(blockCount()-1) >= defaultBlockSize))
        newBlock();
      const UInt idx   = iter->idx;
      index_.at(idx)   = blockIndices_.back();
      block_.at(idx)   = blockCount()-1;
      count_.at(idx)   = iter->names.size();
      parameterNames_.insert(parameterNames_.end(), iter->names.begin(), iter->names.end());
      blockIndices_.back() += iter->names.size();
      blockRank_.back() = iter->rank;
    };
    // ------------------

    block_.resize(parameters.size(), NULLINDEX);
    index_.resize(parameters.size(), NULLINDEX);
    count_.resize(parameters.size(), 0);
    blockIndices_ = {0, 0};
    blockRank_    = {NULLINDEX};
    parameterNames_.reserve(std::accumulate(parameters.begin(), parameters.end(), UInt(0), [](UInt count, const auto &p){return count+p.names.size();}));

    // set process rank of receivers
    for(auto &p : parameters)
      if((p.idRecv != NULLINDEX) && (p.idTrans == NULLINDEX) && recvProcess(p.idRecv))
        p.rank = recvProcess(p.idRecv)-1;

    parameters.sort([](auto &p1, auto &p2)
                    {
                      const Bool isAmbi1 = (p1.idRecv != NULLINDEX) && (p1.idTrans != NULLINDEX);
                      const Bool isAmbi2 = (p2.idRecv != NULLINDEX) && (p2.idTrans != NULLINDEX);
                      if(isAmbi1 != isAmbi2) return isAmbi2;                         // ambiguities always at end
                      if(p1.idEpoch != p2.idEpoch) return (p1.idEpoch < p2.idEpoch); // epoch
                      if(p1.rank    != p2.rank)    return (p1.rank    < p2.rank);    // process rank
                      if(p1.idRecv  != p2.idRecv)  return (p1.idRecv  < p2.idRecv);  // idRecv
                      return (p1.idTrans < p2.idTrans);                              // idTrans
                    });
    auto iter = parameters.begin();

    // epoch parameters
    std::fill(blockCountEpoch_.begin(), blockCountEpoch_.end(), 0);
    for(UInt idEpoch : idEpochs)
    {
      newBlock();
      UInt blockEpochStart = blockIndices_.size();
      UInt idRecv = NULLINDEX;
      UInt countStation = 0;
      while((iter != parameters.end()) && (iter->idEpoch == idEpoch) && ((iter->idRecv == NULLINDEX) || (iter->idTrans == NULLINDEX)))
      {
        if(iter->idRecv != idRecv) // next receiver?
          if(defaultBlockReceiverCount && ((iter->idRecv == NULLINDEX) || ((countStation++ % defaultBlockReceiverCount) == 0)))
            newBlock();
        idRecv = iter->idRecv;

        insert(iter++, defaultBlockSizeEpoch);
      }
      blockCountEpoch_.at(idEpoch) = blockIndices_.size() - blockEpochStart + (blockSize(blockCount()-1) ? 1 : 0);
    }

    // receiver interval parameters
    newBlock();
    blockInterval_ = blockCount()-1;
    UInt countStation = 0;
    while((iter != parameters.end()) && (iter->idEpoch == NULLINDEX) && (iter->idRecv != NULLINDEX) && (iter->idTrans == NULLINDEX))
    {
      const UInt idRecv = iter->idRecv;
      if(defaultBlockReceiverCount && ((countStation++ % defaultBlockReceiverCount) == 0))
        newBlock();

      Bool firstBlock = TRUE;
      while((iter != parameters.end()) && (iter->idEpoch == NULLINDEX) && (iter->idRecv == idRecv) && (iter->idTrans == NULLINDEX))
      {
        insert(iter++, firstBlock ? defaultBlockSizeInterval : 0); // do not split parameters of a receiver
        firstBlock = FALSE;
      }
    }

    // transmitter interval parameters
    newBlock();
    while((iter != parameters.end()) && (iter->idEpoch == NULLINDEX) && (iter->idRecv == NULLINDEX) && (iter->idTrans != NULLINDEX))
    {
      const UInt idTrans = iter->idTrans;
      Bool firstBlock = TRUE;
      while((iter != parameters.end()) && (iter->idEpoch == NULLINDEX) && (iter->idRecv == NULLINDEX) && (iter->idTrans == idTrans))
      {
        insert(iter++, firstBlock ? defaultBlockSizeInterval : 0); // do not split parameters of a transmitter
        firstBlock = FALSE;
      }
    }

    // other interval parameters
    newBlock();
    while((iter != parameters.end()) && (iter->idEpoch == NULLINDEX) && (iter->idRecv == NULLINDEX) && (iter->idTrans == NULLINDEX))
      insert(iter++, defaultBlockSizeInterval);

    // ambiguity parameters
    newBlock();
    blockAmbiguity_ = blockCount()-1;
    while((iter != parameters.end()) && (iter->idRecv != NULLINDEX) && (iter->idTrans != NULLINDEX))
      insert(iter++, defaultBlockSizeAmbiguity);

    // remove possible last empty block
    if(blockCount() && !blockSize(blockCount()-1))
      blockIndices_.pop_back();

    parameters.clear();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt GnssNormalEquationInfo::normalsBlockRank(UInt i, UInt k, UInt commSize)
{
  try
  {
    // if(blockRank_.at(i) != NULLINDEX) return blockRank_.at(i);
    // if(blockRank_.at(k) != NULLINDEX) return blockRank_.at(k);

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
