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
  parameters.push_back(std::make_tuple(idEpoch, idRecv, idTrans, parameters.size(), parameterNames));
  return GnssParameterIndex(parameters.size()-1);
}

/***********************************************/

void GnssNormalEquationInfo::calculateIndex()
{
  try
  {
    auto newBlock = [&]()
    {
      const UInt parameterCount = this->parameterCount();
      if(blockSize(blockCount()-1))
        blockIndices_.push_back(parameterCount);
    };

    auto insert = [&](auto iter, UInt defaultBlockSize)
    {
      if((defaultBlockSize > 0) && (blockSize(blockCount()-1) >= defaultBlockSize))
        newBlock();
      const UInt idx   = std::get<3>(*iter);
      index_.at(idx)   = blockIndices_.back();
      block_.at(idx)   = blockCount()-1;
      count_.at(idx)   = std::get<4>(*iter).size();
      parameterNames_.insert(parameterNames_.end(), std::get<4>(*iter).begin(), std::get<4>(*iter).end());
      blockIndices_.back() += std::get<4>(*iter).size();
    };

    block_.resize(parameters.size(), NULLINDEX);
    index_.resize(parameters.size(), NULLINDEX);
    count_.resize(parameters.size(), 0);
    blockIndices_ = {0, 0};
    parameterNames_.reserve(std::accumulate(parameters.begin(), parameters.end(), UInt(0), [](UInt count, const auto &p){return count+std::get<4>(p).size();}));

    std::stable_sort(parameters.begin(), parameters.end(), [](auto &p1, auto &p2)
                    {
                      const Bool isAmbi1 = (std::get<1>(p1) != NULLINDEX) && (std::get<2>(p1) != NULLINDEX);
                      const Bool isAmbi2 = (std::get<1>(p2) != NULLINDEX) && (std::get<2>(p2) != NULLINDEX);
                      if(isAmbi1 != isAmbi2) return isAmbi2; // ambiguities always at end
                      if(std::get<0>(p1) != std::get<0>(p2)) return (std::get<0>(p1) < std::get<0>(p2)); // epoch
                      if(std::get<1>(p1) != std::get<1>(p2)) return (std::get<1>(p1) < std::get<1>(p2)); // idRecv
                      return (std::get<2>(p1) < std::get<2>(p2));                                        // idTrans
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
      while((iter != parameters.end()) && (std::get<0>(*iter) == idEpoch) && ((std::get<1>(*iter) == NULLINDEX) || (std::get<2>(*iter) == NULLINDEX)))
      {
        if(std::get<1>(*iter) != idRecv) // next receiver?
          if(defaultBlockReceiverCount && ((std::get<1>(*iter) == NULLINDEX) || ((countStation++ % defaultBlockReceiverCount) == 0)))
            newBlock();
        idRecv = std::get<1>(*iter);

        insert(iter, defaultBlockSizeEpoch);
        iter++;
      }
      blockCountEpoch_.at(idEpoch) = blockIndices_.size() - blockEpochStart + (blockSize(blockCount()-1) ? 1 : 0);
    }

    // receiver interval parameters
    newBlock();
    blockInterval_ = blockCount()-1;
    UInt countStation = 0;
    for(UInt idRecv=0; idRecv<estimateReceiver.size(); idRecv++)
    {
      if(defaultBlockReceiverCount && ((countStation++ % defaultBlockReceiverCount) == 0))
        newBlock();
      Bool firstBlock = TRUE;
      while((iter != parameters.end()) && (std::get<0>(*iter) == NULLINDEX) && (std::get<1>(*iter) == idRecv) && (std::get<2>(*iter) == NULLINDEX))
      {
        insert(iter, firstBlock ? defaultBlockSizeInterval : 0); // do not split parameters of a receiver
        firstBlock = FALSE;
        iter++;
      }
    }

    // transmitter interval parameters
    newBlock();
    for(UInt idTrans=0; idTrans<countTransmitter_; idTrans++)
    {
      Bool firstBlock = TRUE;
      while((iter != parameters.end()) && (std::get<0>(*iter) == NULLINDEX) && (std::get<1>(*iter) == NULLINDEX) && (std::get<2>(*iter) == idTrans))
      {
        insert(iter, firstBlock ? defaultBlockSizeInterval : 0); // do not split parameters of a transmitter
        firstBlock = FALSE;
        iter++;
      }
    }

    // other interval parameters
    newBlock();
    while((iter != parameters.end()) && (std::get<0>(*iter) == NULLINDEX) && (std::get<1>(*iter) == NULLINDEX) && (std::get<2>(*iter) == NULLINDEX))
    {
      insert(iter, defaultBlockSizeInterval);
      iter++;
    }

    // ambiguity parameters
    newBlock();
    blockAmbiguity_ = blockCount()-1;
    while((iter != parameters.end()) && (std::get<1>(*iter) != NULLINDEX) && (std::get<2>(*iter) != NULLINDEX))
    {
      insert(iter, defaultBlockSizeAmbiguity);
      iter++;
    }

    // remove possible last empty block
    if(blockCount() && !blockSize(blockCount()-1))
      blockIndices_.pop_back();

    parameters.clear();
    parameters.shrink_to_fit();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
