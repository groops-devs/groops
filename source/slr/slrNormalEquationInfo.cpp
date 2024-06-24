/***********************************************/
/**
* @file slrNormalEquationInfo.cpp
*
* @brief SLR normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "parallel/matrixDistributed.h"
#include "base/parameterName.h"
#include "slrNormalEquationInfo.h"

/***********************************************/

SlrNormalEquationInfo::SlrNormalEquationInfo(UInt countStations, UInt countSatellites) :
    estimateStation(countStations, TRUE), estimateSatellite(countSatellites, TRUE), defaultBlockSize(2048)
{
}

/***********************************************/

void SlrNormalEquationInfo::initNewParameterNames()
{
  try
  {
    indexGravity.clear();
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

SlrParameterIndex SlrNormalEquationInfo::addParameters(UInt idStat, UInt idSat, const std::vector<ParameterName> &parameterNames)
{
  if(!parameterNames.size())
    return SlrParameterIndex(NULLINDEX);
  parameters.push_back(std::make_tuple(idStat, idSat, parameters.size(), parameterNames));
  return SlrParameterIndex(parameters.size()-1);
}

/***********************************************/

void SlrNormalEquationInfo::calculateIndex()
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
      const UInt idx   = std::get<2>(*iter);
      index_.at(idx)   = blockIndices_.back();
      block_.at(idx)   = blockCount()-1;
      count_.at(idx)   = std::get<3>(*iter).size();
      parameterNames_.insert(parameterNames_.end(), std::get<3>(*iter).begin(), std::get<3>(*iter).end());
      blockIndices_.back() += std::get<3>(*iter).size();
    };

    block_.resize(parameters.size(), NULLINDEX);
    index_.resize(parameters.size(), NULLINDEX);
    count_.resize(parameters.size(), 0);
    blockIndices_ = {0, 0};
    parameterNames_.reserve(std::accumulate(parameters.begin(), parameters.end(), UInt(0), [](UInt count, const auto &p){return count+std::get<3>(p).size();}));

    std::stable_sort(parameters.begin(), parameters.end(), [](auto &p1, auto &p2)
                    {
                      if(std::get<0>(p1) != std::get<0>(p2)) return (std::get<0>(p1) < std::get<0>(p2)); // idStat
                      return (std::get<1>(p1) < std::get<1>(p2));                                        // idSat
                    });
    auto iter = parameters.begin();

    // station interval parameters
    newBlock();
    for(UInt idStat=0; idStat<estimateStation.size(); idStat++)
    {
      Bool firstBlock = TRUE;
      while((iter != parameters.end()) && (std::get<0>(*iter) == idStat) && (std::get<1>(*iter) == NULLINDEX))
      {
        insert(iter, firstBlock ? defaultBlockSize : 0); // do not split parameters of a station
        firstBlock = FALSE;
        iter++;
      }
    }

    // satellite interval parameters
    newBlock();
    for(UInt idSat=0; idSat<estimateSatellite.size(); idSat++)
    {
      Bool firstBlock = TRUE;
      while((iter != parameters.end()) && (std::get<0>(*iter) == NULLINDEX) && (std::get<1>(*iter) == idSat))
      {
        insert(iter, firstBlock ? defaultBlockSize : 0); // do not split parameters of a satellite
        firstBlock = FALSE;
        iter++;
      }
    }

    // other interval parameters
    newBlock();
    while((iter != parameters.end()) && (std::get<0>(*iter) == NULLINDEX) && (std::get<1>(*iter) == NULLINDEX))
    {
      insert(iter, defaultBlockSize);
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
