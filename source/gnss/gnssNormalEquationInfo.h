/***********************************************/
/**
* @file gnssNormalEquationInfo.h
*
* @brief GNSS normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2010-08-03
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSNORMALEQUATIONINFO__
#define __GROOPS_GNSSNORMALEQUATIONINFO__

#include <regex>
#include "parallel/parallel.h"
#include "parallel/matrixDistributed.h"
#include "base/parameterName.h"

/** @addtogroup gnssGroup */
/// @{

/***** CLASS ***********************************/

class GnssParameterIndex
{
  UInt index;
public:
  GnssParameterIndex(UInt idx=NULLINDEX) : index(idx) {}
  explicit operator bool() const {return index != NULLINDEX;}

  friend class GnssNormalEquationInfo;
  friend class GnssDesignMatrix;
};

/***** CLASS ***********************************/

class GnssNormalEquationInfo
{
public:
  Parallel::CommunicatorPtr comm;
  std::vector<std::pair<std::regex, Bool>> enableParametrizations; // wildcard matching, enable/disable
  Bool                      isEachReceiverSeparately; // only receiver parameters are allowed
  std::vector<Byte>         estimateReceiver;         // subset of observations
  std::vector<UInt>         idEpochs;                 // epochs to estimate
  UInt                      defaultBlockSizeEpoch;
  UInt                      defaultBlockSizeInterval;
  UInt                      defaultBlockSizeAmbiguity;
  UInt                      defaultBlockReceiverCount;
  UInt                      defaultBlockCountReduction;
  Bool                      keepEpochNormalsInMemory;
  Bool                      accumulateEpochObservations;

  GnssNormalEquationInfo(UInt countEpoch, UInt countReceiver, UInt countTransmitter, Parallel::CommunicatorPtr comm);

  void initNewParameterNames();
  GnssParameterIndex parameterNamesEpochReceiver   (UInt idEpoch, UInt idRecv,  const std::vector<ParameterName> &parameterNames)              {return addParameters(idEpoch,   idRecv,    NULLINDEX, parameterNames);}
  GnssParameterIndex parameterNamesEpochTransmitter(UInt idEpoch, UInt idTrans, const std::vector<ParameterName> &parameterNames)              {return addParameters(idEpoch,   NULLINDEX, idTrans,   parameterNames);}
  GnssParameterIndex parameterNamesEpochOther      (UInt idEpoch, const std::vector<ParameterName> &parameterNames)                            {return addParameters(idEpoch,   NULLINDEX, NULLINDEX, parameterNames);}
  GnssParameterIndex parameterNamesReceiver        (UInt idRecv,  const std::vector<ParameterName> &parameterNames)                            {return addParameters(NULLINDEX, idRecv,    NULLINDEX, parameterNames);}
  GnssParameterIndex parameterNamesTransmitter     (UInt idTrans, const std::vector<ParameterName> &parameterNames)                            {return addParameters(NULLINDEX, NULLINDEX, idTrans,   parameterNames);}
  GnssParameterIndex parameterNamesOther           (const std::vector<ParameterName> &parameterNames)                                          {return addParameters(NULLINDEX, NULLINDEX, NULLINDEX, parameterNames);}
  GnssParameterIndex parameterNamesAmbiguity       (UInt idEpoch, UInt idRecv, UInt idTrans, const std::vector<ParameterName> &parameterNames) {return addParameters(idEpoch,   idRecv,    idTrans,   parameterNames);}
  void calculateIndex(const Vector &recvProcess);

  UInt block(const GnssParameterIndex &index) const {return block_.at(index.index);}
  UInt index(const GnssParameterIndex &index) const {return index_.at(index.index);}
  UInt count(const GnssParameterIndex &index) const {return count_.at(index.index);}

  const std::vector<ParameterName> &parameterNames() const {return parameterNames_;}
  const std::vector<UInt>          &blockIndices()   const {return blockIndices_;}

  UInt parameterCount()   const {return blockIndices_.back();}                      //!< Number of rows/columns (dimension) of distributed matrix
  UInt blockIndex(UInt i) const {return blockIndices_.at(i);}                       //!< Start index of block @a i
  UInt blockSize(UInt i)  const {return blockIndices_.at(i+1)-blockIndices_.at(i);} //!< Size of block @a i
  UInt blockCount()       const {return blockIndices_.size()-1;}                    //!< Number of block rows/columns

  UInt blockCountEpoch(UInt idEpoch) const {return blockCountEpoch_.at(idEpoch);}
  UInt blockInterval()  const {return blockInterval_;}
  UInt blockAmbiguity() const {return blockAmbiguity_;}

  UInt normalsBlockRank(UInt i, UInt k, UInt commSize);

private:
  class Parameter
  {
  public:
    UInt idEpoch, idRecv, idTrans, rank, idx;
    std::vector<ParameterName> names;
  };

  std::list<Parameter>       parameters;
  std::vector<UInt>          block_, index_, count_;  // for each parameter index
  std::vector<ParameterName> parameterNames_;
  std::vector<UInt>          blockIndices_;           // for each block
  std::vector<UInt>          blockRank_;              // for each block
  std::vector<UInt>          blockCountEpoch_;        // for each epoch
  UInt                       blockInterval_, blockAmbiguity_;

  GnssParameterIndex addParameters(UInt idEpoch, UInt idRecv, UInt idTrans, const std::vector<ParameterName> &parameterNames);
};

/// @}

/***********************************************/

#endif
