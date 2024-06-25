/***********************************************/
/**
* @file slrNormalEquationInfo.h
*
* @brief SLR normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRNORMALEQUATIONINFO__
#define __GROOPS_SLRNORMALEQUATIONINFO__

#include <regex>
#include "parallel/parallel.h"
#include "parallel/matrixDistributed.h"
#include "base/parameterName.h"

/** @addtogroup slrGroup */
/// @{

/***** CLASS ***********************************/

class SlrParameterIndex
{
  UInt index;

public:
  SlrParameterIndex(UInt idx=NULLINDEX) : index(idx) {}
  explicit operator bool() const {return index != NULLINDEX;}

  friend class SlrNormalEquationInfo;
  friend class SlrDesignMatrix;
};

/***** CLASS ***********************************/

class SlrNormalEquationInfo
{
public:
  std::vector<std::pair<std::regex, Bool>> enableParametrizations; // wildcard matching, enable/disable
  std::vector<Byte>                        estimateStation;
  std::vector<Byte>                        estimateSatellite;
  UInt                                     defaultBlockSize;
  std::vector<SlrParameterIndex>           indexGravity;

  SlrNormalEquationInfo(UInt countStations, UInt countSatellites);

  void initNewParameterNames();
  SlrParameterIndex parameterNamesStation  (UInt idStat, const std::vector<ParameterName> &parameterNames) {return addParameters(idStat,    NULLINDEX, parameterNames);}
  SlrParameterIndex parameterNamesSatellite(UInt idSat,  const std::vector<ParameterName> &parameterNames) {return addParameters(NULLINDEX, idSat,     parameterNames);}
  SlrParameterIndex parameterNamesOther    (const std::vector<ParameterName> &parameterNames)              {return addParameters(NULLINDEX, NULLINDEX, parameterNames);}
  void calculateIndex();

  UInt block(const SlrParameterIndex &index) const {return block_.at(index.index);}
  UInt index(const SlrParameterIndex &index) const {return index_.at(index.index);}
  UInt count(const SlrParameterIndex &index) const {return count_.at(index.index);}

  const std::vector<ParameterName> &parameterNames() const {return parameterNames_;}
  const std::vector<UInt>          &blockIndices()   const {return blockIndices_;}

  UInt parameterCount()   const {return blockIndices_.back();}                      //!< Number of rows/columns (dimension) of distributed matrix
  UInt blockIndex(UInt i) const {return blockIndices_.at(i);}                       //!< Start index of block @a i
  UInt blockSize(UInt i)  const {return blockIndices_.at(i+1)-blockIndices_.at(i);} //!< Size of block @a i
  UInt blockCount()       const {return blockIndices_.size()-1;}                    //!< Number of block rows/columns

private:
  std::vector<std::tuple<UInt, UInt, UInt, std::vector<ParameterName>>> parameters; // idStat, idSat, idx, name
  std::vector<UInt>          block_, index_, count_;
  std::vector<ParameterName> parameterNames_;
  std::vector<UInt>          blockIndices_;
  std::vector<UInt>          blockCountEpoch_;

  SlrParameterIndex addParameters(UInt idStat, UInt idSat, const std::vector<ParameterName> &parameterNames);
};

/// @}

/***********************************************/

#endif
