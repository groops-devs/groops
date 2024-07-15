/***********************************************/
/**
* @file slrStation.h
*
* @brief SLR station.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRSTATION__
#define __GROOPS_SLRSTATION__

#include "parser/expressionParser.h"
#include "files/fileInstrument.h"
#include "slr/slrObservation.h"
#include "slr/slrPlatform.h"
#include "slr/slrSatellite.h"

/** @addtogroup slrGroup */
/// @{

/***** TYPES ***********************************/

class SlrStation;
typedef std::shared_ptr<SlrStation> SlrStationPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for SLR station. */
class SlrStation : public SlrPlatform
{
  Polynomial            polynomial;
  Transform3d           global2local;
  mutable VariableList  accuracyVariableList;
  ExpressionVariablePtr accuracyExpr;

public:
  std::vector<Time>         times;      // of offset and timeBiases
  Vector3d                  pos;        // regularized marker pos in global system
  Matrix                    offset;     // pos to SRF in local system
  Vector                    timeBiases; // Observed time - corrected time [seconds]
  std::vector<std::string>  preprocessingInfos;
  std::string               disableReason;
  std::vector<std::vector<SlrObservationPtr>> observations; // observations at station (for each satellite, for each pass)

  SlrStation(const Platform &platform, const std::vector<Time> &times, const Vector3d &pos, const_MatrixSliceRef offset,
             ExpressionVariablePtr accuracyExpr, UInt interpolationDegree);

  /// Destructor.
  virtual ~SlrStation() {}

  /** @brief Identify number in the SLR system. */
  UInt idStat() const {return id_;}

  /** @brief Disable station completely. */
  void disable(const std::string &reason) override;

  /** @brief system reference point (SRF) in TRF. */
  Vector3d position(const Time &time) const;

  /** @brief Rotation from terrestrial reference frame (TRF) to local frame (north, east, up). */
  Transform3d global2localFrame(const Time &time) const;

  /** @brief Corrected for system time biases. */
  std::vector<Time> correctedTimes(const std::vector<Time> &times) const;

  /** @brief Direction (and other parameters) dependent standard deviation.
  * @a azmiut and @a elevation must be given in the antenna frame (left-handed). */
  Double accuracy(const Time &time, Double residual, Double accuracy, Double redundancy,
                  Double laserWavelength, Angle azimut, Angle elevation) const;
}; //class SlrStation

/***********************************************/

/// @}

#endif /* __GROOPS___ */
