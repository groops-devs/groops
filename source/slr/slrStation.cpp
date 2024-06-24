/***********************************************/
/**
* @file slrStation.cpp
*
* @brief SLR station.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

// #include <random>
#include "base/import.h"
#include "base/string.h"
#include "parser/expressionParser.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "inputOutput/logging.h"
#include "inputOutput/system.h"
#include "misc/varianceComponentEstimation.h"
#include "slr/slrObservation.h"
#include "slr/slrSatellite.h"
#include "slr/slrStation.h"

/***********************************************/

SlrStation::SlrStation(const Platform &platform, const std::vector<Time> &times, const Vector3d &pos, const_MatrixSliceRef offset,
                       ExpressionVariablePtr accuracyExpr, UInt interpolationDegree)
  : SlrPlatform(platform), polynomial(times, interpolationDegree, FALSE/*throwException*/), accuracyExpr(accuracyExpr), times(times), pos(pos), offset(offset)
{
  global2local = inverse(localNorthEastUp(platform.approxPosition, Ellipsoid()));
}

/***********************************************/

void SlrStation::disable(const std::string &reason)
{
  try
  {
    SlrPlatform::disable(reason);
    if(!reason.empty())
      disableReason = reason;
    observations.clear();
    observations.shrink_to_fit();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d SlrStation::position(const Time &time) const
{
  try
  {
    return pos + global2local.inverseTransform(Vector3d(polynomial.interpolate({time}, offset)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Transform3d SlrStation::global2localFrame(const Time &/*time*/) const
{
  return global2local;
}

/***********************************************/

std::vector<Time> SlrStation::correctedTimes(const std::vector<Time> &times) const
{
  try
  {
    if(!timeBiases.size())
      return times;

    std::vector<Time> timesCorrected = times;
    Vector bias = polynomial.interpolate({times}, timeBiases);
    for(UInt i=0; i<times.size(); i++)
      timesCorrected.at(i) -= seconds2time(bias(i));
    return timesCorrected;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SlrStation::accuracy(const Time &/*time*/, Double residual, Double accuracy, Double redundancy,
                            Double laserWavelength, Angle azimut, Angle elevation) const
{
  try
  {
    accuracyVariableList.setVariable("residual",        residual);
    accuracyVariableList.setVariable("accuracy",        accuracy);
    accuracyVariableList.setVariable("redundancy",      redundancy);
    accuracyVariableList.setVariable("laserWavelength", laserWavelength);
    accuracyVariableList.setVariable("azimut",          azimut);
    accuracyVariableList.setVariable("elevation",       elevation);
    return accuracyExpr->evaluate(accuracyVariableList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
