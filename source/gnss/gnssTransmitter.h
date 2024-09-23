/***********************************************/
/**
* @file gnssTransmitter.h
*
* @brief GNSS transmitter.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2013-06-28
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSMITTER__
#define __GROOPS_GNSSTRANSMITTER__

#include "base/polynomial.h"
#include "base/gnssType.h"
#include "gnss/gnssTransceiver.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssTransmitter;
typedef std::shared_ptr<GnssTransmitter> GnssTransmitterPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS transmitter.
* eg. GPS satellites. */
class GnssTransmitter : public GnssTransceiver
{
  GnssType                 type; // system + PRN
  Polynomial               polynomial;
  std::vector<Double>      clk;
  std::vector<Vector3d>    offset;   // between CoM and ARF in SRF
  std::vector<Transform3d> crf2srf, srf2arf;

public:
  std::vector<Time> timesPosVel;
  Matrix            pos, vel; // CoM in CRF (epoch times (x,y,z))

  GnssTransmitter(GnssType prn, const Platform &platform,
                  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction,
                  const Vector &useableEpochs, const std::vector<Double> &clock, const std::vector<Vector3d> &offset,
                  const std::vector<Transform3d> &crf2srf, const std::vector<Transform3d> &srf2arf,
                  const std::vector<Time> &timesPosVel, const_MatrixSliceRef position, const_MatrixSliceRef velocity, UInt interpolationDegree)
  : GnssTransceiver(platform, noPatternFoundAction, useableEpochs),
    type(prn), polynomial(timesPosVel, interpolationDegree, TRUE/*throwException*/, FALSE/*leastSquares*/, -(interpolationDegree+1.1), -1.1, 1e-7),
    clk(clock), offset(offset), crf2srf(crf2srf), srf2arf(srf2arf), timesPosVel(timesPosVel), pos(position), vel(velocity) {}

  /// Destructor.
  virtual ~GnssTransmitter() {}

  /** @brief Identify number in the GNSS system. */
  UInt idTrans() const {return id_;}

  /** @brief PRN number of satellite.
  *  = prn + GnssType::SYSTEM. */
  GnssType PRN() const {return type;}

  /** @brief Clock error.
  * error = clock time - system time [s] */
  Double clockError(UInt idEpoch) const {return clk.at(idEpoch);}

  /** @brief set clock error.
  * error = observed clock time - system time [s] */
  void updateClockError(UInt idEpoch, Double deltaClock) {clk.at(idEpoch) += deltaClock;}

  /** @brief center of mass in celestial reference frame (CRF). */
  Vector3d positionCoM(const Time &time) const;

  /** @brief antenna reference point in celestial reference frame (CRF). */
  Vector3d position(UInt idEpoch, const Time &time) const {return positionCoM(time) + crf2srf.at(idEpoch).inverseTransform(offset.at(idEpoch));}

  /** @brief velocity in CRF [m/s]. */
  Vector3d velocity(const Time &time) const;

  /** @brief Rotation from celestial reference frame (CRF) to left-handed antenna system. */
  Transform3d celestial2antennaFrame(UInt idEpoch, const Time &/*time*/) const {return srf2arf.at(idEpoch) * crf2srf.at(idEpoch);}
};

/***********************************************/

inline Vector3d GnssTransmitter::positionCoM(const Time &time) const
{
  try
  {
    return Vector3d(polynomial.interpolate({time}, pos));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d GnssTransmitter::velocity(const Time &time) const
{
  try
  {
    return Vector3d(polynomial.interpolate({time}, vel));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
