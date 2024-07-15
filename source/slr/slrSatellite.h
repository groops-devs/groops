/***********************************************/
/**
* @file slrSatellite.h
*
* @brief SLR satellite.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRSATELLITE__
#define __GROOPS_SLRSATELLITE__

#include "base/polynomial.h"
#include "slr/slrPlatform.h"

/** @addtogroup slrGroup */
/// @{

/***** TYPES ***********************************/

class SlrSatellite;
typedef std::shared_ptr<SlrSatellite> SlrSatellitePtr;

/***** CLASS ***********************************/

/** @brief Abstract class for SLR satellite.
* eg. GPS satellites. */
class SlrSatellite : public SlrPlatform
{
public:
  Polynomial        polynomial;
  std::vector<Time> times;
  Matrix            pos, vel; // CoM in CRF (epoch x (x,y,z))
  Matrix            srf2crf;

  SlrSatellite(const Platform &platform, const std::vector<Time> &times,
               const_MatrixSliceRef position, const_MatrixSliceRef velocity, const_MatrixSliceRef srf2crf, UInt interpolationDegree)
    : SlrPlatform(platform), polynomial(times, interpolationDegree, FALSE/*throwException*/), times(times), pos(position), vel(velocity), srf2crf(srf2crf) {}

  /// Destructor.
  virtual ~SlrSatellite() {}

  /** @brief Identify number in the SLR system. */
  UInt idSat() const {return id_;}

  /** @brief center of mass in celestial reference frame (CRF). */
  Vector3d position(const Time &time) const;

  /** @brief velocity in CRF [m/s]. */
  Vector3d velocity(const Time &time) const;

  /** @brief Vector from center of mass to the reflector in CRF [m].
  * With corrections (e.g. range bias).
  * @param time of reflection
  * @param direction towards the reflector in CRF. */
  Vector3d centerOfMass2Reflector(const Time &time, const Vector3d &direction) const;
};

/***********************************************/

inline Vector3d SlrSatellite::position(const Time &time) const
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

inline Vector3d SlrSatellite::velocity(const Time &time) const
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

inline Vector3d SlrSatellite::centerOfMass2Reflector(const Time &time, const Vector3d &direction) const
{
  try
  {
    auto lrr = platform.findEquipment<PlatformLaserRetroReflector>(time);
    if(!lrr)
      throw(Exception("no LRR found for "+platform.name));
    const Rotary3d srf2crf = Rotary3d(polynomial.interpolate({time}, this->srf2crf).trans());
    Vector3d offset  = srf2crf.rotate(lrr->position-platform.referencePoint(time));

    if(lrr->range.size() == 1)
      offset -= lrr->range(0,0) * direction;
    else if(lrr->range.size() > 1)
    {
      const Vector3d d = lrr->platform2reflectorFrame.transform(srf2crf.inverseRotate(direction));
      Double tauA      = std::fmod(Double(d.lambda())/(2*PI)+1,1) * lrr->range.rows();
      Double tauZ      = (PI/2-Double(d.phi()))/Double(lrr->dZenit);
      const UInt idxA1 = static_cast<UInt>(std::floor(tauA));
      const UInt idxA2 = (idxA1+1)%lrr->range.rows();
      const UInt idxZ1 = std::min(static_cast<UInt>(std::floor(tauZ)), lrr->range.columns()-1);
      const UInt idxZ2 = std::min(idxZ1+1, lrr->range.columns()-1);
      tauA -= std::floor(tauA);
      tauZ -= std::floor(tauZ);

      // bilinear interpolation
      offset -= ((1-tauA) * (1-tauZ) * lrr->range(idxA1, idxZ1)
               + (1-tauA) *  (tauZ)  * lrr->range(idxA1, idxZ2)
               +  (tauA)  * (1-tauZ) * lrr->range(idxA2, idxZ1)
               +  (tauA)  *  (tauZ)  * lrr->range(idxA2, idxZ2)) * direction;
    }

    return offset;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
