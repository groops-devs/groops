/***********************************************/
/**
* @file gnssTransceiver.h
*
* @brief GNSS receiver or transmitter.
*
* @author Torsten Mayer-Guerr
* @date 2021-01-30
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSCEIVER__
#define __GROOPS_GNSSTRANSCEIVER__

#include "base/gnssType.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileGnssStationInfo.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssTransceiver;
typedef std::shared_ptr<GnssTransceiver> GnssTransceiverPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS receiver or transmitter.
* eg. GPS satellites. */
class GnssTransceiver
{
   std::string name_;
   Vector      useableEpochs;
   UInt        countUseableEpochs;
   GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;

public:
   UInt            id_; // set by Gnss::init()
   GnssStationInfo info;
   GnssSignalBias  signalBias;

public:
  /// Constructor.
  GnssTransceiver(const std::string &name, const GnssStationInfo &info,
                  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs);

  /// Destructor.
  virtual ~GnssTransceiver() {}

  /** @brief name. */
  std::string name() const {return name_;}

  /** @brief Is the platform usable at given epoch (or all epochs). */
  Bool useable(UInt idEpoch=NULLINDEX) const {return countUseableEpochs && ((idEpoch == NULLINDEX) || useableEpochs(idEpoch));}

  /** @brief Disable given epoch (or all epochs). */
  virtual void disable(UInt idEpoch=NULLINDEX);

  /** @brief Allowed signal types. Empty if no  definition was provided. */
  std::vector<GnssType> definedTypes(const Time &time) const;

  /** @brief Direction dependent corrections.
  * observed range = range (ARPs of transmitter and receiver) + antennaVariations. */
  Vector antennaVariations(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const;

  /** @brief Direction (and other parameters) dependent standard deviation.
  * @a azmiut and @a elevation must be given in the antenna frame (left-handed). */
  Vector accuracy(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const;

  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

inline GnssTransceiver::GnssTransceiver(const std::string &name, const GnssStationInfo &info,
                                        GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs)
  : name_(name), useableEpochs(useableEpochs), countUseableEpochs(sum(useableEpochs)),
    noPatternFoundAction(noPatternFoundAction), info(info) {}

/***********************************************/

inline void GnssTransceiver::disable(UInt idEpoch)
{
  try
  {
    if(idEpoch != NULLINDEX)
    {
      if(useableEpochs(idEpoch))
        countUseableEpochs--;
      useableEpochs(idEpoch) = FALSE;
    }
    else
    {
      countUseableEpochs = 0;
      useableEpochs.setNull();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<GnssType> GnssTransceiver::definedTypes(const Time &time) const
{
  try
  {
    const UInt idRecv = info.findReceiver(time);
    if((idRecv == NULLINDEX) || !info.receiver.at(idRecv).receiverDef)
      return std::vector<GnssType>();
    return info.receiver.at(idRecv).receiverDef->types;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector GnssTransceiver::antennaVariations(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &types) const
{
  try
  {
    Vector corr(types.size());
    corr += info.antennaVariations(time, azimut, elevation, types, noPatternFoundAction);
    corr += signalBias.compute(types);
    return corr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector GnssTransceiver::accuracy(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &types) const
{
  try
  {
    return info.accuracy(time, azimut, elevation, types, noPatternFoundAction);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssTransceiver::save(OutArchive &oa) const
{
  oa<<nameValue("useableEpochs",      useableEpochs);
  oa<<nameValue("countUseableEpochs", countUseableEpochs);
  oa<<nameValue("signalBias",         signalBias);
}

/***********************************************/

inline void GnssTransceiver::load(InArchive  &ia)
{
  ia>>nameValue("useableEpochs",      useableEpochs);
  ia>>nameValue("countUseableEpochs", countUseableEpochs);
  ia>>nameValue("signalBias",         signalBias);
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
