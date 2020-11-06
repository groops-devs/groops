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

#include "gnss/gnss.h"

/** @addtogroup gnssGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS transmitter.
* eg. GPS satellites. */
class Gnss::Transmitter
{
  Gnss   *_gnss;    // is set by Gnss::registerTransmitterAndReceiver
  UInt    _idTrans; // is set by Gnss::registerTransmitterAndReceiver

public:
/// Constructor.
Transmitter();

/// Destructor.
virtual ~Transmitter();

/** @brief is called by Gnss::init(). */
void registerGnss(Gnss *gnss, UInt idTrans) {_gnss = gnss; _idTrans = idTrans;}

/** @brief Reference to the complete GNSS system. */
Gnss &gnss()    const;

/** @brief Identify number in the GNSS system. */
UInt idTrans() const {return _idTrans;}

/** @brief PRN number of satellite.
*  = prn + GnssType::SYSTEM. */
virtual GnssType PRN() const = 0; // prn + GnssType::SYSTEM

/** @brief name of satellite. */
virtual std::string name() const = 0;

/** @brief Is the transmitter usable at given epoch (or all epochs). */
virtual Bool useable(UInt idEpoch=NULLINDEX) const = 0;

/** @brief antenna reference point in CRF. */
virtual Vector3d position(UInt idEpoch, const Time &time) const = 0;

/** @brief velocity in CRF [m/s]. */
virtual Vector3d velocity(UInt idEpoch, const Time &time) const = 0;

/** @brief Rotation from celestial reference frame (CRF) to left-handed antenna system. */
virtual Transform3d celestial2antennaFrame(UInt idEpoch, const Time &time) const = 0;

/** @brief Clock error.
* error = transmitter clock time - system time [s] */
virtual Double clockError(UInt idEpoch, const Time &time) const = 0;

/** @brief Transmitted signal types. Empty if no GNSS receiver definition was provided. */
virtual std::vector<GnssType> definedTypes(UInt idEpoch) const = 0;

/** @brief Direction dependent corrections.
* observed range = range (ARPs of transmitter and receiver)  + antennaVariations. */
virtual Vector antennaVariations(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const = 0;

/** @brief Position and Time of transmitter.
* For given receiver position and time the corresponding transmit time and position is computed */
void transmitTime(UInt idEpoch, const Time &timeRecv, const Vector3d &posRecv, Time &timeTransmit, Vector3d &posTransmit) const;

/** @brief Parameters of transmitter of this observation. */
virtual Bool isDesignMatrixTransmitter(const NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idEpoch) const = 0;

/** @brief Fill in the design matrix with parameters of transmitter. */
virtual void designMatrixTransmitter(const NormalEquationInfo &normalEquationInfo, const ObservationEquation &eqn, DesignMatrix &A) const = 0;

/** @brief Receivers are able to track full cycle integer ambiguities. */
virtual Bool supportsIntegerAmbiguities(const NormalEquationInfo &normalEquationInfo) const = 0;

/** @brief Are transmitter code biases estimated?. */
virtual Bool isCodeBiasEstimated(const NormalEquationInfo &normalEquationInfo) const = 0;

/** @brief Are transmitter (float) phase biases estimated?. */
virtual Bool isPhaseBiasEstimated(const NormalEquationInfo &normalEquationInfo) const = 0;
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
