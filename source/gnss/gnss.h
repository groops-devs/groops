/***********************************************/
/**
* @file gnss.h
*
* @brief global navigation satellite system classes.
*
* @author Torsten Mayer-Guerr
* @date 2010-08-03
*
*/
/***********************************************/

#ifndef __GROOPS_GNSS__
#define __GROOPS_GNSS__

#include "parallel/parallel.h"
#include "base/gnssType.h"
#include "base/parameterName.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssNormalEquationInfo.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class Gnss;
typedef std::shared_ptr<Gnss> GnssPtr;

class MatrixDistributed;
class GnssTransmitterGenerator;
class GnssReceiverGenerator;
class GnssParametrization;
class EarthRotation;
class PlatformSelector;
typedef std::shared_ptr<GnssTransmitterGenerator> GnssTransmitterGeneratorPtr;
typedef std::shared_ptr<GnssReceiverGenerator>    GnssReceiverGeneratorPtr;
typedef std::shared_ptr<GnssParametrization>      GnssParametrizationPtr;
typedef std::shared_ptr<EarthRotation>            EarthRotationPtr;
typedef std::shared_ptr<PlatformSelector>         PlatformSelectorPtr;

/***** CLASS ***********************************/

class Gnss
{
public:
  std::vector<Time>               times;           // Epochs
  std::vector<GnssTransmitterPtr> transmitters;    // GNSS satellites
  std::vector<GnssReceiverPtr>    receivers;       // stations & LEOs
  GnssParametrizationPtr          parametrization; // parameters
  std::function<void(GnssObservationEquation &eqn)> funcReduceModels;
  std::function<Rotary3d(const Time &time)>         funcRotationCrf2Trf;
  Matrix                          eop;             // Matrix eop columns: xp, yp, sp, deltaUT, LOD, X, Y, S
  std::vector<std::vector<std::vector<GnssType>>> typesRecvTrans; // for each receiver and transmitter: used types (receiver types)

  void init(const std::vector<Time> &times, const Time &timeMargin,
            GnssTransmitterGeneratorPtr transmitterGenerator, GnssReceiverGeneratorPtr receiverGenerator,
            EarthRotationPtr earthRotation, GnssParametrizationPtr parametrization, Parallel::CommunicatorPtr comm);

  Rotary3d rotationCrf2Trf(const Time &time) const; // Inertial system (CRF) -> earth fixed system (TRF).
  void     synchronizeTransceivers(Parallel::CommunicatorPtr comm);

  void   initParameter            (GnssNormalEquationInfo &normalEquationInfo);
  Vector aprioriParameter         (const GnssNormalEquationInfo &normalEquationInfo) const;
  Bool   basicObservationEquations(const GnssNormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch, GnssObservationEquation &eqn) const;
  void   designMatrix             (const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const;
  void   constraintsEpoch         (const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;
  void   constraints              (const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;
  Double ambiguityResolve         (const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                                   const std::vector<Byte> &selectedTransmitters, const std::vector<Byte> &selectedReceivers,
                                   const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger);
  Double updateParameter          (const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz);
  void   updateCovariance         (const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance);
  void   writeResults             (const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix="");

  /** @brief sorted list of used types. */
  std::vector<GnssType> types(const GnssType mask=GnssType::ALL) const;

  std::vector<Byte> selectTransmitters(PlatformSelectorPtr selector);
  std::vector<Byte> selectReceivers(PlatformSelectorPtr selector);

  class InfoParameterChange
  {
  public:
    std::string unit;
    UInt        count;
    Double      rms;
    Double      maxChange;
    std::string info;

    InfoParameterChange(const std::string &unit) : unit(unit), count(0), rms(0), maxChange(0) {}
    Bool update(Double change);
    void synchronizeAndPrint(Parallel::CommunicatorPtr comm, Double convertToMeter, Double &maxChangeTotal);
  };
};

/// @}

/***********************************************/

#endif /* __GROOPS___ */
