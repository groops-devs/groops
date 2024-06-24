/***********************************************/
/**
* @file slr.h
*
* @brief global navigation satellite system classes.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLR__
#define __GROOPS_SLR__

#include "parallel/parallel.h"
#include "base/parameterName.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "slr/slrObservation.h"
#include "slr/slrDesignMatrix.h"
#include "slr/slrSatellite.h"
#include "slr/slrStation.h"
#include "slr/slrNormalEquationInfo.h"

/** @addtogroup slrGroup */
/// @{

/***** TYPES ***********************************/

class Slr;
typedef std::shared_ptr<Slr> SlrPtr;

class MatrixDistributed;
class SlrSatelliteGenerator;
class SlrStationGenerator;
class SlrParametrization;
class EarthRotation;
class PlatformSelector;
typedef std::shared_ptr<SlrSatelliteGenerator> SlrSatelliteGeneratorPtr;
typedef std::shared_ptr<SlrStationGenerator>   SlrStationGeneratorPtr;
typedef std::shared_ptr<SlrParametrization>    SlrParametrizationPtr;
typedef std::shared_ptr<EarthRotation>         EarthRotationPtr;
typedef std::shared_ptr<PlatformSelector>      PlatformSelectorPtr;

/***** CLASS ***********************************/

class Slr
{
public:
  std::vector<SlrSatellitePtr>                     satellites;
  std::vector<SlrStationPtr>                       stations;
  SlrParametrizationPtr                            parametrization;
  std::function<void(SlrObservationEquation &eqn)> funcReduceModels;
  std::function<Rotary3d(const Time &time)>        funcRotationCrf2Trf;

  Polynomial        polynomialEop;
  std::vector<Time> times;
  Matrix            eop;             // Matrix eop columns: xp, yp, sp, deltaUT, LOD, X, Y, S

  void init(const std::vector<Time> &times, SlrSatelliteGeneratorPtr satelliteGenerator, SlrStationGeneratorPtr stationGenerator,
            EarthRotationPtr earthRotation, SlrParametrizationPtr parametrization);
  Rotary3d rotationCrf2Trf(const Time &time) const; // Inertial system (CRF) -> earth fixed system (TRF).

  void   initParameter            (SlrNormalEquationInfo &normalEquationInfo);
  Vector aprioriParameter         (const SlrNormalEquationInfo &normalEquationInfo) const;
  Bool   basicObservationEquations(const SlrNormalEquationInfo &normalEquationInfo, UInt idStat, UInt idSat, UInt idPass, SlrObservationEquation &eqn) const;
  void   designMatrix             (const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const;
  void   constraints              (const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;
  Double updateParameter          (const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz);
  void   updateCovariance         (const SlrNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance);
  void   writeResults             (const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix="");

  std::vector<Byte> selectStations(PlatformSelectorPtr selector);
  std::vector<Byte> selectSatellites(PlatformSelectorPtr selector);

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
    void print(Double convertToMeter, Double &maxChangeTotal);
  };
};

/// @}

/***********************************************/

#endif /* __GROOPS___ */
