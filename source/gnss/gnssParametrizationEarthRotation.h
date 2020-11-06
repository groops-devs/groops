/***********************************************/
/**
* @file gnssParametrizationEarthRotation.h
*
* @brief GNSS Earth rotation.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2014-05-25
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONEARTHROTATION__
#define __GROOPS_GNSSPARAMETRIZATIONEARTHROTATION__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationEarthRotation
static const char *docstringGnssParametrizationEarthRotation = R"(
\section{GnssParametrizationEarthRotation}\label{gnssParametrizationEarthRotationType}

Earth rotation parameters (ERPs) can be estimated by defining
\config{estimatePole} ($x_p$, $y_p$) and \config{estimateUT1} (dUT1, LOD).

Estimating length of day (LOD) with the sign according to IGS conventions requires a negative
value in \configClass{parametrizationTemporal:trend:timeStep}{parametrizationTemporalType:trend}.

Constraints on the defined parameters can be added via
\configClass{gnssParametrizationConstraints}{gnssParametrizationConstraintsType}.
An example would be to set up \configClass{estimateUT1:constant}{parametrizationTemporalType:constant}
so the \emph{dUT1} parameter is included in the normal equation system . Since \emph{dUT1} cannot be
determined by GNSS, a hard constraint to its a priori value can then be added.

See also \program{GnssProcessing} and \program{GnssSimulateReceiver}.
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "gnss/gnss.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssParametrizationEarthRotation;
typedef std::shared_ptr<GnssParametrizationEarthRotation> GnssParametrizationEarthRotationPtr;

/***** CLASS ***********************************/

/** @brief GNSS Earth rotation.
* An Instance of this class can be created by @ref readConfig. */
class GnssParametrizationEarthRotation : public Gnss::Parametrization
{
  FileName fileNameTimeSeriesEOP;

  // quaternion interpolation
  // ------------------------
  EarthRotationPtr  earthRotationPtr;
  EarthRotationPtr  earthRotationModelsPtr;
  Polynomial        polynomial;
  std::vector<Time> times;

  // Earth Rotation Parameters (ERP)
  // -------------------------------
  ParametrizationTemporalPtr parametrizationPole;
  ParametrizationTemporalPtr parametrizationUT1;
  ParametrizationTemporalPtr parametrizationNutation;

  Gnss::ParameterIndex indexParameterPole;
  Gnss::ParameterIndex indexParameterUT1;
  Gnss::ParameterIndex indexParameterNutation;

  Matrix eopInitial;
  Matrix eopModels;
  Matrix eop;

  void designMatrixTemporal(ParametrizationTemporalPtr parametrization, const Time &time, const_MatrixSliceRef B,
                            const Gnss::ParameterIndex &index, Gnss::DesignMatrix &A) const;

  void rotaryMatrices(const Time &timeGPS, const_MatrixSliceRef eop, Rotary3d &rotW, Rotary3d &rotERA,
                      Rotary3d &rotQ, Rotary3d &rotS, Rotary3d &rotsp) const;

public:
  /// Constructor.
  GnssParametrizationEarthRotation(Config &config, const std::string &name);

  /// Destructor.
  virtual ~GnssParametrizationEarthRotation() {}

  /** @brief Earth orientation parameter.
  * All models and corrections applied. */
  void earthOrientationParameter(const Time &time, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const;

  /** @brief Rotary matrix.
  * Inertial system (CRF) -> earth fixed system (TRF).
  * @param time modified julian date (MJD) in GPS time system. */
  Rotary3d rotaryMatrix(const Time &time) const;

  /** @brief Reference to EarthRotation.
  * Without estimated parameters. */
  EarthRotationPtr earthRotation() const {return earthRotationPtr;}

  // Realization of Gnss::Parametrization
  // ------------------------------------
  void   initIntervalEarly(Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  void   initParameter(Gnss::NormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  Bool   isDesignMatrix(const Gnss::NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch) const override;
  void   designMatrix(const Gnss::NormalEquationInfo &normalEquationInfo, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const override;
  Double updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics) override;
  void   writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix) override;

  /** @brief creates an derived instance of this class. */
  static GnssParametrizationEarthRotationPtr create(Config &config, const std::string &name) {return std::make_shared<GnssParametrizationEarthRotation>(config, name);}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssParametrizationEarthRotation.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssParametrizationEarthRotation */
template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationEarthRotationPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
