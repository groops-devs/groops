/***********************************************/
/**
* @file gnssParametrizationIonosphere.h
*
* @brief GNSS ionosphere representation.
*
* @author Torsten Mayer-Guerr
* @date 2018-11-20
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONIONOSPHERE__
#define __GROOPS_GNSSPARAMETRIZATIONIONOSPHERE__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationIonosphere
static const char *docstringGnssParametrizationIonosphere = R"(
\section{GnssParametrizationIonosphere}\label{gnssParametrizationIonosphereType}

Definition of settings and constants for ionospheric corrections in GNSS observation equations.

See also \program{GnssProcessing} and \program{GnssSimulateReceiver}.
)";
#endif

/***********************************************/

#include "classes/magnetosphere/magnetosphere.h"
#include "gnss/gnss.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssParametrizationIonosphere;
typedef std::shared_ptr<GnssParametrizationIonosphere> GnssParametrizationIonospherePtr;

/***** CLASS ***********************************/

/** @brief GNSS Earth rotation.
* An Instance of this class can be created by @ref readConfig. */
class GnssParametrizationIonosphere : public Gnss::Parametrization
{
  // Corrections
  MagnetospherePtr magnetosphere;
  Bool             apply1stOrder, apply2ndOrder, apply3rdOrder, applyBending;
  Double           singleLayerHeight;
  Double           mapR, mapH, mapAlpha;

  // Constraint
  Double sigmaSTEC;

  // Parametrization
  // ---------------
  FileName fileNameVTEC;
  std::vector<std::vector<Gnss::ParameterIndex>> indexParameterVTEC; // for each receiver, for each epoch
  Matrix VTEC;

  Vector3d intersection(const Vector3d &posRecv, const Vector3d &posTrans) const;
  Double   mapping(const Vector3d &posRecv, Angle elevation) const;

public:
  /// Constructor.
  GnssParametrizationIonosphere(Config &config, const std::string &name);

  /// Destructor.
  virtual ~GnssParametrizationIonosphere() {}

  // Realization of Gnss::Parametrization
  // ------------------------------------
  void   initIntervalEarly(Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  void   initParameter(Gnss::NormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  Bool   isDesignMatrix(const Gnss::NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch) const override;
  void   designMatrix(const Gnss::NormalEquationInfo &normalEquationInfo, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const override;
  void   eliminateTecParameter(const Gnss::NormalEquationInfo &normalEquationInfo, Gnss::ObservationEquation &eqn) const;
  void   updateAndEliminateTecParameter(const Gnss::NormalEquationInfo &normalEquationInfo, Gnss::ObservationEquation &eqn,
                                        Vector &We, Matrix &AWz, Double &maxChangeTec, std::string &infoMaxChangeTec) const;
  Double updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics) override;
  void   writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix) override;

  Vector slantDelay(const Gnss::ObservationEquation &eqn, Matrix &B) const;

  /** @brief creates an derived instance of this class. */
  static GnssParametrizationIonospherePtr create(Config &config, const std::string &name) {return std::make_shared<GnssParametrizationIonosphere>(config, name);}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssParametrizationIonosphere.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssParametrizationIonosphere */
template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationIonospherePtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
