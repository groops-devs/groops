/***********************************************/
/**
* @file gnssParametrizationGravityField.h
*
* @brief GNSS gravity field.
*
* @author Torsten Mayer-Guerr
* @date 2014-05-26
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONGRAVITYFIELD__
#define __GROOPS_GNSSPARAMETRIZATIONGRAVITYFIELD__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationGravityField
static const char *docstringGnssParametrizationGravityField = R"(
\section{GnssParametrizationGravityField}\label{gnssParametrizationGravityFieldType}

Estimation of low-degree gravity field parameters via GNSS.

See also \program{GnssProcessing}.
)";
#endif

/***********************************************/

#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "gnss/gnss.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssParametrizationGravityField;
typedef std::shared_ptr<GnssParametrizationGravityField> GnssParametrizationGravityFieldPtr;

/***** CLASS ***********************************/

/** @brief gravity field.
* An Instance of this class can be created by @ref readConfig. */
class GnssParametrizationGravityField : public Gnss::Parametrization
{
  ParametrizationGravityPtr parametrizationPtr;
  Vector                    x;
  Gnss::ParameterIndex      index;

public:
  /// Constructor.
  GnssParametrizationGravityField(Config &config, const std::string &name);

  /** @brief Gravity field parametrization. */
  ParametrizationGravityPtr parametrization() const {return parametrizationPtr;}

  /** @brief Index of parameters in the normal equation system.
  * Interval parameters. */
  Gnss::ParameterIndex normalEquationIndex() const {return index;}

  /** @brief Estimated parameters. */
  Vector parameter() const {return x;}

  // Realization of Gnss::Parametrization
  // ------------------------------------
  void   initParameter(Gnss::NormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  Double updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics) override;
  void   writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix) override;

  /** @brief creates an derived instance of this class. */
  static GnssParametrizationGravityFieldPtr create(Config &config, const std::string &name) {return std::make_shared<GnssParametrizationGravityField>(config, name);}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssParametrizationGravityField.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssParametrizationGravityField */
template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationGravityFieldPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
