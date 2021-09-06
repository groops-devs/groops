/***********************************************/
/**
* @file gnssParametrizationIonosphereMap.h
*
* @brief IonosphereMap.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONIONOSPHEREMAP__
#define __GROOPS_GNSSPARAMETRIZATIONIONOSPHEREMAP__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationIonosphereMap = R"(
\subsection{IonosphereMap}\label{gnssParametrizationType:ionosphereMap}
TODO: reading and writing ionosphere maps not implemented yet.
% A priori ionopshere model can be provided with \configFile{inputfileMap}{gnssIonosphereMaps}.

The ionosphere is parametrized as a \configClass{temporal}{parametrizationTemporalType}ly changing
(e.g. hourly linear splines)
spherical harmonics expansion
\begin{equation}
  VTEC(\lambda,\theta,t) = \sum_{n=0}^{n_{max}} \sum_{m=0}^n c_{nm}(t)C_{nm}(\lambda,\theta)+s_{nm}(t)S_{nm}(\lambda,\theta)
\end{equation}
up to \config{maxDegree}=\verb|15| in a solar-geomagentic frame defined
by \configClass{magnetosphere}{magnetosphereType}.
The VTEC values are mapped to STEC values via eq.~\eqref{gnssParametrizationType:IonosphereVTEC:STEC}.

Local and short-term scintillations can be considered by adding constrained
\configClass{parametrization:ionosphereSTEC}{gnssParametrizationType:ionosphereSTEC}.
To account for signal biases add
\configClass{parametrization:tecBiases}{gnssParametrizationType:tecBiases}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/magnetosphere/magnetosphere.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "gnss/gnss.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief IonosphereMap.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationIonosphereMap : public GnssParametrizationBase
{
  Gnss                           *gnss;
  std::string                     name;
  GnssTransceiverSelectorPtr      selectReceivers;
  std::vector<Byte>               selectedReceivers;
  FileName                        fileNameIn, fileNameOut;
  UInt                            maxDegree;
  ParametrizationTemporalPtr      temporal;
  Double                          mapR, mapH, mapAlpha;
  MagnetospherePtr                magnetosphere;
  std::vector<Vector>             x;
  std::vector<GnssParameterIndex> index;

  Vector3d intersection(const Vector3d &posRecv, const Vector3d &posTrans, Angle elevation) const;
  Double   mapping(Angle elevation) const;

public:
  GnssParametrizationIonosphereMap(Config &config);

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
