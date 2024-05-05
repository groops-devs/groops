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
Apriori VTEC maps can be removed from the observations with
\configFile{inputfileGriddedDataTimeSeries}{griddedDataTimeSeries}
(e.g. from \program{GnssIonex2GriddedDataTimeSeries}).

The ionosphere is parametrized in terms of $[TECU]$ in a single layer sphere with
\config{radiusIonosphericLayer} as a \configClass{temporal}{parametrizationTemporalType}ly
changing (e.g. hourly linear splines) spherical harmonics expansion
\begin{equation}
  VTEC(\lambda,\theta,t) = \sum_{n=0}^{n_{max}} \sum_{m=0}^n c_{nm}(t)C_{nm}(\lambda,\theta)+s_{nm}(t)S_{nm}(\lambda,\theta)
\end{equation}
up to \config{maxDegree}=\verb|15| in a solar-geomagentic frame defined
by \configClass{magnetosphere}{magnetosphereType}. The VTEC values are mapped to STEC values
in the observation equations via eq.~\eqref{gnssParametrizationType:IonosphereVTEC:STEC}.

The estimated VTEC inclusive the apriori \configFile{inputfileGriddedDataTimeSeries}{griddedDataTimeSeries}
can be written to \configFile{outputfileGriddedDataTimeSeries}{griddedDataTimeSeries}
evaluated at \configClass{outputGrid}{gridType} and \configClass{outputTimeSeries}{timeSeriesType}.

Local and short-term scintillations should be considered by adding constrained
\configClass{parametrization:ionosphereSTEC}{gnssParametrizationType:ionosphereSTEC}.
To account for signal biases add
\configClass{parametrization:tecBiases}{gnssParametrizationType:tecBiases}.

The \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|VTEC:sphericalHarmonics.c_<degree>_<order>:<temporal>:<interval>|,
\item \verb|VTEC:sphericalHarmonics.s_<degree>_<order>:<temporal>:<interval>|.
\end{itemize}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/griddedData.h"
#include "config/config.h"
#include "classes/magnetosphere/magnetosphere.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief IonosphereMap.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationIonosphereMap : public GnssParametrizationBase
{
  Gnss                           *gnss;
  std::string                     name;
  PlatformSelectorPtr             selectReceivers;
  std::vector<Byte>               selectedReceivers;
  FileName                        fileNameIn, fileNameOut;
  GriddedData                     gridOut;
  std::vector<Time>               timesOut;
  UInt                            maxDegree;
  ParametrizationTemporalPtr      temporal;
  Double                          radiusIono, mapR, mapH, mapAlpha;
  MagnetospherePtr                magnetosphere;
  std::vector<Vector>             x;
  std::vector<GnssParameterIndex> index;

  Vector3d intersection(const Double radiusIono, const Vector3d &posRecv, const Vector3d &posTrans) const;
  Double   mapping(Angle elevation) const;
  Double   sphericalHarmonicSynthesis(const Vector3d &point, const Vector &x) const;

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
