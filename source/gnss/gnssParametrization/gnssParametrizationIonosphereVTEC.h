/***********************************************/
/**
* @file gnssParametrizationIonosphereVTEC.h
*
* @brief IonosphereVTEC.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONIONOSPHEREVTEC__
#define __GROOPS_GNSSPARAMETRIZATIONIONOSPHEREVTEC__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationIonosphereVTEC = R"(
\subsection{IonosphereVTEC}\label{gnssParametrizationType:ionosphereVTEC}
The influence of the ionosphere is modelled by a VTEC parameter (vertical total electron content)
in terms of $[TECU]$ for every selected receiver at each epoch. Optionally, VTEC gradients in the
North (x) and East (y) direction can be estimated via \configClass{gradient}{parametrizationTemporalType}.
The slant TEC is computed based on the VTEC and the optional North and East gradients $\Delta V_x$ and $\Delta V_y$
using the elevation-dependent Modified Single-Layer Model (MSLM) mapping function
\begin{equation}\label{gnssParametrizationType:IonosphereVTEC:STEC}
  STEC = \frac{VTEC + \cos(A) \Delta V_x + \sin(A) \Delta V_y}{\cos z'}
  \qquad\text{with}\qquad
  \sin z'= \left(\frac{R}{R+H}\right)\sin\left(\alpha(\pi/2-E)\right)
\end{equation}
inserted into eq.~\eqref{gnssParametrizationType:IonosphereSTEC:STEC},
where $A$ is the azimuth angle and $E$ is the elevation angle.

The result is written as a \file{times series file}{instrument} at epochs with observations
depending on \configClass{GnssProcessing:processingStep:selectEpochs}{gnssProcessingStepType:selectEpochs}.

This class provides a simplified model of the ionosphere for single receivers
and enables the separation of the TEC and signal biases, meaning
\configClass{parametrization:tecBiases}{gnssParametrizationType:tecBiases} becomes estimable.
Local and short-term scintillations should be considered by adding loosely constrained
\configClass{parametrization:ionosphereSTEC}{gnssParametrizationType:ionosphereSTEC}.

The \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|<station>:VTEC::<time>|,
\item \verb|<station>:VTECGradient.x:<temporal>:<interval>|,
\item \verb|<station>:VTECGradient.y:<temporal>:<interval>|.
\end{itemize}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief IonosphereVTEC.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationIonosphereVTEC : public GnssParametrizationBase
{
  Gnss                                        *gnss;
  std::string                                  name;
  PlatformSelectorPtr                          selectReceivers;
  std::vector<Byte>                            selectedReceivers;
  FileName                                     fileNameVTEC;
  Double                                       mapR, mapH, mapAlpha;
  std::vector<std::vector<GnssParameterIndex>> indexVTEC; // for each receiver, for each epoch
  std::vector<std::vector<Double>>             VTEC;

  ParametrizationTemporalPtr                   parametrizationGradient;
  std::vector<GnssParameterIndex>              indexGradient; // for each receiver
  std::vector<Vector>                          xGradient;
  std::vector<std::vector<Double>>             gradientX, gradientY;

  Double mapping(Angle elevation) const;

public:
  GnssParametrizationIonosphereVTEC(Config &config);

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
