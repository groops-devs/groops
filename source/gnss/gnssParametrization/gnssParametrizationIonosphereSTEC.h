/***********************************************/
/**
* @file gnssParametrizationIonosphereSTEC.h
*
* @brief IonosphereSTEC.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONIONOSPHERESTEC__
#define __GROOPS_GNSSPARAMETRIZATIONIONOSPHERESTEC__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationIonosphereSTEC = R"(
\subsection{IonosphereSTEC}\label{gnssParametrizationType:ionosphereSTEC}
The influence of the ionosphere is modelled by a STEC parameter (slant total electron content)
in terms of $[TECU]$ between each transmitter and receiver at each epoch. These parameters are pre-eliminated
from the observation equations before accumulating the normal equations.
This is similar to using the ionosphere-free linear combination as observations
but only one STEC parameter is needed for an arbitrary number of observation types.

The influence on the code and phase observation is modeled as
\begin{equation}\label{gnssParametrizationType:IonosphereSTEC:STEC}
\begin{split}
\text{ionosphere}([C\nu], STEC) &=  \frac{40.3}{f_{\nu}^2}STEC + \frac{7525\M b^T\M k}{f_{\nu}^3}STEC +  \frac{r}{f_{\nu}^4}STEC^2 \\
\text{ionosphere}([L\nu], STEC) &= -\frac{40.3}{f_{\nu}^2}STEC - \frac{7525\M b^T\M k}{2f_{\nu}^3}STEC - \frac{r}{3f_{\nu}^4}STEC^2 + \text{bending}(E)STEC^2
\end{split}
\end{equation}
The second order term depends on the \configClass{magnetosphere}{magnetosphereType} $\M b$
and the direction of the signal $\M k$.

If further information about the ionosphere is available
(in the form of a prior model or as additional parametrizations
such as \configClass{parametrization:ionosphereMap}{gnssParametrizationType:ionosphereMap} or
\configClass{parametrization:ionosphereVTEC}{gnssParametrizationType:ionosphereVTEC}) the STEC
parameters describe local and short–term scintillations. The STEC parameters are estimated
as additions to the model and it is advised to constrain them towards zero
with a standard deviation of \config{sigmaSTEC}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "classes/magnetosphere/magnetosphere.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief IonosphereSTEC.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationIonosphereSTEC : public GnssParametrizationBase
{
  Gnss                 *gnss;
  std::string           name, nameConstraint;
  Bool                  apply1stOrder, apply2ndOrder, apply3rdOrder, applyBending;
  MagnetospherePtr      magnetosphere;
  Bool                  estimateSTEC, applyConstraint;
  ExpressionVariablePtr exprSigmaSTEC;
  Double                sigmaSTEC; // 0: unconstrained, <0: expression

public:
  GnssParametrizationIonosphereSTEC(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void observationCorrections(GnssObservationEquation &eqn) const override;
  void initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
};

/***********************************************/

#endif
