/***********************************************/
/**
* @file gnssParametrizationReceiverAntennas.h
*
* @brief Antenna center variations.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONRECEIVERANTENNAS__
#define __GROOPS_GNSSPARAMETRIZATIONRECEIVERANTENNAS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationReceiverAntennas = R"(
\subsection{ReceiverAntennas}\label{gnssParametrizationType:receiverAntennas}
This class is for parametrization the antenna for their antenna center offsets (ACO) and
antenna center variations (ACV) by \configClass{antennaCenterVariations}{parametrizationGnssAntennaType}.
The receivers to be estimated can be selected by \configClass{selectReceivers}{gnssTransceiverSelectorType}.

The amount of patterns to be estimated is configurable with a list of \configClass{patternTypes}{gnssType}.
For each added \configClass{patternTypes}{gnssType} a set of parameters will be evaluated. The observations
will be assigned to the first \configClass{patternTypes}{gnssType} that matches their own.
E.g. having the patterns: \verb|***G| and \verb|L1*| would lead to all GPS observations be assigned
to the observation equations of the first pattern. The patterntype \verb|L1*| would then consist
of all other GNSS L1 phase observations. \config{addNonMatchingTypes} will, if activated, create automatically patterns
for \configClass{observations}{gnssType} that are not selected within the list \configClass{patternTypes}{gnssType}.
Furthermore, it is possible to group same antenna build types from different receivers by \config{groupAntennas}.
The grouping by same antenna build ignores antenna serial numbers.

To estimate the antenna variation parameters, a longer period of observations might be necessary
for accurate estimations. Hence one should use this parametrization by
accumulating normal equations from several epochs.
This can be accomplished as the last steps in the \configClass{processing steps}{gnssProcessingStepType}
 by adding \configClass{ReceiverAntennas}{gnssParametrizationType:receiverAntennas}
to current selected parameters with \configClass{GnssProcessing:processingStep:selectParametrizations}{gnssProcessingStepType:selectParametrizations}
and write the normal equation matrix with \configClass{GnssProcessing:processingStep:writeNormalEquations}{gnssProcessingStepType:writeNormalEquations}.
The written normal equations can then be accumulated with \program{NormalsAccumulate} and solved by \program{NormalsSolverVCE}.
Further, one should apply constraints to the normal equations by \program{GnssAntennaNormalsConstraint} since the estimation
 of ACO and ACV can lead to rank deficiencies in the normal equation matrix.
Last the solved normal equation can be parsed to a \file{antenna definition file}{gnssAntennaDefinition}
 with the program \program{ParameterVector2GnssAntennaDefinition}.

As example refering to the cookbook \reference{GNSS satellite orbit determination and station network analysis}{cookbook.gnssNetwork},
one could add additionally \configClass{receiverAntennas}{gnssParametrizationType:receiverAntennas} as parametrization.
Since the estimations are done on a daily basis for each receiver we add an additional
\configClass{selectParametrizations}{gnssProcessingStepType:selectParametrizations} which
disables \verb|parameter.receiverAntenna|. After all stations are processed together with all parameters, one
adds \verb|parameter.receiverAntenna| with \configClass{selectParametrizations}{gnssProcessingStepType:selectParametrizations}
 to the current selected parametrizations.
The last \configClass{processingStep}{gnssProcessingStepType} is \configClass{GnssProcessing:processingStep:writeNormalEquations}{gnssProcessingStepType:writeNormalEquations}
to write the daily normal equations including the parametrization \configClass{receiverAntennas}{gnssParametrizationType:receiverAntennas} into files.
These normal equation files are then processed with the programms:

\begin{itemize}
  \item \program{NormalsAccumulate}: accumulates normal equations.
  \item \program{GnssAntennaNormalsConstraint}: apply constraint to the normal equations.
  \item \program{NormalsSolverVCE}: solves the normal equations.
  \item \program{ParameterVector2GnssAntennaDefinition}: writes the solution into a \file{antenna definition file}{gnssAntennaDefinition}
\end{itemize}

Note that the apriori value $\M x_0$ for this parametrization is always zero and never updated
according to eq.~\eqref{gnssParametrizationType:update}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Antenna center variations.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationReceiverAntennas : public GnssParametrizationBase
{
  Gnss                                        *gnss;
  std::string                                  name;
  GnssTransceiverSelectorPtr                   selectReceivers;
  ParametrizationGnssAntennaPtr                parametrization;
  std::vector<GnssType>                        typesPattern;
  Bool                                         addNonMatchingTypes;
  Bool                                         ignoreSerial;
  std::vector<UInt>                            receiver2antenna;
  std::vector<std::vector<GnssParameterIndex>> index; // for each antenna and pattern
  std::vector<std::vector<GnssType>>           types; // for each antenna and pattern

public:
  GnssParametrizationReceiverAntennas(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
};

/***********************************************/

#endif
