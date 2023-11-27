/***********************************************/
/**
* @file observationStationLoading.h
*
* @brief Loading from station observations.
*
* @author Torsten Mayer-Guerr
* @date 2014-07-17
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONSTATIONLOADING__
#define __GROOPS_OBSERVATIONSTATIONLOADING__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationStationLoading = R"(
\subsection{StationLoading}\label{observationType:stationLoading}
Observation equations for displacements of a list of stations
due to the effect of time variable loading masses. The displacement~$\M u$ of a station is calculated according to
\begin{equation}
\M u(\M r) = \frac{1}{\gamma}\sum_{n=0}^\infty \left[\frac{h_n}{1+k_n}V_n(\M r)\,\M e_{up}
+ R\frac{l_n}{1+k_n}\left(
 \frac{\partial V_n(\M r)}{\partial \M e_{north}}\M e_{north}
+\frac{\partial V_n(\M r)}{\partial \M e_{east}} \M e_{east}\right)\right],
\end{equation}
where $\gamma$ is the normal gravity, the load Love and Shida numbers $h_n,l_n$ are given by
\configFile{inputfileDeformationLoadLoveNumber}{matrix} and the load Love numbers $k_n$ are given by
\configFile{inputfilePotentialLoadLoveNumber}{matrix}.
The $V_n$ are the spherical harmonics expansion of degree $n$ of the full time variable
gravitational potential (potential of the loading mass + deformation potential)
parametrized by \configClass{parametrizationGravity}{parametrizationGravityType}.
Additional parameters can be setup to estimate the realization of the reference frame
of the station coordinates (\config{estimateTranslation},
\config{estimateRotation}, and \config{estimateScale}).

The observations at stations coordinates are calculated from
\configFile{inputfileGriddedData}{griddedData}.
The input columns are enumerated by \verb|data0|,~\verb|data1|,~\ldots,
see~\reference{dataVariables}{general.parser:dataVariables}.

The ellipsoid parameters \config{R} and \config{inverseFlattening} are used
to define the local frame (north, east, up).

The following parameters with \file{parameter names}{parameterName} are set up:
\begin{itemize}
\item \verb|*:<parametrizationGravity>:*:*|,
\item \verb|*:translation.x:*:*|,
\item \verb|*:translation.y:*:*|,
\item \verb|*:translation.z:*:*|,
\item \verb|*:scale:*:*|,
\item \verb|*:rotation.x:*:*|,
\item \verb|*:rotation.y:*:*|,
\item \verb|*:rotation.z:*:*|.
\end{itemize}

See also \program{Gravityfield2DisplacementTimeSeries}.

Reference:
Rietbroek (2014): Retrieval of Sea Level and Surface Loading Variations from Geodetic Observations
and Model Simulations: an Integrated Approach, Bonn, 2014. - Dissertation,
\url{https://nbn-resolving.org/urn:nbn:de:hbz:5n-35460}
)";
#endif

/***********************************************/

#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/observation/observation.h"

/***** CLASS ***********************************/

/** @brief Loading from station observations.
* @ingroup observationGroup
* Observed function values of Gravityfield.
* @see Observation */
class ObservationStationLoading : public Observation
{
  std::vector<Vector3d>     points;
  std::vector<Double>       valuesNorth, valuesEast, valuesUp;
  std::vector<Double>       sigmasNorth, sigmasEast, sigmasUp;
  Bool                      inGlobalFrame;
  Time                      time;
  Vector                    hn, ln;
  GravityfieldPtr           referencefield;
  ParametrizationGravityPtr parametrization;
  Bool                      estimateTranslation, estimateScale, estimateRotation;
  Ellipsoid                 ellipsoid;

public:
  ObservationStationLoading(Config &config);
 ~ObservationStationLoading() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override {return parametrization->setInterval(timeStart, timeEnd);}
  UInt parameterCount()          const override {return gravityParameterCount() + 3*estimateTranslation + estimateScale + 3*estimateRotation;}
  UInt gravityParameterCount()   const override {return parametrization->parameterCount();}
  UInt rightSideCount()          const override {return 1;}
  UInt arcCount()                const override {return 1;}
  void parameterName(std::vector<ParameterName> &name) const override;

  void observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B) override;
};

/***********************************************/

#endif /* __GROOPS_OBSERVATION__ */
