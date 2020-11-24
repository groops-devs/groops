/***********************************************/
/**
* @file observationTerrestrial.h
*
* @brief Terrestrial observations (point measurements).
* Observed function values of Gravityfield.
*
* @author Torsten Mayer-Guerr
* @date 2005-03-22
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONTERRESTRIAL__
#define __GROOPS_OBSERVATIONTERRESTRIAL__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationTerrestrial = R"(
\subsection{Terrestrial}\label{observationType:terrestrial}
The gravity field is estimated from point wise measurements.
The gravity field parametrization is given by \configClass{parametrizationGravity}{parametrizationGravityType}.
There is no need to have the data regular distributed or given on a sphere or ellipsoid.
The type of the gridded data (e.g gravity anomalies or geoid heights)
must be set with \configClass{kernel}{kernelType}.
A \configClass{referencefield}{gravityfieldType} can be reduced beforehand.

The observations at given positions are calculated from
\configFile{inputfileGriddedData}{griddedData}.
The input columns are enumerated by \verb|data0|,~\verb|data1|,~\ldots,
see~\reference{dataVariables}{general.parser:dataVariables}.

The observations can be divided into small blocks for parallelization.
With \config{blockingSize} set the maximum count of observations in each block.
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/observation/observation.h"

/***** CLASS ***********************************/

/** @brief Terrestrial observations (point measurements).
* @ingroup observationGroup
* Observed function values of Gravityfield.
* @see Observation */
class ObservationTerrestrial : public Observation
{
  std::vector<Vector3d>     points;
  std::vector<Double>       values;
  std::vector<Double>       sigmas;
  GravityfieldPtr           referencefield;
  KernelPtr                 kernel;
  ParametrizationGravityPtr parametrization;
  Time                      time;
  UInt                      obsPerArc;

  UInt idx(UInt arc, UInt obs)    const {return arc*obsPerArc+obs;}
  UInt observationCount(UInt arc) const {return (arc<arcCount()-1) ? (obsPerArc) : (points.size()-arc*obsPerArc);}

public:
  ObservationTerrestrial(Config &config);
 ~ObservationTerrestrial() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override {return parametrization->setInterval(timeStart, timeEnd);}
  UInt parameterCount()          const override {return parametrization->parameterCount();}
  UInt gravityParameterCount()   const override {return parametrization->parameterCount();}
  UInt rightSideCount()          const override {return 1;}
  UInt arcCount()                const override {return (points.size()+obsPerArc-1)/obsPerArc;}
  void parameterName(std::vector<ParameterName> &name) const override {parametrization->parameterName(name);}

  void observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B) override;
};

/***********************************************/

#endif /* __GROOPS_OBSERVATION__ */
