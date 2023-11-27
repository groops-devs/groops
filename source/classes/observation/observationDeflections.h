/***********************************************/
/**
* @file observationDeflections.h
*
* @brief point measareuments of xi and eta.
*
* @author Christian Pock
* @date 2012-05-30
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONDEFLECTIONS__
#define __GROOPS_OBSERVATIONDEFLECTIONS__

// Latex documentation
#ifdef DOCSTRING_Observation
static const char *docstringObservationDeflections = R"(
\subsection{Deflections}\label{observationType:deflections}
The gravity field parametrized by \configClass{parametrizationGravity}{parametrizationGravityType}
is estimated from deflections of the vertical measurements.
A \configClass{referencefield}{gravityfieldType} can be reduced beforehand.

The observations $\xi$ in north direction and $\eta$ in east direction
at given positions are calculated from
\configFile{inputfileGriddedData}{griddedData}.
The input columns are enumerated by \verb|data0|,~\verb|data1|,~\ldots,
see~\reference{dataVariables}{general.parser:dataVariables}.

The ellipsoid parameters \config{R} and \config{inverseFlattening} are used
to define the local normal direction.

The observations can be divided into small blocks for parallelization.
With \config{blockingSize} set the maximum count of observations in each block.

The following parameters with \file{parameter names}{parameterName} are set up:
\verb|*:<parametrizationGravity>:*:*|.
)";
#endif

/***********************************************/

#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/observation/observation.h"

/***** CLASS ***********************************/

/**
* @brief point measurements of xi and eta
* @see Observation
*/
class ObservationDeflections : public Observation
{
  std::vector<Vector3d>     points;
  std::vector<Double>       xi, eta;
  std::vector<Double>       sigmasXi, sigmasEta;
  GravityfieldPtr           referencefield;
  ParametrizationGravityPtr parametrization;
  Time                      time;
  Ellipsoid                 ellipsoid;
  UInt                      obsPerArc;

  UInt idx(UInt arc, UInt obs)    const {return arc*obsPerArc+obs;}
  UInt observationCount(UInt arc) const {return (arc<arcCount()-1) ? (obsPerArc) : (points.size()-arc*obsPerArc);}

public:
  ObservationDeflections(Config &config);
 ~ObservationDeflections() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override {return parametrization->setInterval(timeStart, timeEnd);}
  UInt parameterCount()        const override {return parametrization->parameterCount();}
  UInt gravityParameterCount() const override {return parametrization->parameterCount();}
  UInt rightSideCount()        const override {return 1;}
  UInt arcCount()              const override {return (points.size()+obsPerArc-1)/obsPerArc;}
  void parameterName(std::vector<ParameterName> &name) const override {parametrization->parameterName(name);}

  void observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B) override;
};

/***********************************************/

#endif /* __GROOPS_OBSERVATION__ */
