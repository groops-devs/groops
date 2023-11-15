/***********************************************/
/**
* @file parametrizationGnssAntennaRadialBasis.h
*
* @brief Parametrization of antenna center variations.
*
* @author Norbert Zehentner
* @date 2018-07-27
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGNSSANTENNARADIALBASIS__
#define __GROOPS_PARAMETRIZATIONGNSSANTENNARADIALBASIS__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGnssAntenna
static const char *docstringParametrizationGnssAntennaRadialBasis = R"(
\subsection{RadialBasis}
Parametrization of antenna center variations with radial basis functions
\begin{equation}
  ACV(\M x(A, E)) = \sum_i a_i \Phi(\M x\cdot\M x_i)
\end{equation}
where $a_i$ in $[m]$ the coefficients which has to be estimated and $\Phi$ are the basis
functions
\begin{equation}
  \Phi(\cos\psi) = \sum_n \sqrt{2n+1}P_n(\cos\psi).
\end{equation}

The \file{parameter names}{parameterName} are
\verb|*:antennaCenterVariations.radialBasis.<index>.<total count>:*:*|.

\fig{!hb}{0.4}{parametrizationGnssAntennaRadialBasis}{fig:parametrizationGnssAntennaRadialBasis}{Nodal points of the basis functions
using a Reuter grid for transmitting satellites (view angle of 18 deg). The red line indicates the view angle of 14 deg of ground stations.}
)";
#endif

/***********************************************/

#include "base/legendrePolynomial.h"
#include "classes/grid/grid.h"
#include "parametrizationGnssAntenna.h"

/***** CLASS ***********************************/

class ParametrizationGnssAntennaRadialBasis : public ParametrizationGnssAntennaBase
{
  std::vector<Vector3d> point;
  Vector                coeff;

public:
  ParametrizationGnssAntennaRadialBasis(Config &config);
  UInt parameterCount() const override {return point.size();}
  void parameterName(std::vector<ParameterName> &name) const override;
  void designMatrix(Angle azimut, Angle elevation, MatrixSliceRef A) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationGnssAntennaRadialBasis::ParametrizationGnssAntennaRadialBasis(Config &config)
{
  try
  {
    GridPtr grid;
    UInt    minDegree, maxDegree;

    readConfig(config, "grid",      grid,      Config::MUSTSET, "",  "nodal points of shannon kernels");
    readConfig(config, "minDegree", minDegree, Config::MUSTSET, "2", "min degree of shannon kernel");
    readConfig(config, "maxDegree", maxDegree, Config::MUSTSET, "",  "max degree of shannon kernel");
    if(isCreateSchema(config)) return;

    // init nodal points
    point = grid->points();
    for(UInt i=0; i<point.size(); i++)
      point.at(i).normalize();

    // init kernel coefficients
    coeff = Vector(maxDegree+1);
    for(UInt n=minDegree; n<=maxDegree; n++)
      coeff(n) = 1.0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationGnssAntennaRadialBasis::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt i=0; i<point.size(); i++)
      name.push_back(ParameterName("", "antennaCenterVariations.radialBasis."+i%"%i"s+"."+point.size()%"%i"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationGnssAntennaRadialBasis::designMatrix(Angle azimut, Angle elevation, MatrixSliceRef A)
{
  try
  {
    const Vector3d e12 = polar(azimut, elevation, 1.);
    for(UInt i=0; i<point.size(); i++)
      A(0,i) = LegendrePolynomial::sum(inner(e12, point.at(i)), coeff, coeff.size()-1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS___ */
