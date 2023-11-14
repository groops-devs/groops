/***********************************************/
/**
* @file parametrizationGnssAntennaSphericalHarmonics.h
*
* @brief Parametrization of antenna center variations.
*
* @author Torsten Mayer-Guerr
* @date 2012-12-28
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGNSSANTENNASPHERICALHARMONICS__
#define __GROOPS_PARAMETRIZATIONGNSSANTENNASPHERICALHARMONICS__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGnssAntenna
static const char *docstringParametrizationGnssAntennaSphericalHarmonics = R"(
\subsection{SphericalHarmonics}
Parametrization of antenna center variations in $[m]$ in terms of spherical harmonics.
As usually only data above the horizon are observed only the even spherical harmonics
(degree/order $m+n$ even), which are symmetric to the equator, are setup.

The total count of parameters is $((n_{max}+1)(n_{max}+2)-n_{min}(n_{min}+1)/2$ and
the \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|*:antennaCenterVariations.sphericalHarmonics.c_<degree>_<order>:*:*|,
\item \verb|*:antennaCenterVariations.sphericalHarmonics.s_<degree>_<order>:*:*|.
\end{itemize}

)";
#endif

/***********************************************/

#include "base/sphericalHarmonics.h"
#include "parametrizationGnssAntenna.h"

/***** CLASS ***********************************/

class ParametrizationGnssAntennaSphericalHarmonics : public ParametrizationGnssAntennaBase
{
  UInt minDegree, maxDegree;
  UInt parameterCount_;

public:
  ParametrizationGnssAntennaSphericalHarmonics(Config &config);
  UInt parameterCount() const override {return parameterCount_;}
  void parameterName(std::vector<ParameterName> &name) const override;
  void designMatrix(Angle azimut, Angle elevation, MatrixSliceRef A) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationGnssAntennaSphericalHarmonics::ParametrizationGnssAntennaSphericalHarmonics(Config &config)
{
  try
  {
    readConfig(config, "minDegree", minDegree, Config::MUSTSET, "2", "min degree");
    readConfig(config, "maxDegree", maxDegree, Config::MUSTSET, "",  "max degree");;
    if(isCreateSchema(config)) return;

    parameterCount_ = ((maxDegree+1)*(maxDegree+2)-minDegree*(minDegree+1))/2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationGnssAntennaSphericalHarmonics::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      if(n%2 == 0) // even?
        name.push_back(ParameterName("", "antennaCenterVariations.sphericalHarmonics.c_"+n%"%i"s+"_0"));
      for(UInt m=2-n%2; m<=n; m+=2)
      {
        name.push_back(ParameterName("", "antennaCenterVariations.sphericalHarmonics.c_"+n%"%i"s+"_"+m%"%i"s));
        name.push_back(ParameterName("", "antennaCenterVariations.sphericalHarmonics.s_"+n%"%i"s+"_"+m%"%i"s));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationGnssAntennaSphericalHarmonics::designMatrix(Angle azimut, Angle elevation, MatrixSliceRef A)
{
  try
  {
    const Vector3d e12 = polar(azimut, elevation, 1.);
    Matrix Cnm, Snm;
    SphericalHarmonics::CnmSnm(e12, maxDegree, Cnm, Snm);
    // only the functions which are symmetric to equator
    UInt idx = 0;
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      if((n%2) == 0) // even?
        A(0,idx++) = Cnm(n,0);
      for(UInt m=2-n%2; m<=n; m+=2)
      {
        A(0,idx++) = Cnm(n,m);
        A(0,idx++) = Snm(n,m);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS___ */
