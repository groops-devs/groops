/***********************************************/
/**
* @file parametrizationGnssAntennaCenter.h
*
* @brief Parametrization of antenna center variations.
*
* @author Torsten Mayer-Guerr
* @date 2018-07-27
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGNSSANTENNACENTER__
#define __GROOPS_PARAMETRIZATIONGNSSANTENNACENTER__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGnssAntenna
static const char *docstringParametrizationGnssAntennaCenter = R"(
\subsection{Center}\label{parametrizationGnssAntennaType:center}
Antenna center or, if setup for a specific \configClass{gnssType}{gnssType},
phase/code center offset (e.g. \verb|*1*G| for GPS L1 phase center offset) in $[m]$.

The \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|*:antennaCenter.x:*:*|,
\item \verb|*:antennaCenter.y:*:*|,
\item \verb|*:antennaCenter.z:*:*|.
\end{itemize}
)";
#endif

/***********************************************/

#include "parametrizationGnssAntenna.h"

/***** CLASS ***********************************/

class ParametrizationGnssAntennaCenter : public ParametrizationGnssAntennaBase
{
  Bool estimateX, estimateY, estimateZ;

public:
  ParametrizationGnssAntennaCenter(Config &config);
  UInt parameterCount() const override {return estimateX+estimateY+estimateZ;}
  void parameterName(std::vector<ParameterName> &name) const override;
  void designMatrix(Angle azimut, Angle elevation, MatrixSliceRef A) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationGnssAntennaCenter::ParametrizationGnssAntennaCenter(Config &config)
{
  try
  {
    readConfig(config, "estimateX", estimateX,  Config::DEFAULT,  "1", "");
    readConfig(config, "estimateY", estimateY,  Config::DEFAULT,  "1", "");
    readConfig(config, "estimateZ", estimateZ,  Config::DEFAULT,  "1", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationGnssAntennaCenter::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    if(estimateX) name.push_back(ParameterName("", "antennaCenter.x"));
    if(estimateY) name.push_back(ParameterName("", "antennaCenter.y"));
    if(estimateZ) name.push_back(ParameterName("", "antennaCenter.z"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationGnssAntennaCenter::designMatrix(Angle azimut, Angle elevation, MatrixSliceRef A)
{
  try
  {
    UInt idxAxis = 0;
    if(estimateX) A(0, idxAxis++) = -std::cos(elevation) * std::cos(azimut);
    if(estimateY) A(0, idxAxis++) = -std::cos(elevation) * std::sin(azimut);
    if(estimateZ) A(0, idxAxis++) = -std::sin(elevation);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS___ */
