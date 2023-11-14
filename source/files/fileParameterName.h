/***********************************************/
/**
* @file fileParameterName.h
*
* @brief Read/write ParameterName.
*
* @author Torsten Mayer-Guerr
* @date 2017-07-31
*
*/
/***********************************************/

#ifndef __GROOPS_FILEPARAMETERNAME__
#define __GROOPS_FILEPARAMETERNAME__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_ParameterName
static const char *docstringParameterName = R"(
Name of parameters of a system of \file{normal equations}{normalEquation} or \file{solution vector}{matrix}.

A parameter name is a string \verb|<object>:<type>:<temporal>:<interval>| containg four parts divided by \verb|:|
\begin{enumerate}
\item object: Object this parameter refers to, e.g. \verb|graceA|, \verb|G023|, \verb|earth|, \ldots
\item type: Type of this parameter, e.g. \verb|accBias|, \verb|position.x|, \ldots
\item temporal: Temporal representation of this parameter, e.g. \verb|trend|, \verb|polynomial.degree1|, \ldots
\item interval: Interval/epoch this parameter represents, e.g. \verb|2017-01-01_00-00-00_2017-01-02_00-00-00|, \verb|2018-01-01_00-00-00|.
\end{enumerate}
In the documentation a star (\verb|*|) in the name means this part is untouched and useally set by other classes.
Times are written as \verb|yyyy-mm-dd_hh-mm-ss| and intervals (if not empty) as \verb|<timeStart>_<timeEnd>|.

See \program{ParameterNamesCreate}.

\begin{verbatim}
groops parameterName version=20200123
# object:type:temporal:interval
# =============================
      10080 # number of parameters
 karr:position.x::2018-06-01_00-00-00_2018-06-02_00-00-00
 karr:position.y::2018-06-01_00-00-00_2018-06-02_00-00-00
 karr:position.z::2018-06-01_00-00-00_2018-06-02_00-00-00
 karr:troposphereWet:spline.n1:2018-06-01_00-00-00_2018-06-01_02-00-00
 karr:troposphereWet:spline.n1:2018-06-01_00-00-00_2018-06-01_04-00-00
 karr:troposphereWet:spline.n1:2018-06-01_02-00-00_2018-06-01_06-00-00
 karr:troposphereWet:spline.n1:2018-06-01_04-00-00_2018-06-01_08-00-00
 karr:troposphereWet:spline.n1:2018-06-01_06-00-00_2018-06-01_10-00-00
 karr:troposphereWet:spline.n1:2018-06-01_08-00-00_2018-06-01_12-00-00
 karr:troposphereWet:spline.n1:2018-06-01_10-00-00_2018-06-01_14-00-00
 karr:troposphereWet:spline.n1:2018-06-01_12-00-00_2018-06-01_16-00-00
 karr:troposphereWet:spline.n1:2018-06-01_14-00-00_2018-06-01_18-00-00
 karr:troposphereWet:spline.n1:2018-06-01_16-00-00_2018-06-01_20-00-00
 karr:troposphereWet:spline.n1:2018-06-01_18-00-00_2018-06-01_22-00-00
 karr:troposphereWet:spline.n1:2018-06-01_20-00-00_2018-06-02_00-00-00
 karr:troposphereWet:spline.n1:2018-06-01_22-00-00_2018-06-02_00-00-00
 karr:troposphereGradient.x:spline.n1:2018-06-01_00-00-00_2018-06-02_00-00-00
 karr:troposphereGradient.y:spline.n1:2018-06-01_00-00-00_2018-06-02_00-00-00
 karr:troposphereGradient.x:spline.n1:2018-06-01_00-00-00_2018-06-02_00-00-00
 karr:troposphereGradient.y:spline.n1:2018-06-01_00-00-00_2018-06-02_00-00-00
 karr:signalBias01(+1.00L1CG**)::
 karr:signalBias02(+1.00L2WG**)::
 karr:signalBias03(+1.00L2XG**)::
 G01:solarRadiationPressure.ECOM.D0::
 G01:solarRadiationPressure.ECOM.DC2::
 G01:solarRadiationPressure.ECOM.DS2::
 G01:solarRadiationPressure.ECOM.Y0::
 G01:solarRadiationPressure.ECOM.B0::
 G01:solarRadiationPressure.ECOM.BC1::
 G01:solarRadiationPressure.ECOM.BS1::
 G01:stochasticPulse.x::2018-06-01_12-00-00
 G01:stochasticPulse.y::2018-06-01_12-00-00
 G01:stochasticPulse.z::2018-06-01_12-00-00
 G01:arc0.position0.x::
 G01:arc0.position0.y::
 G01:arc0.position0.z::
 G01:arc0.velocity0.x::
 G01:arc0.velocity0.y::
 G01:arc0.velocity0.z::
 G01:signalBias01(-1.00C1CG01)::
 G01:signalBias02(+1.00L1*G01)::
 G01:signalBias03(+1.00L2*G01)::
\end{verbatim}
)";
#endif

/***********************************************/

#include "base/exception.h"
#include "base/parameterName.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_PARAMETERNAME_TYPE = "parameterName";

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const ParameterName &x);
template<> void load(InArchive  &ar, ParameterName &x);

/** @brief Write into a ParameterName file. */
void writeFileParameterName(const FileName &fileName, const std::vector<ParameterName> &x);

/** @brief Read from a ParameterName file. */
void readFileParameterName(const FileName &fileName, std::vector<ParameterName> &x);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
