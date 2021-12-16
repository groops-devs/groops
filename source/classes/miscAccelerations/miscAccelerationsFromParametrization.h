/***********************************************/
/**
* @file miscAccelerationsFromParametrization.h
*
* @brief .
* @see MiscAccelerations
*
* @author Torsten Mayer-Guerr
* @date 2021-11-10
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONSFROMPARAMETRIZATION__
#define __GROOPS_MISCACCELERATIONSFROMPARAMETRIZATION__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerationsFromParametrization = R"(
\subsection{FromParametrization}\label{miscAccelerationsType:fromParametrization}
Reads a solution vector from file \configFile{inputfileSolution}{matrix}
which may be computed by a least squares adjustment (e.g. by \program{NormalsSolverVCE}).
The coefficients of the vector are interpreted from position \config{indexStart}
(counting from zero) with help of \configClass{parametrization}{parametrizationAccelerationType}.
If the solution file contains solution of several right hand sides you can choose
one with number \config{rightSide} (counting from zero).

The computed result is multiplied with \config{factor}.
)";
#endif

/***********************************************/

#include "files/fileMatrix.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief .
* @ingroup miscAccelerationsGroup
* @see MiscAccelerations */
class MiscAccelerationsFromParametrization : public MiscAccelerationsBase
{
  ParametrizationAccelerationPtr parametrization;
  Vector x;
  Double factor;

public:
  MiscAccelerationsFromParametrization(Config &config);

  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) override;
};

/***********************************************/

inline MiscAccelerationsFromParametrization::MiscAccelerationsFromParametrization(Config &config)
{
  try
  {
    FileName fileNamex;
    UInt     rightSide, indexStart;

    readConfig(config, "parametrization",   parametrization, Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileSolution", fileNamex,       Config::MUSTSET,  "",    "solution vector");
    readConfig(config, "indexStart",        indexStart,      Config::DEFAULT,  "0",   "position in the solution vector");
    readConfig(config, "rightSide",         rightSide,       Config::DEFAULT,  "0",   "if solution contains several right hand sides, select one");
    readConfig(config, "factor",            factor,          Config::DEFAULT,  "1.0", "the result is multiplied by this factor, set -1 to subtract the field");
    if(isCreateSchema(config)) return;

    Matrix mx;
    readFileMatrix(fileNamex, mx);
    x = mx.slice(indexStart, rightSide, parametrization->parameterCount(), 1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsFromParametrization::acceleration(SatelliteModelPtr satellite, const Time &time,
                                                                   const Vector3d &pos, const Vector3d &vel,
                                                                   const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides)
{
  try
  {
    Matrix A(3, parametrization->parameterCount());
    Matrix B(3, parametrization->parameterCountArc());
    parametrization->compute(satellite, time, pos, vel, rotSat, rotEarth, ephemerides, A, B);
    Vector a(3);
    if(A.columns()) matMult(factor, A, x.row(0, A.columns()), a);
    if(B.columns()) matMult(factor, B, x.row(A.columns(), B.columns()), a);
    return Vector3d(a);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
