/***********************************************/
/**
* @file gravityfieldTimeSplines.h
*
* @brief TimeSplines.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2004-04-14
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDTIMESPLINES__
#define __GROOPS_GRAVITYFIELDTIMESPLINES__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldTimeSplines = R"(
\subsection{TimeSplines}\label{gravityfieldType:timeSplines}
Read a time variable gravity field from file
\configFile{inputfileTimeSplinesGravityfield}{timeSplinesGravityField}
represented by a spherical harmonics expansion in the spatial domain and spline functions
in the time domain. If set the expansion is limited in the range between
\config{minDegree} and \config{maxDegree} inclusivly.

This file can be created for example by \program{Gravityfield2TimeSplines} or
\program{PotentialCoefficients2BlockMeanTimeSplines}.

The computed result is multiplied with \config{factor}.
)";
#endif

/***********************************************/

#include "files/fileTimeSplinesGravityfield.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief TimeSplines.
* @ingroup gravityfieldGroup
* @see Gravityfield */
class GravityfieldTimeSplines : public GravityfieldBase
{
  Bool    hasCovariance;
  Double  factor;
  UInt    minDegree;
  UInt    maxDegree;
  mutable InFileTimeSplinesGravityfield splinesFile;
  mutable InFileTimeSplinesCovariance   covarianceFile;

public:
  GravityfieldTimeSplines(Config &config);

  // Einfluss des Referenzfieldes im Aufpunkt:
  Double   potential      (const Time &time, const Vector3d &point) const;
  Double   radialGradient (const Time &time, const Vector3d &point) const;
  Double   field          (const Time &time, const Vector3d &point, const Kernel &kernel) const;
  Vector3d gravity        (const Time &time, const Vector3d &point) const;
  Tensor3d gravityGradient(const Time &time, const Vector3d &point) const;
  Vector3d deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const;
  void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                           const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const;

  // Umrechnung der Parameter in Potentialkoeffizienten
  SphericalHarmonics sphericalHarmonics(const Time &time, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const;
  Matrix sphericalHarmonicsCovariance  (const Time &time, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const;

  void variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const;
};

/***********************************************/

#endif /* __GROOPS_GRAVITYFIELD__ */
