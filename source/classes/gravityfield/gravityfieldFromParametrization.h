/***********************************************/
/**
* @file gravityfieldFromParametrization.h
*
* @brief Gravity field from parameter vector.
* @see Gravityfield
* @see ParametrizationGravity
*
* @author Torsten Mayer-Guerr
* @date 2003-10-31
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDFROMPARAMETRIZATION__
#define __GROOPS_GRAVITYFIELDFROMPARAMETRIZATION__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldFromParametrization = R"(
\subsection{FromParametrization}\label{gravityfieldType:fromParametrization}
Reads a solution vector from file \configFile{inputfileSolution}{matrix}
which may be computed by a least squares adjustment (e.g. by \program{NormalsSolverVCE}).
The coefficients of the vector are interpreted from position \config{indexStart}
(counting from zero) with help of \configClass{parametrizationGravity}{parametrizationGravityType}.
If the solution file contains solution of several right hand sides you can choose
one with number \config{rightSide} (counting from zero).
You can also read a vector from file \configFile{inputfileSigmax}{matrix}
containing the accuracies of the coefficients.

The computed result is multiplied with \config{factor}.
)";
#endif

/***********************************************/

#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Gravity field from parameter vector.
* @ingroup gravityfieldGroup
* @see Gravityfield
* @see ParametrizationGravity  */
class GravityfieldFromParametrization : public GravityfieldBase
{
  ParametrizationGravityPtr parametrization;
  Vector  x;
  Double  factor;
  Matrix  C; // full covariance matrix
  Vector  sigma2x;

public:
  GravityfieldFromParametrization(Config &config);

  Double   potential      (const Time &time, const Vector3d &point) const override;
  Double   radialGradient (const Time &time, const Vector3d &point) const override;
  Double   field          (const Time &time, const Vector3d &point, const Kernel &kernel) const override;
  Vector3d gravity        (const Time &time, const Vector3d &point) const override;
  Tensor3d gravityGradient(const Time &time, const Vector3d &point) const override;
  Vector3d deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const override;
  void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                           const Vector &hn, const Vector &ln, std::vector< std::vector<Vector3d> > &disp) const override;
  SphericalHarmonics sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const override;
  void variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const override;
};

/***********************************************/

#endif /* __GROOPS_GRAVITYFIELD__ */
