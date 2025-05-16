/***********************************************/
/**
* @file sphericalHarmonics.h
*
* @brief Spherical harmonics and functions of the gravity field.
*
* (4Pi normalized).
*
* @author Torsten Mayer-Guerr
* @date 2001-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICS__
#define __GROOPS_SPHERICALHARMONICS__

#include "base/importStd.h"
#include "base/matrix.h"
#include "base/vector3d.h"

/***** CLASS ***********************************/

class Rotary3d;
class Tensor3d;

/***** CLASS ***********************************/

/** @brief Spherical harmonics and functions of the gravity field.
* @ingroup base
* (4Pi normalized).
* Implements LinearSpace: SphericalHarmonics can be added/subtracted and be scaled with a Double. */
class SphericalHarmonics
{
  Double _GM, _R;
  UInt   _maxDegree;
  Matrix _cnm, _snm, _sigma2cnm, _sigma2snm;
  Bool   _interior;

  static Matrix factor1, factor2;

  // factors needed for recursion formular
  static void computeFactors(UInt degree);

public:
  /// Default Constructor.
  explicit SphericalHarmonics(Bool interior=FALSE);

  /// Constructor from potential coefficients.
  SphericalHarmonics(Double GM, Double R, const const_MatrixSlice &cnm, const const_MatrixSlice &snm, Bool interior=FALSE);

  /// Constructor from potential coefficients and variances.
  SphericalHarmonics(Double GM, Double R, const const_MatrixSlice &cnm, const const_MatrixSlice &snm, const const_MatrixSlice &sigma2cnm, const const_MatrixSlice &sigma2snm, Bool interior=FALSE);

  /** @brief Solid spherical harmonics (Cnm und Snm).
  * Basis functions (4Pi normalized).
  * @f[ C_{nm}(\lambda,\vartheta,r) = \frac{1}{r^{n+1}} \cos(m\lambda)P_n^m(\cos\vartheta) @f]
  * @f[ S_{nm}(\lambda,\vartheta,r) = \frac{1}{r^{n+1}} \sin(m\lambda)P_n^m(\cos\vartheta) @f]
  * if interior:
  * @f[ C_{nm}(\lambda,\vartheta,r) = r^n \cos(m\lambda)P_n^m(\cos\vartheta) @f]
  * @f[ S_{nm}(\lambda,\vartheta,r) = r^n \sin(m\lambda)P_n^m(\cos\vartheta) @f] */
  static void CnmSnm(const Vector3d &point, UInt maxDegree, Matrix &Cnm, Matrix &Snm, Bool interior=FALSE);

  /** @brief Solid Legendre functions (Pnm).
  * (4Pi normalized).
  * @f[ P_{nm}(\lambda,\vartheta,r) = \frac{1}{r^{n+1}} P_n^m(\cos\vartheta) @f]
  * if interior
  * @f[ P_{nm}(\lambda,\vartheta,r) = r^n P_n^m(\cos\vartheta) @f] */
  static Matrix Pnm(Angle theta, Double r, UInt degree, Bool interior=FALSE);

  /** @brief Interpret coefficients as inner/outer space harmonics. */
  void setInterior(Bool interior=TRUE) {_interior = interior;}

  /** @brief Interpret coefficients as inner space harmonics? */
  Bool isInterior() const {return _interior;}

  /** @brief Convert spherical harmonics.
  * Based on given GM and R.
  * @param maxDegree maximum degree.
  * @param minDegree minimum degree.
  * @param GM Specific gravitational constant, 0: the internal value is used.
  * @param R  Reference radius, 0: the internal value is used. */
  SphericalHarmonics get(UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const;

  /** @brief Rotate spherical harmonics.
  * Variances are not propagated. */
  SphericalHarmonics rotate(const Rotary3d &rotary) const;

  /** @brief Laplace spherical harmonics.
  * @f[ Y_n(P) = \frac{GM}{R} \sum_{m=0}^n  c_{nm}C^r_{nm}(\frac{1}{R}P) + s_{nm}S^r_{nm}(\frac{1}{R}P) @f] */
  Vector Yn(const Vector3d &point, UInt maxDegree=INFINITYDEGREE) const;

  /** @brief Potential.
  * @f[ V(P) = \sum_{n=minDegree}^{maxDegree} Y_n(P) @f]
  * @see SphericalHarmonics::Yn */
  Double potential(const Vector3d &point, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0) const;

  /** @brief Radial derivative of the potential.
  * @f[ \frac{\partial V(P)}{\partial r} @f] */
  Double radialGradient(const Vector3d &point, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0) const;

  /** @brief Gravitational acceleration.
  * @f[ \mathbf{g}(P) = \nabla V(P) @f] */
  Vector3d gravity(const Vector3d &point, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0) const;

  /** @brief Gravitational gradient (Tensor).
  * @f[ T(P) = \nabla \nabla V(P) @f] */
  Tensor3d gravityGradient(const Vector3d &point, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0) const;

  /** @brief Deformation due to loading.
  * @param  point station position [m]
  * @param  gravity local gravity at station [m/s**2]
  * @param  hn vertical load love numbers
  * @param  ln horizontal load love numbers
  * @param  maxDegree maximum degree (inclusive)
  * @param  minDegree minimum degree (inclusive)
  * @return deformation in TRF [m] */
  Vector3d deformation(const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0) const;

  /// Reference radius.
  Double R() const {return _R;}

  /// Gravitational constant * mass.
  Double GM() const {return _GM;}

  /// Maximum degree.
  UInt maxDegree() const {return _maxDegree;}

        Matrix &cnm()       {return _cnm;} //!< Potential coefficients.
        Matrix &snm()       {return _snm;} //!< Potential coefficients.
  const Matrix &cnm() const {return _cnm;} //!< Potential coefficients.
  const Matrix &snm() const {return _snm;} //!< Potential coefficients.

        Matrix &sigma2cnm()       {return _sigma2cnm;} //!< Variances of the potential coefficients.
        Matrix &sigma2snm()       {return _sigma2snm;} //!< Variances of the potential coefficients.
  const Matrix &sigma2cnm() const {return _sigma2cnm;} //!< Variances of the potential coefficients.
  const Matrix &sigma2snm() const {return _sigma2snm;} //!< Variances of the potential coefficients.


  /** @brief potential coefficients given as vector.
  * the potential coefficients are given in a
  * degree wise sequence with alternating sin, cos:
  * ..., c20,c21,s21,c22,s22,...
  * beginning with c00. */
  Vector x() const;

  /** @brief variances of potential coefficients given as vector.
  * the variances of the potential coefficients are given in a
  * degree wise sequence with alternating sin, cos:
  * ..., c20,c21,s21,c22,s22,...
  * beginning with c00. */
  Vector sigma2x() const;

  // linearSpace Operations
  SphericalHarmonics &operator+= (const SphericalHarmonics &harm);
  SphericalHarmonics &operator-= (const SphericalHarmonics &harm);
  SphericalHarmonics &operator*= (Double c);
}; // end SphericalHarmonics

/***********************************************/

inline SphericalHarmonics operator- (const SphericalHarmonics &t)                  {return SphericalHarmonics(t)  *= -1;}
inline SphericalHarmonics operator+ (const SphericalHarmonics &t1, const SphericalHarmonics &t2) {return SphericalHarmonics(t1) += t2;}
inline SphericalHarmonics operator- (const SphericalHarmonics &t1, const SphericalHarmonics &t2) {return SphericalHarmonics(t1) -= t2;}
inline SphericalHarmonics operator* (Double c, const SphericalHarmonics &t)        {return SphericalHarmonics(t)  *=c;}
inline SphericalHarmonics operator* (const SphericalHarmonics &t, Double c)        {return SphericalHarmonics(t)  *=c;}

/***********************************************/

#endif /* __GROOPS_SPHERICALHARMONICS__ */



