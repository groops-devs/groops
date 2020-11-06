/***********************************************/
/**
* @file legendrePolynomial.h
*
* @brief Legendre polynomials.
* fully normalized.
*
* @author Torsten Mayer-Guerr
* @author Annette Eicker
* @date 2001-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_LEGEDNREPOLYNOMIAL__
#define __GROOPS_LEGEDNREPOLYNOMIAL__

#include "base/importStd.h"
#include "base/matrix.h"

/***** CLASS ***********************************/

/** @brief Legendre polynomials.
* @ingroup base
* fully normalized.
*/
class LegendrePolynomial
{
  static Vector factor1,            factor2;
  static Vector factor1Derivate,    factor2Derivate;
  static Vector factor1Derivate2nd, factor2Derivate2nd;
  static Vector factor1Integral,    factor2Integral;

  static void computeFactors(UInt degree);
  static void computeFactorsDerivate(UInt degree);
  static void computeFactorsDerivate2nd(UInt degree);
  static void computeFactorsIntegral(UInt degree);

public:
  /** @brief  Legendre polynomials.
  * (fully normalized).
  * @f[ P_n(t) \textnormal{ with } t=\cos(\psi),\quad n=0\ldots\textnormal{degree} @f]
  * @return vector with size @a degree+1. */
  static const Vector compute(Double t, UInt degree);

  /** @brief Derivative of Legendre polynomials.
  * (fully normalized).
  * @f[ P'_n(t) \textnormal{ with } t=\cos(\psi),\quad n=0\ldots\textnormal{degree} @f]
  * @return vector with size @a degree+1. */
  static const Vector derivative(Double t, UInt degree);

  /** @brief 2nd derivative of Legendre polynomials.
  * (fully normalized).
  * @f[ P''_n(t) \textnormal{ with } t=\cos(\psi),\quad n=0\ldots\textnormal{degree} @f]
  * @return vector with size @a degree+1. */
  static const Vector derivative2nd(Double t, UInt degree);

  /** @brief Integral over Legendre polynomials.
  * For a spherical cap (from 1 to t=cos(psi)).
  * (fully normalized).
  * @f[ \int_t^1 P_n(t)\,dt \textnormal{ with } t=\cos(\psi),\quad n=0\ldots\textnormal{degree} @f]
  * @return vector with size @a degree+1. */
  static const Vector integral(Double t, UInt degree);

  /** @brief Sum of Legendre polynomials.
  * (fully normalized).
  * Clenshaw Algorithm.
  * @f[ \sum_{n=0}^N k_n P_n(t) \textnormal{ with } t=\cos(\psi) @f]
  * @param t aperture angle \f$\cos(\psi)\f$.
  * @param coeff coefficients \f$k_n\f$
  * @param degree degree \f$N\f$ */
  static Double sum(Double t, const Vector &coeff, UInt degree);

  /** @brief Sum of the derivative of Legendre polynomials.
  * (fully normalized).
  * Clenshaw Algorithm.
  * @f[ \sum_{n=0}^N k_n P'_n(t) \textnormal{ with } t=\cos(\psi) @f]
  * @param t aperture angle \f$\cos(\psi)\f$.
  * @param coeff coefficients \f$k_n\f$
  * @param degree degree \f$N\f$
  */
  static Double sumDerivative(Double t, const Vector &coeff, UInt degree);

  /** @brief Sum of the 2nd derivative of Legendre polynomials.
  * (fully normalized).
  * Clenshaw Algorithm.
  * @f[ \sum_{n=0}^N k_n P''_n(t) \textnormal{ with } t=\cos(\psi) @f]
  * @param t aperture angle \f$\cos(\psi)\f$.
  * @param coeff coefficients \f$k_n\f$
  * @param degree degree \f$N\f$ */
  static Double sumDerivative2nd(Double t, const Vector &coeff, UInt degree);

  /** @brief Zero crossing sof Legendre polynomials.
  * The weights for Gauss-Legendre Integration are computed additionally.
  * @param degree degree \f$N\f$
  * @param[out] zeros zero crossings in a vector of size @a degree
  * @param[out] weights weights in a vector of size @a degree */
  static void zeros(UInt degree, Vector &zeros, Vector &weights);
}; // end LegendrePolynomial


/***********************************************/

#endif /* __GROOPS_LEGENDREPOLYNOMIAL__ */



