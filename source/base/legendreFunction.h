/***********************************************/
/**
* @file legendreFunction.h
*
* @brief Associated Legendre functions.
* (fully normalized).
*
* @author Torsten Mayer-Guerr
* @author Annette Eicker
* @date 2001-05-31
*
*/
/***********************************************/

#ifndef __LEGEDNREFUNCTION__
#define __LEGEDNREFUNCTION__

#include "base/importStd.h"
#include "base/matrix.h"

/***** CLASS ***********************************/

/** @brief Associated Legendre functions.
* @ingroup base
* fully normalized. */
class LegendreFunction
{
  static Matrix factor1, factor2;
  static Matrix factor1Integral, factor2Integral;
  static Vector factorSmall;

  static void computeFactors(UInt degree);
  static void computeFactorsIntegral(UInt degree);

public:
  /** @brief Legendre functions.
  * (fully normalized).
  * @f[ P_n^m(t) \textnormal{ with } t=\cos(\psi),\quad n,m=0\ldots\textnormal{degree} @f]
  * @return lower triangular matrix with dimension @a degree+1.
  */
  static const Matrix compute(Double t, UInt degree);

  /** @brief integral over Legendre functions.
  * (fully normalized).
  * @f[ I_n^m = \int_{t_1}^{t_2} P_n^m(t)\,dt \textnormal{ with } t=\cos(\psi),\quad n,m=0\ldots\textnormal{degree} @f]
  * @return lower triangular matrix with dimension @a degree+1.
  */
  static const Matrix integral(Double t1, Double t2, UInt degree);
}; // end LegendreFunction

/***********************************************/

#endif /* __LEGENDREFUNCTION__ */



