/***********************************************/
/**
* @file basisSplines.h
*
* @brief Basis Splines.
*
* @author Beate Klinger
* @date 2015-05-30
*
*/
/***********************************************/

#ifndef __GROOPS_BASISSPLINES__
#define __GROOPS_BASISSPLINES__

#include "base/importStd.h"
#include "matrix.h"

/***** CLASS ***********************************/

/** @brief Basis Splines.
* @ingroup base */
namespace BasisSplines
{
  /** @brief Basis Splines.
  * @param t in [0,1).
  * @param degree spline degree.
  * @return vector(degree+1) with splines coefficients. */
  inline Vector compute(Double t, UInt degree);
}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

inline Vector BasisSplines::compute(Double t, UInt degree)
{
  Vector b(degree+1);
  b(0) = 1.;
  for(UInt n=1; n<=degree; n++)
  {
    b(n) = t/n * b(n-1);
    for(UInt i=n; i-->1;)
      b(i) = (t+n-i)/n * b(i-1) - (t-i-1)/n * b(i);
    b(0) *= (1-t)/n;
  }

  return b;
}

/***********************************************/

#endif /* __GROOPS__ */
