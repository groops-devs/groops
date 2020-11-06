/***********************************************/
/**
* @file sphericalHarmonicsNumberingOrderNonAlternating.h
*
* @brief Numbering schema of spherical harmonics coefficients.
* Numbering order by order with cnm, snm non-alternating:
* c20, c30, c40, [...], c21, c31, c41, [...], s21, s31, s41, [...]
* @see SphericalHarmonicsNumbering
*
* @author Judith Schall
* @date 2009-08-24
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICSNUMBERINGORDERNONALTERNATING__
#define __GROOPS_SPHERICALHARMONICSNUMBERINGORDERNONALTERNATING__

#ifdef DOCSTRING_SphericalHarmonicsNumbering
static const char *docstringSphericalHarmonicsNumberingOrderNonAlternating = R"(
\subsection{OrderNonAlternating}
Numbering order by order with cnm, snm non-alternating:
\[ c20, c30, c40, \ldots, c21, c31, c41, \ldots, s21, s31, s41, \]
)";
#endif

/***********************************************/

#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Numbering schema of spherical harmonics coefficients.
* @ingroup sphericalHarmonicsNumberingGroup
* Numbering order by order with cnm, snm non-alternating:
* c20, c30, c40, [...], c21, c31, c41, [...], s21, s31, s41, [...]
* @see SphericalHarmonicsNumbering */
class SphericalHarmonicsNumberingOrderNonAlternating : public SphericalHarmonicsNumbering
{
public:
  SphericalHarmonicsNumberingOrderNonAlternating(Config &/*config*/) {}

  virtual UInt parameterCount(UInt maxDegree, UInt minDegree) const override {return (maxDegree+1)*(maxDegree+1) - minDegree*minDegree;}
  virtual void numbering(UInt maxDegree, UInt minDegree, std::vector<std::vector<UInt>> &Cnm, std::vector<std::vector<UInt>> &Snm) const override;
};

/***********************************************/

inline void SphericalHarmonicsNumberingOrderNonAlternating::numbering(UInt maxDegree, UInt minDegree, std::vector<std::vector<UInt>> &Cnm, std::vector<std::vector<UInt>> &Snm) const
{
  Cnm.resize(0); Cnm.resize(maxDegree+1);
  Snm.resize(0); Snm.resize(maxDegree+1);
  for(UInt n=0; n<=maxDegree; n++)
  {
    Cnm[n].resize(n+1, NULLINDEX);
    Snm[n].resize(n+1, NULLINDEX);
  }

  UInt idx = 0;
  for(UInt n=minDegree; n<=maxDegree; n++)
    Cnm[n][0] = idx++;
  for(UInt m=1; m<=maxDegree; m++)
  {
    for(UInt n=std::max(minDegree,m); n<=maxDegree; n++)
      Cnm[n][m] = idx++;
    for(UInt n=std::max(minDegree,m); n<=maxDegree; n++)
      Snm[n][m] = idx++;
  }
}

/***********************************************/

#endif
