/***********************************************/
/**
* @file sphericalHarmonicsNumberingOrder.h
*
* @brief Numbering schema of spherical harmonics coefficients.
* Numbering order by order:
* c20, c30, c40, ..., c21, s21, c31, s31, ..., c22, s22
* @see SphericalHarmonicsNumbering
*
* @author Torsten Mayer-Guerr
* @date 2009-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICSNUMBERINGORDER__
#define __GROOPS_SPHERICALHARMONICSNUMBERINGORDER__

#ifdef DOCSTRING_SphericalHarmonicsNumbering
static const char *docstringSphericalHarmonicsNumberingOrder = R"(
\subsection{Order}
Numbering order by order:
\[ c20, c30, c40, \ldots, c21, s21, c31, s31, \ldots, c22, s22 \]
)";
#endif

/***********************************************/

#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Numbering schema of spherical harmonics coefficients.
* @ingroup sphericalHarmonicsNumberingGroup
* Numbering order by order:
* c20, c30, c40, ..., c21, s21, c31, s31, ..., c22, s22
* @see SphericalHarmonicsNumbering */
class SphericalHarmonicsNumberingOrder : public SphericalHarmonicsNumbering
{
public:
  SphericalHarmonicsNumberingOrder(Config &/*config*/) {}

  virtual UInt parameterCount(UInt maxDegree, UInt minDegree) const override {return (maxDegree+1)*(maxDegree+1) - minDegree*minDegree;}
  virtual void numbering(UInt maxDegree, UInt minDegree, std::vector<std::vector<UInt>> &Cnm, std::vector<std::vector<UInt>> &Snm) const override;
};

/***********************************************/

inline void SphericalHarmonicsNumberingOrder::numbering(UInt maxDegree, UInt minDegree, std::vector<std::vector<UInt>> &Cnm, std::vector<std::vector<UInt>> &Snm) const
{
  Cnm.clear(); Cnm.resize(maxDegree+1);
  Snm.clear(); Snm.resize(maxDegree+1);
  for(UInt n=0; n<=maxDegree; n++)
  {
    Cnm[n].resize(n+1, NULLINDEX);
    Snm[n].resize(n+1, NULLINDEX);
  }

  UInt idx = 0;
  for(UInt n=minDegree; n<=maxDegree; n++)
    Cnm[n][0] = idx++;
  for(UInt m=1; m<=maxDegree; m++)
    for(UInt n=std::max(minDegree,m); n<=maxDegree; n++)
    {
      Cnm[n][m] = idx++;
      Snm[n][m] = idx++;
    }
}

/***********************************************/

#endif
