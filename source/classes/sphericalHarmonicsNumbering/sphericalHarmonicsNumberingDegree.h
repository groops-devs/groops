/***********************************************/
/**
* @file sphericalHarmonicsNumberingDegree.h
*
* @brief Numbering schema of spherical harmonics coefficients.
* Numbering degree by degree:
* c20, c21, s21, c22, s22, c30, c31, s31, c32, s32, ...
* @see SphericalHarmonicsNumbering
*
* @author Torsten Mayer-Guerr
* @date 2009-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICSNUMBERINGDEGREE__
#define __GROOPS_SPHERICALHARMONICSNUMBERINGDEGREE__

// Latex documentation
#ifdef DOCSTRING_SphericalHarmonicsNumbering
static const char *docstringSphericalHarmonicsNumberingDegree = R"(
\subsection{Degree}
Numbering degree by degree:
\[ c20, c21, s21, c22, s22, c30, c31, s31, c32, s32,\ldots \]
)";
#endif

/***********************************************/

#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Numbering schema of spherical harmonics coefficients.
* @ingroup sphericalHarmonicsNumberingGroup
* Numbering degree by degree:
* c20, c21, s21, c22, s22, c30, c31, s31, c32, s32, ...
* @see SphericalHarmonicsNumbering */
class SphericalHarmonicsNumberingDegree : public SphericalHarmonicsNumbering
{
public:
  SphericalHarmonicsNumberingDegree(Config &/*config*/) {}

  virtual UInt parameterCount(UInt maxDegree, UInt minDegree) const override {return (maxDegree+1)*(maxDegree+1) - minDegree*minDegree;}
  virtual void numbering(UInt maxDegree, UInt minDegree, std::vector<std::vector<UInt>> &Cnm, std::vector<std::vector<UInt>> &Snm) const override;
};

/***********************************************/

inline void SphericalHarmonicsNumberingDegree::numbering(UInt maxDegree, UInt minDegree, std::vector<std::vector<UInt>> &Cnm, std::vector<std::vector<UInt>> &Snm) const
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
  {
    Cnm[n][0] = idx++;
    for(UInt m=1; m<=n; m++)
    {
      Cnm[n][m] = idx++;
      Snm[n][m] = idx++;
    }
  }
}

/***********************************************/

#endif
