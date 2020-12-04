/***********************************************/
/**
* @file sphericalHarmonicsFilterMatrix.h
*
* @brief Smoothing by a filter matrix.
* @see SphericalHarmonicsFilter
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2017-02-21
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICSFILTERMATRIX__
#define __GROOPS_SPHERICALHARMONICSFILTERMATRIX__

// Latex documentation
#ifdef DOCSTRING_SphericalHarmonicsFilter
static const char *docstringSphericalHarmonicsFilterMatrix = R"(
\subsection{Matrix}
Filtering the spherical harmonics expansion with a matrix filter.
)";
#endif

/***********************************************/

#include "files/fileMatrix.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"
#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilter.h"

/***** CLASS ***********************************/

/** @brief Smoothing by a filter matrix.
* @ingroup sphericalHarmonicsFilterGroup
* @see SphericalHarmonicsFilter */
class SphericalHarmonicsFilterMatrix : public SphericalHarmonicsFilterBase
{
  Matrix A;
  UInt   minDegree, maxDegree;
  std::vector< std::vector<UInt> > idxC, idxS;

public:
  SphericalHarmonicsFilterMatrix(Config &config);

  SphericalHarmonics filter(const SphericalHarmonics &harm) const;
};

/***********************************************/

inline SphericalHarmonicsFilterMatrix::SphericalHarmonicsFilterMatrix(Config &config)
{
  try
  {
    FileName fileNameMatrix;
    SphericalHarmonicsNumberingPtr numbering;

    readConfig(config, "inputfileMatrix", fileNameMatrix, Config::MUSTSET, "",  "");
    readConfig(config, "minDegree",       minDegree,      Config::MUSTSET, "2", "of matrix");
    readConfig(config, "maxDegree",       maxDegree,      Config::MUSTSET, "",  "of matrix");
    readConfig(config, "numbering",       numbering,      Config::MUSTSET, "",  "numbering scheme of the matrix");
    if(isCreateSchema(config)) return;

    readFileMatrix(fileNameMatrix, A);

    numbering->numbering(maxDegree, minDegree, idxC, idxS);
    const UInt dim = numbering->parameterCount(maxDegree, minDegree);
    if((A.rows() != dim) || (A.columns() != dim))
      throw(Exception("Matrix dimension and numbering disagree"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline SphericalHarmonics SphericalHarmonicsFilterMatrix::filter(const SphericalHarmonics &harm) const
{
  try
  {
    // sort into matrix ordering
    const Vector xDegreeWise = harm.get(maxDegree).x(); // degree wise
    Vector x(A.columns());
    UInt idx = 0;
    for(UInt n=0; n<=maxDegree; n++)
    {
      if(idxC[n][0]!=NULLINDEX) x(idxC[n][0]) = xDegreeWise(idx);
      idx++;
      for(UInt m=1; m<=n; m++)
      {
        if(idxC[n][m]!=NULLINDEX) x(idxC[n][m]) = xDegreeWise(idx);
        if(idxS[n][m]!=NULLINDEX) x(idxS[n][m]) = xDegreeWise(idx + 1);
        idx += 2;
      }
    }

    // filter
    x = A*x;

    // sort back
    Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    for(UInt n=minDegree; n<=maxDegree; n++)
    {
      if(idxC[n][0]!=NULLINDEX) cnm(n,0) = x(idxC[n][0]);
      for(UInt m=1; m<=n; m++)
      {
        if(idxC[n][m]!=NULLINDEX) cnm(n,m) = x(idxC[n][m]);
        if(idxS[n][m]!=NULLINDEX) snm(n,m) = x(idxS[n][m]);
      }
    }

    return SphericalHarmonics(harm.GM(), harm.R(), cnm, snm).get(harm.maxDegree());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
