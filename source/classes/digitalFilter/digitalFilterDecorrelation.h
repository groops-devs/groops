/***********************************************/
/**
* @file DigitalFilterDecorrelation.h
*
* @brief Decorrelation filter.
*
* @author Andreas Kvas
* @date 2016-06-21
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERDECORRELATION__
#define __GROOPS_DIGITALFILTERDECORRELATION__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterDecorrelation = R"(
\subsection{Decorrelation}
Moving average decorrelation filter based on eigendecomposition of a Toeplitz covariance matrix.
)";
#endif

/***********************************************/

#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Decorrelation filter.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterDecorrelation : public DigitalFilterARMA
{
public:
  DigitalFilterDecorrelation(Config &config);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterDecorrelation::DigitalFilterDecorrelation(Config &config)
{
  try
  {
    FileName fileCovarianceFunction;

    readConfig(config, "inputfileCovarianceFunction",            fileCovarianceFunction,            Config::MUSTSET, "", "covariance function of time series");
    readConfig(config, "inFrequencyDomain", inFrequencyDomain, Config::DEFAULT,  "0", "apply filter in frequency domain");
    readConfig(config, "padType",           padType,           Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    Matrix cov;
    readFileMatrix(fileCovarianceFunction, cov);

    Matrix C(cov.rows(), Matrix::SYMMETRIC, Matrix::UPPER);
    for(UInt r = 0; r<cov.rows(); r++)
      for(UInt c = r; c<cov.rows(); c++)
        C(r, c) = cov(c-r,1);

    Vector e = eigenValueDecomposition(C);
    bn = Vector(cov.rows());
    bnStartIndex = bn.rows()/2;

    Vector qi = C.row(bnStartIndex).trans();
    for(UInt k=0; k<qi.rows(); k++)
      qi(k) *= 1./std::sqrt(e(k));

    matMult(1.0, C, qi, bn);
    for(UInt k = 0; k<bnStartIndex; k++)
      std::swap(bn(k), bn(bn.rows()-1-k));

    an = Vector(1, 1.);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
