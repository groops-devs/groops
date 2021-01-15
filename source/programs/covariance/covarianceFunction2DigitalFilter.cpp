/***********************************************/
/**
* @file covarianceFunction2DigitalFilter.cpp
*
* @brief Digital filter coefficients from covariance functions.
*
* @author Andreas Kvas
* @date 2018-03-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes digital filter coefficients for a \configClass{digital filter}{digitalFilterType} of given degree and
order. The filter coefficients are computed by fitting them to an approximated
impulse response represented by the cholesky factor of the covariance matrix.

The parameter \config{warmup} determines from which element of the cholesky matrix the
coefficients (default: half the covariance length) are fitted.

Per default, the program computes filter coefficients which generate colored noise
when applied to a white noise sequence. When \config{decorrelationFilter} is set,
a decorrelation filter is computed which yields white noise when applied to colored noise.
)";

/***********************************************/

#include "programs/program.h"
#include "base/fourier.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Digital filter coefficients from covariance functions.
* @ingroup programsGroup */
class CovarianceFunction2DigitalFilter
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(CovarianceFunction2DigitalFilter, SINGLEPROCESS, "Digital filter coefficients from covariance function.", Covariance)

/***********************************************/

void CovarianceFunction2DigitalFilter::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outFileFilter, inputFileCovariance;
    UInt column, p, q;
    Bool decorrelate;
    UInt start = MAX_UINT;

    readConfig(config, "outputfileFilter",     outFileFilter,       Config::MUSTSET, "", "filter coefficients");
    readConfig(config, "inputfileCovariance",  inputFileCovariance, Config::MUSTSET, "", "first column: time steps, following columns: covariance functions");
    readConfig(config, "column",               column,              Config::DEFAULT,  "1", "Column with covariance function to be fitted");
    readConfig(config, "warmup",               start,               Config::OPTIONAL, "", "number of samples until diagonal of Cholesky factor is flat (default: half covariance length)");
    readConfig(config, "numeratorDegree",      p,                   Config::DEFAULT,  "3", "Maximum degree of numerator polynomial (MA constituent)");
    readConfig(config, "denominatorDegree",    q,                   Config::DEFAULT,  "3", "Maximum degree of denominator polynomial (AR constitutent)");
    readConfig(config, "decorrelationFilter",  decorrelate,         Config::DEFAULT,  "0", "compute a decorrelation filter");
    if(isCreateSchema(config)) return;

    // read covariance function
    // ------------------------
    logStatus<<"read covariance function file <"<<inputFileCovariance<<">"<<Log::endl;
    Matrix covFunc;
    readFileMatrix(inputFileCovariance, covFunc);
    Vector c = covFunc.column(column);

    Matrix W(c.rows(), Matrix::SYMMETRIC, Matrix::UPPER);
    for(UInt i = 0; i<W.rows(); i++)
      for(UInt j = i; j<W.rows(); j++)
        W(i, j) = c(j-i);

    cholesky(W);

    if(decorrelate)
      inverse(W);

    // compute best fitting filter
    // ---------------------------
    if(start == MAX_UINT)
      start = W.rows()/2;

    Vector b, a;
    Matrix H_approx = W.trans().slice(start, start, W.rows()-start, W.rows()-start);

    Vector l = H_approx.slice(p, 0, H_approx.rows()-p, 1);
    Matrix A = H_approx.slice(p, 1, H_approx.rows()-p, q);

    // AR part
    Vector ar_hat = leastSquares(A, l);  // "alpha" AR coefficients
    a = Vector(ar_hat.rows()+1, 1.0);
    copy(-ar_hat, a.row(1, a.rows()-1)); // "beta" AR coefficients

    // MA part
    b = H_approx.slice(0, 0, p, q+1)*a;

    logStatus << "Writing filter to <" << outFileFilter << ">"<< Log::endl;
    Matrix F(std::max(b.rows(), a.rows()), 3);
    for(UInt i=0; i<F.rows(); i++)
      F(i, 0) = static_cast<Double>(i);
    copy(b, F.slice(0, 1, b.rows(), 1));
    copy(a, F.slice(0, 2, a.rows(), 1));
    writeFileMatrix(outFileFilter, F);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
