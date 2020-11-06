/***********************************************/
/**
* @file covarianceFunction2PowerSpectralDensity.cpp
*
* @brief Power Spectral Density (PSD) from covariance functions.
*
* @author Torsten Mayer-Guerr
* @date 2011-07-03
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
One sided Power Spectral Density (PSD) from a covariance function. The first column of \configFile{inputfileCovarianceFunction}{matrix}
should contain the time lag in seconds.
Multiple covariance functions (in the following column)s are supported.
The output is a \file{matrix}{matrix} with first column contains the frequency $[Hz]$ and the other columns the PSD $[unit^2/Hz]$.

Conversion between covariance function $c_j$ and PSD $p_k$ is performed by discrete cosine transformation:
\begin{equation}
p_k = 2\Delta t\left(c_0 + c_{n-1} (-1)^k + \sum_{j=1}^{n-2} 2 c_j \cos(\pi jk/(n-1))\right).
\end{equation}

See also \program{PowerSpectralDensity2CovarianceFunction}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/fourier.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Power Spectral Density (PSD) from covariance functions.
* @ingroup programsGroup */
class CovarianceFunction2PowerSpectralDensity
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(CovarianceFunction2PowerSpectralDensity, SINGLEPROCESS, "Power Spectral Density (PSD) from covariance functions.", Covariance)
GROOPS_RENAMED_PROGRAM(CovarianceFunction2PSD, CovarianceFunction2PowerSpectralDensity, date2time(2020, 5, 24))

/***********************************************/

void CovarianceFunction2PowerSpectralDensity::run(Config &config)
{
  try
  {
    FileName psdName, covName;

    readConfig(config, "outputfilePSD",                psdName,   Config::MUSTSET,  "",  "first column: frequency [Hz], other columns PSD [unit^2/Hz]");
    readConfig(config, "inputfileCovarianceFunction",  covName,   Config::MUSTSET,  "",  "first column: time steps, following columns: covariance functions");
    if(isCreateSchema(config)) return;

    logStatus<<"read covariance function file <"<<covName<<">"<<Log::endl;
    Matrix covFunc;
    readFileMatrix(covName, covFunc);
    const Double sampling = covFunc(1,0) - covFunc(0,0);

    Matrix A(covFunc.rows(), covFunc.columns());
    for(UInt i=0; i<covFunc.rows(); i++)
      A(i,0) = i/(2*sampling*(A.rows()-1));
    for(UInt k=1; k<covFunc.columns(); k++)
      copy(Fourier::covariance2psd(covFunc.column(k), sampling), A.column(k));

    logStatus<<"write PSD file <"<<psdName<<">"<<Log::endl;
    writeFileMatrix(psdName, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
