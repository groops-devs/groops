/***********************************************/
/**
* @file powerSpectralDensity2CovarianceFunction.cpp
*
* @brief Covariance functions from Power Spectral Density (PSD).
*
* @author Matthias Ellmer
* @date 2013-08-27
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Covariance function from Power Spectral Density (PSD).
The \configFile{inputfilePSD}{matrix} contains in the first column the frequency $[Hz]$, followed by (possibly multiple) PSDs $[unit^2/Hz]$.
The output is a \file{matrix}{matrix}, the first column containing time lag $[s]$ and the other columns the covariance functions $[unit^2]$.
Conversion between PSD $p_j$ and covariance function $c_k$ is performed by discrete cosine transformation:
\begin{equation}
c_k = \frac{1}{4\Delta t (n-1)}\left(p_0 + p_{n-1} (-1)^k + \sum_{j=1}^{n-2} 2 p_j \cos(\pi jk/(n-1))\right).
\end{equation}

See also \program{CovarianceFunction2PowerSpectralDensity}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/fourier.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Covariance functions from Power Spectral Density (PSD).
* @ingroup programsGroup */
class PowerSpectralDensity2CovarianceFunction
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(PowerSpectralDensity2CovarianceFunction, SINGLEPROCESS, "Covariance functions from Power Spectral Density (PSD).", Covariance)
GROOPS_RENAMED_PROGRAM(CovariancePsd2CovarianceFunction, PowerSpectralDensity2CovarianceFunction, date2time(2020, 5, 24))

/***********************************************/

void PowerSpectralDensity2CovarianceFunction::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameCov, fileNamePsd;

    readConfig(config, "outputfileCovarianceFunction", fileNameCov, Config::MUSTSET, "", "first column: time steps [seconds], following columns: covariance functions");
    readConfig(config, "inputfilePSD",                 fileNamePsd, Config::MUSTSET, "", "first column: frequency [Hz], following columns PSD [unit^2/Hz]");
    if(isCreateSchema(config)) return;

    logStatus<<"read PSD file <"<<fileNamePsd<<">"<<Log::endl;
    Matrix psd;
    readFileMatrix(fileNamePsd, psd);

    Matrix cov(psd.rows(), psd.columns());
    const Double dt = 1./(2*psd(psd.rows()-1, 0));
    for(UInt i=0; i<cov.rows(); i++)
      cov(i,0) = i*dt;

    for(UInt k=1; k<cov.columns(); k++)
      copy(Fourier::psd2covariance(psd.column(k), dt), cov.column(k));

    logStatus<<"write covariance file <"<<fileNameCov<<">"<<Log::endl;
    writeFileMatrix(fileNameCov, cov);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
