/***********************************************/
/**
* @file instrumentEstimateEmpiricalCovariance.cpp
*
* @brief Estimate the empirical covariance matrix of an instrument file.
**
* @author Andreas Kvas
* @date 2017-06-06
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates the empirical auto- and cross-covariance of selected data columns per arc
of \configFile{inputfileInstrument}{instrument}.
The maximum computed lag is determined by the number of \configFile{outputfileCovarianceMatrix}{matrix} specified
(for a single output file only the auto-covariance is determined, for two output files auto- and cross-covariance is computed and so on).

Stationarity is assumed for the input time series, which means the temporal covariance matrix has Toeplitz structure.
\begin{equation}
\begin{bmatrix}
\Sigma & \Sigma_{\Delta_1} & \Sigma_{\Delta_2} & \Sigma_{\Delta_3} & \Sigma_{\Delta_4} \\
       & \Sigma            & \Sigma_{\Delta_1} & \Sigma_{\Delta_2} & \Sigma_{\Delta_3} \\
       &                   & \Sigma            & \Sigma_{\Delta_1} & \Sigma_{\Delta_2} \\
       &                   &                   & \Sigma            & \Sigma_{\Delta_1} \\
       &                   &                   &                   & \Sigma            \\
\end{bmatrix}
\end{equation}

The matrix for lag $h$ describes the covariance between $x_{t-h}$ and $x_{t}$, i.e. $\Sigma(t-h, t)$.

To get a reliable estimate, \program{InstrumentDetrend} should be called first.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Estimate the empirical covariance matrix of an instrument file.
* @ingroup programsGroup */
class InstrumentEstimateEmpiricalCovariance
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentEstimateEmpiricalCovariance, PARALLEL, "Estimate the empirical covariance matrix of an instrument file", Instrument, Covariance)

/***********************************************/

void InstrumentEstimateEmpiricalCovariance::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    std::vector<FileName> fileNameOut;
    FileName              fileNameIn;
    UInt                  startData, countData = MAX_UINT;

    readConfig(config, "outputfileCovarianceMatrix", fileNameOut, Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileInstrument",        fileNameIn,  Config::MUSTSET,  "",  "");
    readConfig(config, "startDataFields",            startData,   Config::DEFAULT,  "0", "start");
    readConfig(config, "countDataFields",            countData,   Config::OPTIONAL, "",  "number of data fields (default: all after start)");
    if(isCreateSchema(config)) return;

    InstrumentFile instrumentFile(fileNameIn);
    const UInt dim = std::min(countData, instrumentFile.dataCount(TRUE/*mustDefined*/)-startData);

    Double countEpochs = 0;
    std::vector<Matrix> covarianceMatrix(fileNameOut.size(), Matrix(dim, dim));
    covarianceMatrix.front().setType(Matrix::SYMMETRIC);
    Parallel::forEach(instrumentFile.arcCount(), [&](UInt arcNo)
    {
      Matrix data = instrumentFile.readArc(arcNo).matrix();
      if(data.rows() < fileNameOut.size()) return;
      countEpochs += data.rows();

      rankKUpdate(1, data.column(startData+1, dim), covarianceMatrix.at(0)); // x_{t-h}*x_{t}^T
      for(UInt h=1; h<fileNameOut.size(); h++)
        matMult(1., data.slice(0, startData+1, data.rows()-h, dim).trans(), data.slice(h, startData+1, data.rows()-h, dim), covarianceMatrix.at(h)); // x_{t-h}*x_{t}^T
    }, comm); // forEach

    Parallel::reduceSum(countEpochs, 0, comm);
    for(UInt h=0; h<covarianceMatrix.size(); h++)
      Parallel::reduceSum(covarianceMatrix.at(h), 0, comm);

    if(Parallel::isMaster(comm))
    {
      for(UInt h=0; h<fileNameOut.size(); h++)
      {
        logStatus<<"write covariance matrix to <"<<fileNameOut.at(h)<<">."<<Log::endl;
        writeFileMatrix(fileNameOut.at(h), 1./(countEpochs-instrumentFile.arcCount()*h) * covarianceMatrix.at(h));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
