/***********************************************/
/**
* @file instrument2CovarianceFunctionVCE.cpp
*
* @brief Covariance functions via variance component estimation.
*
* @author Torsten Mayer-Guerr
* @date 2017-07-09
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This estimates a covariance function of \configFile{inputfileInstrument}{instrument}
for all selected columns with \config{startDataFields} and \config{countDataFields}.
The estimation is performed robustly via variance component estimation.
Bad arcs are downweigthed and the accuracies can be written with \configFile{outputfileSigmasPerArc}{matrix}.
The length of the covariance functions are determined by the longest arc.
Additionaly the data can be detrended with \configClass{parameter}{parametrizationTemporalType}
and \configClass{parameterPerArc}{parametrizationTemporalType}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "misc/varianceComponentEstimation.h"

/***** CLASS ***********************************/

/** @brief Covariance functions via variance component estimation.
* @ingroup programsGroup */
class Instrument2CovarianceFunctionVCE
{
  Epoch::Type                    arcType;
  std::vector<std::vector<Time>> arcTimes;
  std::vector<Matrix>            arcData;
  UInt                           startData, countData;
  ParametrizationTemporalPtr     parameter;
  ParametrizationTemporalPtr     parameterPerArc;

  Matrix N;        // normal equation matrix
  Vector n;        // right hand sides
  Vector x;        // solution
  Matrix Wz;       // monte carlo vector for redundancy computation
  Double lPl;      // =l'Pl, weighted norm of the observations
  UInt   obsCount; // number of observations

  std::vector<Matrix> W;                // cholesky of the covariance matrix (per data column)
  Vector              sigma, sigmaNew;  // per arc
  Matrix              CosTransform;
  Matrix              covFunc;
  Matrix              Psd;
  Matrix              ePe, redundancy;  // one row for each frequency, one column for each component

  void computeObservationEquation(UInt arcNo, std::vector<Matrix> &W, Vector &Wl, Matrix &WA, Matrix &WB);
  void buildNormals(UInt arcNo);
  void computeRedundancies(UInt arcNo);
  Arc  computeResiduals(UInt arcNo);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Instrument2CovarianceFunctionVCE, PARALLEL, "Covariance functions via variance component estimation", Instrument, Covariance)

/***********************************************/

void Instrument2CovarianceFunctionVCE::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameCov, fileNameArcSigma, fileNameResiduals, fileNameSolution;
    FileName fileNameIn;
    UInt     iterCount;
    countData = MAX_UINT;

    readConfig(config, "outputfileCovarianceFunction", fileNameCov,       Config::MUSTSET,   "", "covariance functions");
    readConfig(config, "outputfileSigmasPerArc",       fileNameArcSigma,  Config::OPTIONAL,  "", "accuracies of each arc");
    readConfig(config, "outputfileResiduals",          fileNameResiduals, Config::OPTIONAL,  "", "");
    readConfig(config, "outputfileSolution",           fileNameSolution,  Config::OPTIONAL,  "", "estimated parameter vector (global part only)");
    readConfig(config, "inputfileInstrument",          fileNameIn,        Config::MUSTSET,   "",  "");
    readConfig(config, "startDataFields",              startData,         Config::DEFAULT,   "0", "start");
    readConfig(config, "countDataFields",              countData,         Config::OPTIONAL,  "",  "number of data fields (default: all after start)");
    readConfig(config, "parameter",                    parameter,         Config::DEFAULT,   "",  "data is reduced by temporal representation");
    readConfig(config, "parameterPerArc",              parameterPerArc,   Config::DEFAULT,   "",  "data is reduced by temporal representation");
    readConfig(config, "iterationCount",               iterCount,         Config::DEFAULT,   "5",  "number of iterations for the estimation");
    if(isCreateSchema(config)) return;

    // =============================================

    // Determine max. length of ovariance functions
    // --------------------------------------------
    InstrumentFile instrumentFile(fileNameIn);
    arcType = instrumentFile.getType();
    const UInt arcCount = instrumentFile.arcCount();
    if(arcCount < 5)
      throw(Exception("Need at least 5 arcs for a reliable covariance estimation"));

    UInt   covLength = 0;
    Double sampling  = 0;
    if(Parallel::isMaster(comm))
    {
      arcTimes.resize(arcCount);
      arcData.resize(arcCount);
      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      {
        Arc arc = instrumentFile.readArc(arcNo);
        arcTimes.at(arcNo) = arc.times();
        arcData.at(arcNo)  = arc.matrix();
        if(arc.size() > covLength)
          sampling = medianSampling(arcTimes.at(arcNo)).seconds();
        covLength = std::max(covLength, arc.size());
        countData = std::min(countData, arcData.at(arcNo).columns()-1-startData);
      }
    }
    Parallel::broadCast(arcTimes,  0, comm);
    Parallel::broadCast(arcData,   0, comm);
    Parallel::broadCast(countData, 0, comm);
    Parallel::broadCast(covLength, 0, comm);
    Parallel::broadCast(sampling,  0, comm);

    // init arc sigmas
    // ---------------
    sigma = Vector(arcCount);
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      sigma(arcNo) = 1.0;

    // init covariance function
    // ------------------------
    CosTransform = Vce::cosTransform(covLength);
    covFunc = Vce::readCovarianceFunction(FileName(), covLength, countData, sampling);
    Psd = CosTransform * covFunc.column(1, countData);

    W.resize(countData, Matrix(covLength, Matrix::SYMMETRIC));
    for(UInt idData=0; idData<countData; idData++)
    {
      for(UInt z=0; z<covLength; z++)
        for(UInt s=z; s<covLength; s++)
          W.at(idData)(z,s) += covFunc(s-z, 1+idData);
      cholesky(W.at(idData));
    }

    // =============================================

    // Iteration
    // ---------
    for(UInt iter=0; iter<iterCount; iter++)
    {
      logStatus<<"starting iteration "<<iter+1<<Log::endl;

      // solve normal equations
      // ----------------------
      if(parameter->parameterCount())
      {
        logStatus<<"accumulate system of normal equations"<<Log::endl;
        N = Matrix(countData*parameter->parameterCount(), Matrix::SYMMETRIC);
        n = Vector(countData*parameter->parameterCount());
        lPl      = 0;
        obsCount = 0;
        Parallel::forEach(arcCount, [this](UInt arcNo) {buildNormals(arcNo);}, comm);

        // collect system of normal equations
        // ----------------------------------
        logStatus<<"collect system of normal equations"<<Log::endl;
        Parallel::reduceSum(N,        0, comm);
        Parallel::reduceSum(n,        0, comm);
        Parallel::reduceSum(obsCount, 0, comm);
        Parallel::reduceSum(lPl,      0, comm);

        logStatus<<"solve system of normal equations"<<Log::endl;
        if(Parallel::isMaster(comm))
        {
          x = solve(N, n);
          logInfo<<"  aposteriori sigma = "<<sqrt((lPl-inner(x, n))/(obsCount-x.rows()))<<Log::endl;

          // N contains now the cholesky decomposition
          Wz = Vce::monteCarlo(x.rows(), 100); // monte carlo vector for VCE
          triangularSolve(1., N, Wz);
        }
        Parallel::broadCast(x,  0, comm);
        Parallel::broadCast(Wz, 0, comm);

        if(Parallel::isMaster(comm) && !fileNameSolution.empty())
        {
          logStatus<<"write solution to <"<<fileNameSolution<<">"<<Log::endl;
          writeFileMatrix(fileNameSolution, x);
        }
      } // if(parameter->parameterCount())
      Parallel::barrier(comm);

      if(!fileNameResiduals.empty())
      {
        logStatus<<"compute residuals"<<Log::endl;
        std::vector<Arc> arcs(arcCount);
        Parallel::forEach(arcs, [this](UInt arcNo) {return computeResiduals(arcNo);}, comm);

        if(Parallel::isMaster(comm))
        {
          logStatus<<"write residual file <"<<fileNameResiduals<<">"<<Log::endl;
          InstrumentFile::write(fileNameResiduals, arcs);
        }
      }

      logStatus<<"compute redundancies"<<Log::endl;
      sigmaNew = Vector(arcCount);
      ePe = redundancy = Matrix(covLength, countData);
      Parallel::forEach(arcCount, [this](UInt arcNo) {computeRedundancies(arcNo);}, comm);

      // sigmas per arc
      // --------------
      Parallel::reduceSum(sigmaNew, 0, comm);
      if(Parallel::isMaster(comm))
      {
        sigma = sigmaNew;
        const Double sigma0 = Vce::meanSigma(sigma);
        sigma *= 1./sigma0;
        logInfo<<"   sigma per arc (median): "<<sigma0<<Log::endl;

        if(!fileNameArcSigma.empty())
        {
          logStatus<<"write arc sigma file <"<<fileNameArcSigma<<">"<<Log::endl;
          writeFileMatrix(fileNameArcSigma, sigma);
        }
      }
      Parallel::broadCast(sigma, 0, comm);

      // estimate new PSD through variance component estimation
      // ------------------------------------------------------
      logStatus<<"compute psd through variance component estimation"<<Log::endl;
      Parallel::reduceSum(ePe, 0, comm);
      Parallel::reduceSum(redundancy, 0, comm);
      if(Parallel::isMaster(comm))
      {
        Double maxFactor = 0;
        Vce::estimatePsd(ePe, redundancy, Psd, maxFactor);
        logInfo<<"  max. PSD adjustment factor: "<<maxFactor<<Log::endl;
      } // if(Parallel::isMaster(comm))
      Parallel::broadCast(Psd, 0, comm);
      copy(CosTransform*Psd, covFunc.column(1, Psd.columns())); // compute new covariance function

      for(UInt idData=0; idData<countData; idData++)
      {
        W.at(idData) = Matrix(covLength, Matrix::SYMMETRIC);
        for(UInt z=0; z<covLength; z++)
          for(UInt s=z; s<covLength; s++)
              W.at(idData)(z,s) += covFunc(s-z, 1+idData);
        cholesky(W.at(idData));
      }

      if(Parallel::isMaster(comm) && !fileNameCov.empty())
      {
        logStatus<<"write covariance function file <"<<fileNameCov<<">"<<Log::endl;
        writeFileMatrix(fileNameCov, covFunc);
      }
    } // for(iter)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Instrument2CovarianceFunctionVCE::computeObservationEquation(UInt arcNo, std::vector<Matrix> &W, Vector &Wl, Matrix &WA, Matrix &WB)
{
  try
  {
    const UInt countEpoch = arcData.at(arcNo).rows();
    parameterPerArc->setInterval(arcTimes.at(arcNo).front(), arcTimes.at(arcNo).back()+medianSampling(arcTimes.at(arcNo)), TRUE);

    W.resize(countData);
    Wl = Vector(countEpoch*countData);
    WA = Matrix(Wl.rows(), countData*parameter->parameterCount());
    WB = Matrix(Wl.rows(), countData*parameterPerArc->parameterCount());

    for(UInt idData=0; idData<countData; idData++)
    {
      // observations
      copy(arcData.at(arcNo).column(1+startData+idData), Wl.row(idData*countEpoch, countEpoch));

      // design matrix
      for(UInt idEpoch=0; idEpoch<arcTimes.at(arcNo).size(); idEpoch++)
      {
        if(parameter->parameterCount())
          copy(parameter->factors(arcTimes.at(arcNo).at(idEpoch)).trans(),
               WA.slice(idData*countEpoch+idEpoch, idData*parameter->parameterCount(), 1, parameter->parameterCount()));
        if(parameterPerArc->parameterCount())
          copy(parameterPerArc->factors(arcTimes.at(arcNo).at(idEpoch)).trans(),
               WB.slice(idData*countEpoch+idEpoch, idData*parameterPerArc->parameterCount(), 1, parameterPerArc->parameterCount()));
      }

      // decorrelation
      W.at(idData) = sigma(arcNo) * this->W.at(idData).slice(0, 0, countEpoch, countEpoch);
      if(Wl.size()) triangularSolve(1., W.at(idData).trans(), Wl.row(idData*countEpoch, countEpoch));
      if(WA.size()) triangularSolve(1., W.at(idData).trans(), WA.row(idData*countEpoch, countEpoch));
      if(WB.size()) triangularSolve(1., W.at(idData).trans(), WB.row(idData*countEpoch, countEpoch));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Instrument2CovarianceFunctionVCE::buildNormals(UInt arcNo)
{
  try
  {
    std::vector<Matrix> W;
    Vector Wl;
    Matrix WA, WB;
    computeObservationEquation(arcNo, W, Wl, WA, WB);

    // eliminate arc dependent parameters
    Vector tau;
    if(WB.size())
    {
      tau = QR_decomposition(WB);
      QTransMult(WB, tau, Wl); // transform observations: l:= Q'l
      QTransMult(WB, tau, WA); // transform design matrix A:=Q'A
    }
    // use only nullspace of design matrix WB
    MatrixSlice A_bar( WA.row(WB.columns(), WA.rows()-WB.columns()) );
    MatrixSlice l_bar( Wl.row(WB.columns(), Wl.rows()-WB.columns()) );

    // build normals
    this->lPl      += quadsum(l_bar);
    this->obsCount += l_bar.rows();
    matMult(1., A_bar.trans(), l_bar, n);
    rankKUpdate(1., A_bar, N);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Instrument2CovarianceFunctionVCE::computeRedundancies(UInt arcNo)
{
  try
  {
    std::vector<Matrix> W;
    Vector Wl;
    Matrix WA, WB;
    computeObservationEquation(arcNo, W, Wl, WA, WB);

    // eliminate arc dependent parameters
    // ----------------------------------
    if(WB.size())
    {
      Vector tau = QR_decomposition(WB);
      QTransMult(WB, tau, Wl);           // transform observations: l:= Q'l
      Wl.row(0, WB.columns()).setNull(); // residuals: remove WB*x
      QMult(WB, tau, Wl);                // back transformation

      if(WA.size())
      {
        QTransMult(WB, tau, WA); // transform design matrix A:=Q'A
        WA.row(0, WB.columns()).setNull(); // residuals: remove WB*x
        QMult(WB, tau, WA); // back transformation
      }

      generateQ(WB, tau);
      WB.setType(Matrix::GENERAL);
    }

    // decorrelated residuals
    // ----------------------
    Matrix We = Wl;
    Matrix WAz(Wl.rows(), Wz.columns());
    if(WA.size()) matMult(-1., WA,  x, We);
    if(WA.size()) matMult( 1., WA, Wz, WAz);

    // Variance component estimation
    // -----------------------------
    const UInt countEpoch = Wl.rows()/countData;
    std::vector<UInt> index(countEpoch);
    for(UInt i=0; i<index.size(); i++)
      index.at(i) = i;
    Double ePeSum=0, redundancySum=0;

    for(UInt idData=0; idData<countData; idData++)
    {
      Matrix R;
      Vector WWe;
      Vce::redundancy(W.at(idData),
                      We.row(idData*countEpoch, countEpoch),
                      WAz.row(idData*countEpoch, countEpoch),
                      WB.row(idData*countEpoch, countEpoch),
                      R, WWe);
      Vce::psd(R, WWe, index, sigma(arcNo), CosTransform, Psd.column(idData),
               ePe.column(idData), redundancy.column(idData), ePeSum, redundancySum);
    }
    sigmaNew(arcNo) = sqrt(ePeSum/redundancySum) * sigma(arcNo);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Arc Instrument2CovarianceFunctionVCE::computeResiduals(UInt arcNo)
{
  try
  {
    std::vector<Matrix> W;
    Vector Wl;
    Matrix WA, WB;
    computeObservationEquation(arcNo, W, Wl, WA, WB);

    // decorrelated residuals
    if(WA.size())
      matMult(-1., WA, x, Wl);
    if(WB.size())
      reduceLeastSquaresFit(WB, Wl);

    // remove decorrelation
    const UInt countEpoch = Wl.rows()/countData;
    for(UInt idData=0; idData<countData; idData++)
      triangularMult(1., W.at(idData).trans(), Wl.row(idData*countEpoch, countEpoch));

    // observations
    Matrix data = arcData.at(arcNo);
    for(UInt idData=0; idData<countData; idData++)
      copy(Wl.row(idData*countEpoch, countEpoch), data.column(1+startData+idData));

    return Arc(arcTimes.at(arcNo), arcData.at(arcNo), arcType);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
