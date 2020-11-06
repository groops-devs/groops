/***********************************************/
/**
* @file gravityfield2EmpiricalCovariance.cpp
*
* @brief Estimate a empircal covariance function from time series.
*
* @author Torsten Mayer-Guerr
* @date 2009-03-30
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates an spatial and temporal covariance matrix from
a time series of gravity fields.

Firstly for every time step $t_i$
a spherical harmonics vector $\M x_i$ from the time variable gravity field
is generated. The coefficients of the spherical harmonics expansion are
in the sequence given by \configClass{numbering}{sphericalHarmonicsNumberingType}.
If set the expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusivly. The coefficients are related to the
reference radius~\config{R} and the Earth gravitational constant \config{GM}.

In the next step the empirical covariance matrix is estimated according to
\begin{equation}
\M\Sigma(\Delta i)_{full} = \frac{1}{N}\sum_{i=1}^N \M x_i \M x_{i+\Delta i}^T,
\end{equation}
where $\Delta i$ is given by \config{differenceStep}.

From the diagonal elements of $\M\Sigma(\Delta i)$ the isotropic accuracies
are computed
\begin{equation}
\sigma_n^2 = \frac{1}{2n+1}\sum_{m=0}^n \sigma_{cnm}^2+\sigma_{snm}^2,
\end{equation}
and a diagonal matrix is constructed $\Sigma_{iso} = \text{diag}(\sigma_2^2,\ldots,\sigma_N^2)$.
The result is computed:
\begin{equation}
\M\Sigma(\Delta i) = \alpha_{full}\M\Sigma(\Delta i)_{full}+\alpha_{iso}\M\Sigma(\Delta i)_{iso}.
\end{equation}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Estimate a empircal covariance function from time series.
* @ingroup programsGroup */
class Gravityfield2EmpiricalCovariance
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2EmpiricalCovariance, SINGLEPROCESS, "Estimate a empircal covariance function from time series", Gravityfield, Covariance)

/***********************************************/

void Gravityfield2EmpiricalCovariance::run(Config &config)
{
  try
  {
    FileName        outputName, fileNamePotCoeff;
    UInt            minDegree, maxDegree = INFINITYDEGREE;
    Time            time;
    Double          GM, R;
    GravityfieldPtr gravityfield;
    TimeSeriesPtr   timeSeriesPtr, timesIntervalPtr;
    Bool            removeMean;
    UInt            delay;
    SphericalHarmonicsNumberingPtr numbering;
    Double          factorFullMatrix, factorIsotropic;

    readConfig(config, "outputfileCovarianceMatrix",      outputName,       Config::MUSTSET,  "",    "");
    readConfig(config, "outputfilePotentialCoefficients", fileNamePotCoeff, Config::OPTIONAL, "",    "");
    readConfig(config, "gravityfield",                    gravityfield,     Config::MUSTSET,  "",    "");
    readConfig(config, "minDegree",                       minDegree,        Config::DEFAULT,  "2",   "");
    readConfig(config, "maxDegree",                       maxDegree,        Config::MUSTSET,  "",    "");
    readConfig(config, "GM",                              GM,               Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                               R,                Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "numbering",                       numbering,        Config::MUSTSET,  "",    "numbering scheme for solution vector");
    readConfig(config, "removeMean",                      removeMean,       Config::DEFAULT,  "1",   "");
    readConfig(config, "timeSeries",                      timeSeriesPtr,    Config::MUSTSET,  "",    "sampling of the gravityfield");
    readConfig(config, "differenceStep",                  delay,            Config::DEFAULT,  "0",   "choose dt for: x,i(t) - x,j(t+dt)");
    readConfig(config, "factorFullMatrixPart",            factorFullMatrix, Config::DEFAULT,  "1",   "");
    readConfig(config, "factorIsotropicPart",             factorIsotropic,  Config::DEFAULT,  "0.1", "");
    readConfig(config, "intervals",                       timesIntervalPtr, Config::DEFAULT,  "",    "");
    if(isCreateSchema(config)) return;

    // ============================

    // init time intervals
    // -------------------
    std::vector<Time> timeSeries = timeSeriesPtr->times();
    std::vector<Time> timesInterval;
    if(timesIntervalPtr)
      timesInterval = timesIntervalPtr->times();
    if(timesInterval.size()==0)
    {
      timesInterval.push_back(timeSeries.at(0));
      timesInterval.push_back(timeSeries.back()+seconds2time(1));
    }

    // ============================

    // numbering of covariance matrix
    // ------------------------------
    std::vector< std::vector<UInt> > idxC, idxS;
    numbering->numbering(maxDegree, minDegree, idxC, idxS);
    const UInt dim = numbering->parameterCount(maxDegree, minDegree);
    Matrix CovFull(dim, dim);

    logStatus<<"estimate covariance covariance function in each interval"<<Log::endl;
    UInt countTotal = 0;
    logTimerStart;
    for(UInt idInterval=0; idInterval<timesInterval.size()-1; idInterval++)
    {
      logTimerLoop(idInterval, timesInterval.size()-1);

      // init time series
      // ----------------
      std::vector<Time> times;
      for(UInt i=0; i<timeSeries.size(); i++)
        if(timeSeries.at(i).isInInterval(timesInterval.at(idInterval), timesInterval.at(idInterval+1)))
          times.push_back(timeSeries.at(i));

      // create time series of potential coefficients
      // --------------------------------------------
      UInt count = times.size();
      std::vector<SphericalHarmonics> harm(count);
      for(UInt i=0; i<count; i++)
      {
        harm.at(i) = gravityfield->sphericalHarmonics(times.at(i), maxDegree, minDegree, GM, R);
        maxDegree = harm.at(i).maxDegree();
      }

      // remove temporal mean
      // --------------------
      if(removeMean)
      {
        SphericalHarmonics harmMean;
        for(UInt i=0; i<count; i++)
          harmMean += harm.at(i);
        harmMean *= 1./count;
        for(UInt i=0; i<count; i++)
          harm.at(i) -= harmMean;
      }

      // covariance function
      // -------------------
      count = harm.size()-delay;
      for(UInt i=0; i<count; i++)
      {
        Matrix cnm1 = harm.at(i).cnm();
        Matrix snm1 = harm.at(i).snm();
        Matrix cnm2 = harm.at(i+delay).cnm();
        Matrix snm2 = harm.at(i+delay).snm();
        // sort coefficients into vectors
        Vector x1(dim), x2(dim);
        for(UInt n=0; n<=maxDegree; n++)
        {
          if(idxC[n][0]!=NULLINDEX) x1(idxC[n][0]) = cnm1(n,0);
          if(idxC[n][0]!=NULLINDEX) x2(idxC[n][0]) = cnm2(n,0);
          for(UInt m=1; m<=n; m++)
          {
            if(idxC[n][m]!=NULLINDEX)  x1(idxC[n][m]) = cnm1(n,m);
            if(idxC[n][m]!=NULLINDEX)  x2(idxC[n][m]) = cnm2(n,m);
            if(idxS[n][m]!=NULLINDEX)  x1(idxS[n][m]) = snm1(n,m);
            if(idxS[n][m]!=NULLINDEX)  x2(idxS[n][m]) = snm2(n,m);
          }
        }

        matMult(1., x1, x2.trans(), CovFull);
      }

      countTotal += count;
      if(removeMean)
        countTotal -= 1;
    } // for(idInterval)
    logTimerLoopEnd(timesInterval.size()-1);

    // ============================

    if(delay==0)
      CovFull.setType(Matrix::SYMMETRIC);
    CovFull *= 1./countTotal;
    Matrix Cov = CovFull;
    Cov *= factorFullMatrix;

    if(factorIsotropic!=0)
    {
      Vector kn(maxDegree+1);
      for(UInt n=0; n<=maxDegree; n++)
      {
        if(idxC[n][0]!=NULLINDEX) kn(n) += CovFull(idxC[n][0],idxC[n][0]);
        for(UInt m=1; m<=n; m++)
        {
          if(idxC[n][m]!=NULLINDEX) kn(n) += CovFull(idxC[n][m],idxC[n][m]);
          if(idxS[n][m]!=NULLINDEX) kn(n) += CovFull(idxS[n][m],idxS[n][m]);
        }
      }
      for(UInt n=0; n<=maxDegree; n++)
      {
        if(idxC[n][0]!=NULLINDEX) Cov(idxC[n][0],idxC[n][0]) += factorIsotropic * kn(n)/(2*n+1);
        for(UInt m=1; m<=n; m++)
        {
          if(idxC[n][m]!=NULLINDEX) Cov(idxC[n][m],idxC[n][m]) += factorIsotropic * kn(n)/(2*n+1);
          if(idxS[n][m]!=NULLINDEX) Cov(idxS[n][m],idxS[n][m]) += factorIsotropic * kn(n)/(2*n+1);
        }
      }
    }

    // ============================

    // Write
    // -----
    logStatus << "save covariance matrix to <"<<outputName<<">"<< Log::endl;
    if(delay == 0)
      Cov.setType(Matrix::SYMMETRIC);
    writeFileMatrix(outputName, Cov);

    logInfo<<"  number of fields     : "<<countTotal<<Log::endl;
    logInfo<<"  minDegree            : "<<minDegree<< Log::endl;
    logInfo<<"  maxDegree            : "<<maxDegree<< Log::endl;
    logInfo<<"  size of state vector : "<<dim<< Log::endl;
    logInfo<<"  number of unknowns   : "<<dim*dim<< Log::endl;
    logInfo<<"  redundancy           : "<<countTotal/static_cast<Double>(dim) << Log::endl;

    // Spherical Harmonics
    // -------------------
    if(!fileNamePotCoeff.empty())
    {
      logStatus<<"save potential coefficients to <"<<fileNamePotCoeff<<">"<< Log::endl;
      Matrix cnm(maxDegree+1, Matrix::SYMMETRIC);
      Matrix snm(maxDegree+1, Matrix::SYMMETRIC);
      Matrix sigma2cnm(maxDegree+1, Matrix::SYMMETRIC);
      Matrix sigma2snm(maxDegree+1, Matrix::SYMMETRIC);

      for(UInt n=0; n<=maxDegree; n++)
      {
        if(idxC[n][0]!=NULLINDEX) cnm(n,0) = sqrt(Cov(idxC[n][0],idxC[n][0]));
        if(idxC[n][0]!=NULLINDEX) sigma2cnm(n,0) = Cov(idxC[n][0],idxC[n][0]);
        for(UInt m=1; m<=n; m++)
        {
          if(idxC[n][m]!=NULLINDEX) cnm(n,m) = sqrt(Cov(idxC[n][m],idxC[n][m]));
          if(idxS[n][m]!=NULLINDEX) snm(n,m) = sqrt(Cov(idxS[n][m],idxS[n][m]));
          if(idxC[n][m]!=NULLINDEX) sigma2cnm(n,m) = Cov(idxC[n][m],idxC[n][m]);
          if(idxS[n][m]!=NULLINDEX) sigma2snm(n,m) = Cov(idxS[n][m],idxS[n][m]);
        }
      }

      writeFileSphericalHarmonics(fileNamePotCoeff, SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm));
    }

    // Test for positive definite
    // --------------------------
    try
    {
      if(delay==0)
        cholesky(Cov);
    }
    catch(std::exception &/*e*/)
    {
      logWarning<<"covariance matrix is not positive definite!"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
