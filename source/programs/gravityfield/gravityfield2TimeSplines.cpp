/***********************************************/
/**
* @file gravityfield2TimeSplines.cpp
*
* @brief Estimate splines in time domain from a time variable gravity field.
*
* @author Torsten Mayer-Guerr
* @date 2006-09-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates splines in time domain from a time variable gravity field
and writes \configFile{outputfileTimeSplines}{timeSplinesGravityField}.
The \configClass{gravityfield}{gravityfieldType} is sampled at \configClass{sampling}{timeSeriesType}, converted to potential coefficients
in the range between \config{minDegree} and \config{maxDegree} inclusively.
The time series of spherical harmonics can be temporal filtered with \configClass{temporalFilter}{digitalFilterType}.

In the next step temporal splines with \config{splineDegree} and nodal points given
at \configClass{splineTimeSeries}{timeSeriesType} are adjusted to the time series in a least squares sense.
This is very fast for block means (splineDegree = 0) but for other splines a large systems of equations
must be solved. In the adjustment process the time series of gravity fields can be interpreted as samples
at the given times or as continuous function with linear behaviour between sampled points (\config{linearInterpolation}).

To combine a series of potential coefficients to a spline file with block means (splineDegree = 0)
use the fast \program{PotentialCoefficients2BlockMeanTimeSplines} instead.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileTimeSplinesGravityfield.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Estimate splines in time domain from a time variable gravity field.
* @ingroup programsGroup */
class Gravityfield2TimeSplines
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2TimeSplines, SINGLEPROCESS, "Estimate splines in time domain from a time variable gravity field", Gravityfield, TimeSplines)

/***********************************************/

void Gravityfield2TimeSplines::run(Config &config)
{
  try
  {
    FileName         outputName;
    UInt             minDegree, maxDegree = INFINITYDEGREE;
    Time             time;
    Double           GM, R;
    GravityfieldPtr  gravityfield;
    TimeSeriesPtr    timeSeriesSplines, timeSeriesObs;
    Bool             removeMean, interpolateLinear;
    UInt             splineDegree;
    DigitalFilterPtr digitalFilter;

    readConfig(config, "outputfileTimeSplines", outputName,        Config::MUSTSET,  "",  "");
    readConfig(config, "gravityfield",          gravityfield,      Config::MUSTSET,  "",  "");
    readConfig(config, "temporalFilter",        digitalFilter,     Config::OPTIONAL, "",  "filter sampled gravity field in time");
    readConfig(config, "minDegree",             minDegree,         Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",             maxDegree,         Config::OPTIONAL, "",  "");
    readConfig(config, "GM",                    GM,                Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                     R,                 Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "sampling",              timeSeriesObs,     Config::MUSTSET,  "",  "gravity field is sampled at these times");
    readConfig(config, "removeMean",            removeMean,        Config::DEFAULT,  "0", "remove the temporal mean of the series before estimating the splines");
    readConfig(config, "linearInterpolation",   interpolateLinear, Config::DEFAULT,  "0", "assume linear behavior between sampled points");
    readConfig(config, "splineDegree",          splineDegree,      Config::MUSTSET,  "",  "degree of splines");
    readConfig(config, "splineTimeSeries",      timeSeriesSplines, Config::MUSTSET,  "",  "nodal points of splines in time domain");
    if(isCreateSchema(config)) return;

    // ============================================

    // sample gravity field to create time series of potential coefficients
    // --------------------------------------------------------------------
    logStatus<<"generate time series of potential coefficients"<<Log::endl;
    std::vector<Time> timesObs = timeSeriesObs->times();

    SphericalHarmonics shc = gravityfield->sphericalHarmonics(timesObs.front(), maxDegree, minDegree, GM, R);
    maxDegree = shc.maxDegree();
    Matrix sample(timesObs.size(), (maxDegree+1)*(maxDegree+1));

    logTimerStart;
    for(UInt i=0; i<timesObs.size(); i++)
    {
      logTimerLoop(i,timesObs.size());
      copy(gravityfield->sphericalHarmonics(timesObs.at(i), maxDegree, minDegree, GM, R).x().trans(), sample.row(i));
    }
    logTimerLoopEnd(timesObs.size());

    // remove temporal mean
    if(removeMean)
    {
      logStatus<<"remove temporal mean"<<Log::endl;
      for(UInt k=0; k<sample.columns(); k++)
        sample.column(k) -= mean(sample.column(k));
    }

    // filter time series
    if(digitalFilter)
    {
      logStatus<<"apply temporal filter"<<Log::endl;
      sample = digitalFilter->filter(sample);
    }

    // ============================================

    // test time series
    // ----------------
    std::vector<Time> timesSplines = timeSeriesSplines->times();
    const UInt nodeCount = timesSplines.size()-1+splineDegree;
    if(timesSplines.size()<2)
      throw(Exception("2 points in time must be given at least"));

    // result
    // ------
    Matrix x(nodeCount, sample.columns());

    // Test special cases
    // ------------------
    if((splineDegree==0) && (!interpolateLinear)) // block mean time splines
    {
      UInt idxObs  = 0;
      for(UInt i=0; i<nodeCount; i++)
      {
        UInt count = 0;
        for(; (idxObs<timesObs.size()) && (timesObs.at(idxObs)<timesSplines.at(i+1)); idxObs++)
        {
          if(timesObs.at(idxObs)<timesSplines.at(i))
            continue;
          x.row(i) += sample.row(idxObs);
          count++;
        }
        x.row(i) *= 1./count;
      }
    }
    else if((splineDegree==0) && interpolateLinear) // block mean time splines
    {
      UInt idxObs  = 0;
      for(UInt i=0; i<nodeCount; i++)
      {
        const Double T = (timesSplines.at(i+1)-timesSplines.at(i)).mjd();
        for(; (idxObs+1<timesObs.size()) && (timesObs.at(idxObs)<timesSplines.at(i+1)); idxObs++)
        {
          const Double tau1    = (timesObs.at(idxObs)  -timesSplines.at(i)).mjd()/T;
          const Double tau2    = (timesObs.at(idxObs+1)-timesSplines.at(i)).mjd()/T;
          const Double dtau    = tau2-tau1;
          const Double factor1 = 1./dtau * ((tau2*std::min(1.,tau2) - 0.5*pow(std::min(1.,tau2),2)) - (tau2*std::max(0.,tau1) - 0.5*pow(std::max(0.,tau1),2)));
          const Double factor2 = 1./dtau * ((0.5*pow(std::min(1.,tau2),2) - tau1*std::min(1.,tau2)) - (0.5*pow(std::max(0.,tau1),2) - tau1*std::max(0.,tau1)));
          x.row(i) += factor1 * sample.row(idxObs) + factor2 * sample.row(idxObs+1);
        }
      }
    }
    else if((splineDegree==1) && (timesObs.size() == timesSplines.size())) // linear splines -> sampling points == splines values
    {
      x = sample;
    }
    else // high order splines (general case)
    {
      // coefficients of splines polynomials phi_i(tau) = sum_n coeff_{i,n}*tau^n
      // ------------------------------------------------------------------------
      Matrix coeff(splineDegree+1, splineDegree+1);
      switch(splineDegree)
      {
        case 0: // Blockmittel
          coeff(0,0) = 1.0;
          break;
        case 1: // Linear
          coeff(1,1) =  1.0;                   // t
          coeff(0,1) = -1.0; coeff(0,0) = 1.0; // 1-t
          break;
        case 2: // Quadratic
          coeff(2,2) =  0.5;                                       //  0.5*t^2
          coeff(1,2) = -1.0; coeff(1,1) =  1.0; coeff(1,0) = 0.5;  // -1.0*t^2 + t + 0.5
          coeff(0,2) =  0.5; coeff(0,1) = -1.0; coeff(0,0) = 0.5;  //  0.5*t^2 - t + 0.5
          break;
        case 3: // Cubic
          coeff(3,3) = +1./6.;                                                                 // +1./6.*t^3
          coeff(2,3) = -3./6.; coeff(2,2) = +3./6.; coeff(2,1) = +3./6.; coeff(2,0) = +1./6.;  // -3./6.*t^3 + 3./6.*t^2 + 3./6.*t + 1./6.;
          coeff(1,3) = +3./6.; coeff(1,2) = -6./6.; coeff(1,1) = +0./6.; coeff(1,0) = +4./6.;  // +3./6.*t^3 - 6./6.*t^2 + 0./6.*t + 4./6.;
          coeff(0,3) = -1./6.; coeff(0,2) = +3./6.; coeff(0,1) = -3./6.; coeff(0,0) = +1./6.;  // -1./6.*t^3 - 3./6.*t^2 - 3./6.*t + 1./6.;
          break;
        default:
          throw(Exception("degree of spline not implemented"));
      }

      logStatus<<"least squares adjustment"<<Log::endl;
      logInfo<<"  size of design matrix "<<timesObs.size()<<" x "<<nodeCount<<" = "<<8.*timesObs.size()*nodeCount/1024./1024.<<" MB"<<Log::endl;
      Matrix At(nodeCount, timesObs.size());  // transpose of design matrix
      Matrix N(nodeCount, Matrix::SYMMETRIC); // normal matrix

      if(!interpolateLinear)
      {
        UInt idxObs  = 0;
        for(UInt i=0; i<timesSplines.size()-1; i++)
        {
          const Double T = (timesSplines.at(i+1)-timesSplines.at(i)).mjd();
          for(; (idxObs<timesObs.size()) && (timesObs.at(idxObs)<timesSplines.at(i+1)); idxObs++)
          {
            if(timesObs.at(idxObs)<timesSplines.at(i))
              continue;
            const Double tau = (timesObs.at(idxObs)-timesSplines.at(i)).mjd()/T;
            for(UInt k=0; k<=splineDegree; k++)
              for(UInt n=0; n<=splineDegree; n++)
                At(i+k, idxObs) += coeff(k,n) * std::pow(tau, n);

            // normal matrix
            for(UInt k1=0; k1<=splineDegree; k1++)
              for(UInt k2=0; k2<=splineDegree; k2++)
                N(i+k1, i+k2) += At(i+k1,idxObs) * At(i+k2,idxObs);
          }
        }
      }
      else
      {
        // interpolateLinear
        UInt idxObs = 0;
        for(UInt i=0; i<timesSplines.size()-1; i++)
        {
          const Double T = (timesSplines.at(i+1)-timesSplines.at(i)).mjd();
          while((idxObs+1<timesObs.size()) && (timesObs.at(idxObs)<timesSplines.at(i+1)))
          {
            const Double tau1 = (timesObs.at(idxObs)  -timesSplines.at(i)).mjd()/T;
            const Double tau2 = (timesObs.at(idxObs+1)-timesSplines.at(i)).mjd()/T;
            const Double dtau = tau2-tau1;

            // linear interpolation between sampling points
            // and integration of the scalar product of the line and the basis functions phi=sum_n coeff_n*tau^n
            for(UInt z=0; z<=splineDegree; z++)
              for(UInt n=0; n<=splineDegree; n++)
              {
                At(i+z,idxObs)   += T/dtau * coeff(z,n) * ((tau2*pow(std::min(1.,tau2),n+1)/(n+1) - pow(std::min(1.,tau2),n+2)/(n+2)) - (tau2*pow(std::max(0.,tau1),n+1)/(n+1) - pow(std::max(0.,tau1),n+2)/(n+2)));
                At(i+z,idxObs+1) += T/dtau * coeff(z,n) * ((pow(std::min(1.,tau2),n+2)/(n+2) - tau1*pow(std::min(1.,tau2),n+1)/(n+1)) - (pow(std::max(0.,tau1),n+2)/(n+2) - tau1*pow(std::max(0.,tau1),n+1)/(n+1)));
              }

            if(timesObs.at(idxObs+1)>timesSplines.at(i+1))
              break;
            idxObs++;
          }

          // normal matrix
          // Integral of the product of the basis functions <phi_z, phi_s>
          for(UInt z=0; z<=splineDegree; z++)
            for(UInt s=0; s<=splineDegree; s++)
              for(UInt n1=0; n1<=splineDegree; n1++)
                for(UInt n2=0; n2<=splineDegree; n2++)
                  N(i+z, i+s) += T/(n1+n2+1) * coeff(z,n1) * coeff(s,n2);
        }
      } // if(interpolateLinear)

      // Solve the equation system
      // -----------------------------
      logStatus<<"solve the equation system"<<Log::endl;
      solveInPlace(N, At);
      x = At * sample;
    } // if(splineDegree>0)

    // ============================================

    // sort into coefficient triangles
    // -------------------------------
    std::vector<Matrix> cnm(nodeCount, Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
    std::vector<Matrix> snm(nodeCount, Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
    for(UInt i=0; i<nodeCount; i++)
    {
      UInt idx = 0;
      for(UInt n=0; n<=maxDegree; n++)
      {
        cnm.at(i)(n,0) = x(i, idx++);
        for(UInt m=1; m<=n; m++)
        {
          cnm.at(i)(n,m) = x(i, idx++);
          snm.at(i)(n,m) = x(i, idx++);
        }
      }
    }

    // write time splines
    // ------------------
    logStatus<<"write time splines to file <"<<outputName<<">"<<Log::endl;
    writeFileTimeSplinesGravityfield(outputName, GM, R, splineDegree, timesSplines, cnm, snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
