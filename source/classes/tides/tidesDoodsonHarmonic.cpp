/***********************************************/
/**
* @file tidesDoodsonHarmonic.cpp
*
* @brief Ocean or atmospheric tides.
* In terms of spherical harmonics.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @author Daniel Rieser
* @date 2002-12-13.
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "base/doodson.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "files/fileAdmittance.h"
#include "files/fileDoodsonHarmonic.h"
#include "classes/tides/tidesDoodsonHarmonic.h"

/***********************************************/

TidesDoodsonHarmonic::TidesDoodsonHarmonic(Config &config)
{
  try
  {
    FileName    tidesName, admittanceName;
    UInt        minDegree;
    UInt        maxDegree = INFINITYDEGREE;
    Double      factor;
    std::vector<Doodson> selectDoodson;

    renameDeprecatedConfig(config, "inputfileOcean", "inputfileTides", date2time(2020, 10, 2));

    readConfig(config, "inputfileTides",      tidesName, Config::MUSTSET,  "{groopsDataDir}/tides/oceanTide_fes2014b_n180_version20170520.dat", "");
    readConfig(config, "inputfileAdmittance", admittanceName, Config::OPTIONAL, "{groopsDataDir}/tides/oceanTide_fes2014b_admittance_linear_linear.txt", "interpolation of minor constituents");
    readConfig(config, "selectDoodson",       selectDoodson,  Config::OPTIONAL, "",    "consider only these constituents, code number (e.g. 255.555) or darwin name (e.g. M2)");
    readConfig(config, "minDegree",           minDegree,      Config::DEFAULT,  "2",   "");
    readConfig(config, "maxDegree",           maxDegree,      Config::OPTIONAL, "",    "");
    readConfig(config, "nodeCorr",            nCorr,          Config::DEFAULT,  "0",   "nodal corrections: 0-no corr, 1-IHO, 2-Schureman");
    readConfig(config, "factor",              factor,         Config::DEFAULT,  "1.0", "the result is multplied by this factor, set -1 to substract the field");
    if(isCreateSchema(config)) return;

    // read tide file
    // --------------------
    DoodsonHarmonic d;
    readFileDoodsonHarmonic(tidesName, d);
    GM = d.GM;
    R  = d.R;
    cnmCos.resize(d.doodson.size());
    snmCos.resize(d.doodson.size());
    cnmSin.resize(d.doodson.size());
    snmSin.resize(d.doodson.size());
    for(UInt i=0; i<d.doodson.size(); i++)
    {
      SphericalHarmonics harmCos = SphericalHarmonics(d.GM, d.R, d.cnmCos.at(i), d.snmCos.at(i)).get(maxDegree, minDegree);
      SphericalHarmonics harmSin = SphericalHarmonics(d.GM, d.R, d.cnmSin.at(i), d.snmSin.at(i)).get(maxDegree, minDegree);
      cnmCos.at(i) = harmCos.cnm();
      snmCos.at(i) = harmCos.snm();
      cnmSin.at(i) = harmSin.cnm();
      snmSin.at(i) = harmSin.snm();
    }

    // read admittace file
    // -------------------
    if(!admittanceName.empty())
    {
      Admittance admit;
      readFileAdmittance(admittanceName, admit);
      doodson    = admit.doodsonMinor;
      admittance = admit.admittance;

      for(UInt i=0; i<d.doodson.size(); i++)
        if(admit.doodsonMajor.at(i) != d.doodson.at(i))
          throw(Exception("different MajorTides in <"+tidesName.str()+"> and <"+admittanceName.str()+">"));
    }
    else
    {
      // unity matrix
      doodson    = d.doodson;
      admittance = Matrix(doodson.size(), doodson.size());
      for(UInt i=0; i<admittance.rows(); i++)
        admittance(i,i) = 1.0;
    }

    // Matrix with Doodson multiplicators
    doodsonMatrix = Doodson::matrix(doodson);

    admittance *= factor;

    // calculate only selected tides
    // -----------------------------
    if(selectDoodson.size()!=0)
    {
      std::vector<Bool> found(selectDoodson.size(), FALSE);
      Matrix I(admittance.rows(), Matrix::SYMMETRIC);
      for(UInt i=0; i<doodson.size(); i++)
      {
        Bool remove = TRUE;
        for(UInt k=0; k<selectDoodson.size(); k++)
          if(doodson.at(i) == selectDoodson.at(k))
          {
            found.at(k) = TRUE;
            remove      = FALSE;
            break;
          }
        if(remove)
          admittance.column(i) *= 0.;
      }
      for(UInt k=0; k<selectDoodson.size(); k++)
        if(!found.at(k))
          logWarning<<selectDoodson.at(k).code()<<" not found in admittance/tides"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix TidesDoodsonHarmonic::interpolationFactors(const Time &time) const
{
  try
  {
    Vector thetaf = doodsonMatrix * Doodson::arguments(time);
    Matrix csMinor(thetaf.rows(),2);

    if(nCorr!=0) //apply nodal corrections here
    {
      Matrix fu = Doodson::nodeCorr(doodson, time, nCorr);
      for(UInt i=0; i<thetaf.rows(); i++)
      {
        csMinor(i,0) = fu(i,0) * cos(thetaf(i)+fu(i,1));
        csMinor(i,1) = fu(i,0) * sin(thetaf(i)+fu(i,1));
      }
    }
    else
    {
      for(UInt i=0; i<thetaf.rows(); i++)
      {
        csMinor(i,0) = cos(thetaf(i));
        csMinor(i,1) = sin(thetaf(i));
      }
    }

    return ((admittance.size()) ? (admittance * csMinor) : csMinor);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics TidesDoodsonHarmonic::sphericalHarmonics(const Time &time, const Rotary3d &/*rotEarth*/, EarthRotationPtr /*rotation*/, EphemeridesPtr /*ephemerides*/, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    const Matrix csMajor = interpolationFactors(time);

    Matrix cnm(cnmCos.at(0).rows(), Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(snmCos.at(0).rows(), Matrix::TRIANGULAR, Matrix::LOWER);
    for(UInt i=0; i<csMajor.rows(); i++)
    {
      axpy(csMajor(i,0), cnmCos.at(i), cnm);
      axpy(csMajor(i,1), cnmSin.at(i), cnm);
      axpy(csMajor(i,0), snmCos.at(i), snm);
      axpy(csMajor(i,1), snmSin.at(i), snm);
    }
    return SphericalHarmonics(this->GM, this->R, cnm, snm).get(maxDegree, minDegree, GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TidesDoodsonHarmonic::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Rotary3d> &/*rotEarth*/,
                                       EarthRotationPtr /*rotation*/, EphemeridesPtr /*ephemerides*/,
                                       const std::vector<Double> &gravity, const Vector &hn, const Vector &ln,
                                       std::vector<std::vector<Vector3d>> &disp) const
{
  try
  {
    if((time.size()==0) || (point.size()==0))
      return;

    Matrix A = deformationMatrix(point, gravity, hn, ln, GM, R, cnmCos.at(0).rows()-1);
    std::vector<Matrix> xCos(cnmCos.size());
    std::vector<Matrix> xSin(snmCos.size());
    for(UInt i=0; i<xCos.size(); i++)
    {
      xCos.at(i) = A * SphericalHarmonics(GM, R, cnmCos.at(i), snmCos.at(i)).x();
      xSin.at(i) = A * SphericalHarmonics(GM, R, cnmSin.at(i), snmSin.at(i)).x();
    }

    for(UInt idEpoch=0; idEpoch<time.size(); idEpoch++)
    {
      const Matrix csMajor = interpolationFactors(time.at(idEpoch));

      Vector x(3*point.size());
      for(UInt i=0; i<csMajor.rows(); i++)
      {
        axpy(csMajor(i,0), xCos.at(i), x);
        axpy(csMajor(i,1), xSin.at(i), x);
      }

      for(UInt k=0; k<point.size(); k++)
      {
        disp.at(k).at(idEpoch).x() += x(3*k+0);
        disp.at(k).at(idEpoch).y() += x(3*k+1);
        disp.at(k).at(idEpoch).z() += x(3*k+2);
      }
    } // for(idEpoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
