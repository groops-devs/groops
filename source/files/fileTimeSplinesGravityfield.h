/***********************************************/
/**
* @file fileTimeSplinesGravityfield.h
*
* @brief time variable gravity field represented by splines in time domain.
*
* @author Torsten Mayer-Guerr
* @date 2004-11-29
*
*/
/***********************************************/

#ifndef __GROOPS_FILETIMESPLINESGRAVITYFIELD__
#define __GROOPS_FILETIMESPLINESGRAVITYFIELD__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_TimeSplinesGravityField
static const char *docstringTimeSplinesGravityField = R"(
Temporal changing gravity field, parametrized as spherical harmonics in the spatial domain and
parametrized as basis splines in the time domain (see~\reference{basis splines}{fundamentals.basisSplines}).
It is evaluated with \configClass{gravityfield:timeSplines}{gravityfieldType:timeSplines}.

See also: \program{Gravityfield2TimeSplines}, \program{PotentialCoefficients2BlockMeanTimeSplines}.
)";
#endif

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_TimeSplinesCovariance
static const char *docstringTimeSplinesCovariance = R"(
Stores covariance information for \file{TimeSplinesGravityField}{timeSplinesGravityField}.
It can be the variances of the potential coefficients or the full covariance matrix for each
temporal nodal point.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_TIMESPLINESGRAVITYFIELD_TYPE = "timeSplinesGravityField";
const char *const FILE_TIMESPLINESCOVARIANCE_TYPE   = "timeSplinesCovariance";

/***** CLASS ***********************************/

/** @brief Time variable gravity field represented by splines in time domain. */
class InFileTimeSplinesGravityfield
{
  InFileArchive      file;
  Double             GM, R;
  UInt               count, index;
  UInt               degree;
  UInt               _maxDegree, _minDegree;
  std::vector<Time>  times;
  std::vector<SphericalHarmonics> harmonics;
  std::streampos     dataStart;
  std::streamoff     dataEpochSize;

public:
  InFileTimeSplinesGravityfield() {}
  InFileTimeSplinesGravityfield(const FileName &name, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0) {open(name, maxDegree, minDegree);}
 ~InFileTimeSplinesGravityfield() {close();}

  void  open(const FileName &name, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0);
  void  close();

  SphericalHarmonics sphericalHarmonics(const Time &time, Double factor=1.0);

  UInt maxDegree() const {return _maxDegree;}
  UInt minDegree() const {return _minDegree;}
};

/***** FUNCTIONS *******************************/

/** @brief Write a TimeSplinesGravityfieldFile. */
void writeFileTimeSplinesGravityfield(const FileName &fileName,
                                      Double GM, Double R, UInt splineDegree,
                                      const std::vector<Time> &times,
                                      const std::vector<Matrix> &cnm,
                                      const std::vector<Matrix> &snm);

/** @brief Read a TimeSplinesGravityfieldFile. */
void readFileTimeSplinesGravityfield(const FileName &fileName,
                                     Double &GM, Double &R, UInt &splineDegree,
                                     std::vector<Time> &times,
                                     std::vector<Matrix> &cnm,
                                     std::vector<Matrix> &snm);

/***** CLASS ***********************************/

/** @brief Time variable covariances represented by splines in time domain. */
class InFileTimeSplinesCovariance
{
  InFileArchive       file;
  Double              _GM, _R;
  UInt                _maxDegree, _minDegree;
  Bool                mustCut;
  UInt                count, index;
  UInt                degree;
  std::vector<Time>   times;
  std::vector<Matrix> cov;

public:
  InFileTimeSplinesCovariance() {}
  InFileTimeSplinesCovariance(const FileName &name, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0) {open(name, maxDegree, minDegree);}
 ~InFileTimeSplinesCovariance() {close();}

  void  open(const FileName &name, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0);
  void  close();

  /** @brief interpolated covariances at time t.
  * The result is a full covariance matrix or a vector containing only the variances
  * of a spherical harmonics expansion.
  * The coefficients are given in a degree wise sequence with alternating sin, cos:
  * ..., c20,c21,s21,c22,s22,..., beginning with c00.  */
  Matrix covariance(const Time &time, Double factor, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0);

  Double GM()        const {return _GM;}
  Double R()         const {return _R;}
  UInt   maxDegree() const {return _maxDegree;}
  UInt   minDegree() const {return _minDegree;}
};

/***** FUNCTIONS *******************************/

/** @brief Write a TimeSplinesCovarianceFile. */
void writeFileTimeSplinesCovariance(const FileName &fileName,
                                    Double GM, Double R, UInt minDegree, UInt maxDegree, UInt splineDegree,
                                    const std::vector<Time> &times,
                                    const std::vector<Matrix> &sigma2);


/** @brief Read a TimeSplinesCovarianceFile. */
void readFileTimeSplinesCovariance(const FileName &fileName,
                                   Double &GM, Double &R, UInt &minDegree, UInt &maxDegree, UInt &splineDegree,
                                   std::vector<Time> &times,
                                   std::vector<Matrix> &sigma2);
/// @}

/***********************************************/

#endif /* __GROOPS TIMESPLINESGRAVITYFIELD */
