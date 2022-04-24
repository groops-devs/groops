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
  InFileArchive       file;
  UInt                splineDegree_;
  std::vector<Time>   times_;
  Double              GM, R;
  UInt                maxDegree_, minDegree_;
  UInt                indexFile, indexData;
  std::streampos      seekStart;
  std::streamoff      seekSize;
  std::vector<SphericalHarmonics> harmonics;

public:
  InFileTimeSplinesGravityfield() {}
  InFileTimeSplinesGravityfield(const FileName &name, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0) {open(name, maxDegree, minDegree);}
 ~InFileTimeSplinesGravityfield() {close();}

  void  open(const FileName &name, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0);
  void  close();

  UInt splineDegree() const {return splineDegree_;}
  UInt nodeCount()    const {return times_.size()+splineDegree_-1;}  //!< number of spline nodal points (agree with times().size() for spline degree 1)
  const std::vector<Time> &times() const {return times_;}

  /// at spline nodal point
  SphericalHarmonics sphericalHarmonics(UInt idNode);

  /// interpolated at @p time.
  SphericalHarmonics sphericalHarmonics(const Time &time, Double factor=1.0);
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
  UInt                splineDegree_;
  std::vector<Time>   times_;
  Double              GM_, R_;
  Bool                mustCut;
  UInt                maxDegree_, minDegree_;
  UInt                indexFile, indexData;
  std::streampos      seekStart;
  std::streamoff      seekSize;
  std::vector<Matrix> cov;

public:
  InFileTimeSplinesCovariance() {}
  InFileTimeSplinesCovariance(const FileName &name, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0) {open(name, maxDegree, minDegree);}
 ~InFileTimeSplinesCovariance() {close();}

  void  open(const FileName &name, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0);
  void  close();

  UInt   splineDegree() const {return splineDegree_;}
  UInt   nodeCount()    const {return times_.size()+splineDegree_-1;}  //!< number of spline nodal points (agree with times().size() for spline degree 1)
  Double GM()           const {return GM_;}
  Double R()            const {return R_;}
  UInt   maxDegree()    const {return maxDegree_;}
  UInt   minDegree()    const {return minDegree_;}
  const std::vector<Time> &times() const {return times_;}

  /// at spline nodal point
  Matrix covariance(UInt idNode);

  /** @brief interpolated covariances at time t.
  * The result is a full covariance matrix or a vector containing only the variances
  * of a spherical harmonics expansion.
  * The coefficients are given in a degree wise sequence with alternating sin, cos:
  * ..., c20,c21,s21,c22,s22,..., beginning with c00.  */
  Matrix covariance(const Time &time, Double factor, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0);
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
