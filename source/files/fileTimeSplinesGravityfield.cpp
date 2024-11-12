/***********************************************/
/**
* @file fileTimeSplinesGravityfield.cpp
*
* @brief time variable gravity field represented by splines in time domain.
*
* @author Torsten Mayer-Guerr
* @date 2004-11-29
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_TimeSplinesGravityField
#define DOCSTRING_FILEFORMAT_TimeSplinesCovariance

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "base/basisSplines.h"
#include "inputOutput/fileArchive.h"
#include "inputOutput/logging.h"
#include "files/fileFormatRegister.h"
#include "files/fileTimeSplinesGravityfield.h"

GROOPS_REGISTER_FILEFORMAT(TimeSplinesGravityField, FILE_TIMESPLINESGRAVITYFIELD_TYPE)
GROOPS_REGISTER_FILEFORMAT(TimeSplinesCovariance,   FILE_TIMESPLINESCOVARIANCE_TYPE)

/***** FUNCTIONS *******************************/

void InFileTimeSplinesGravityfield::open(const FileName &name, UInt maxDegree, UInt minDegree)
{
  try
  {
    maxDegree_ = maxDegree;
    minDegree_ = minDegree;

    close();
    if(!name.empty())
    {
      UInt epochCount;
      file.open(name, FILE_TIMESPLINESGRAVITYFIELD_TYPE, FILE_TIMESPLINESGRAVITYFIELD_VERSION);
      file>>nameValue("GM",        GM);
      file>>nameValue("R",         R);
      file>>nameValue("degree",    splineDegree_);
      file>>nameValue("timeCount", epochCount);
      times_.resize(epochCount);
      for(UInt i=0; i<times_.size(); i++)
        file>>nameValue("time", times_.at(i));

      // If we have a binary file, we can efficiently seek to the needed position in the file.
      if(file.canSeek())
        seekStart = file.position();
      seekSize = 0;

      indexFile = 0;
      indexData = nodeCount();
      harmonics.resize(splineDegree_+1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InFileTimeSplinesGravityfield::close()
{
  file.close();
  times_.clear();
  splineDegree_ = 0;
}

/***********************************************/

SphericalHarmonics InFileTimeSplinesGravityfield::sphericalHarmonics(UInt idNode)
{
  try
  {
    if(file.fileName().empty())
      throw(Exception("no file open"));

    // must restart?
    if(indexFile > idNode)
    {
      if(file.canSeek() && seekSize)
        indexFile = 0;
      else
        open(file.fileName());
    }

    SphericalHarmonics harmonics;
    while(indexFile <= idNode)
    {
      // Seek to appropriate interval.
      if(file.canSeek() && seekSize)
      {
        file.seek(seekStart + static_cast<std::streamoff>(idNode * seekSize));
        indexFile = idNode;
      }

      Matrix cnm, snm;
      file>>beginGroup("node");
      file>>nameValue("cnm", cnm);
      file>>nameValue("snm", snm);
      file>>endGroup("node");
      harmonics = SphericalHarmonics(GM, R, cnm, snm).get(maxDegree_, minDegree_);

      if(file.canSeek() && (indexFile == 0))
        seekSize = file.position() - seekStart;
      indexFile++;
    }

    return harmonics;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics InFileTimeSplinesGravityfield::sphericalHarmonics(const Time &time, Double factor)
{
  try
  {
    if(file.fileName().empty())
      throw(Exception("no file open"));
    if((time < times_.front()) || (time > times_.back()))
      throw(Exception(time.dateTimeStr()+" outside interval ["+times_.front().dateTimeStr()+", "+times_.back().dateTimeStr()+"] of <"+file.fileName().str()+">"));

    // find time interval and read missing nodes
    const UInt idx   = std::min(static_cast<UInt>(std::distance(times_.begin(), std::upper_bound(times_.begin(), times_.end(), time))), times_.size()-1)-1;
    const UInt start = (idx > indexData) ? std::max(idx, indexData+harmonics.size()) : idx;
    const UInt end   = (idx > indexData) ? idx+harmonics.size() : std::min(idx+harmonics.size(), indexData);
    for(UInt i=start; i<end; i++)
      harmonics.at(i%harmonics.size()) = sphericalHarmonics(i);
    indexData = idx;

    // spline interpolation in interval
    const Double       t     = (time-times_.at(idx)).mjd()/(times_.at(idx+1)-times_.at(idx)).mjd();
    const Vector       coeff = BasisSplines::compute(t, splineDegree_);
    SphericalHarmonics harm  = factor * coeff(0) * harmonics.at(indexData%harmonics.size());
    for(UInt i=1; i<coeff.rows(); i++)
      harm += factor * coeff(i) * harmonics.at((indexData+i)%harmonics.size());
    return harm;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileTimeSplinesGravityfield(const FileName &fileName,
                                      Double GM, Double R, UInt splineDegree,
                                      const std::vector<Time> &times,
                                      const std::vector<Matrix> &cnm,
                                      const std::vector<Matrix> &snm)
{
  try
  {
    UInt nodeCount = times.size()+splineDegree-1;
    if((nodeCount != cnm.size()) || (nodeCount != snm.size()))
      throw(Exception("size of series of cnm,snm does not mach to number of times."));

    OutFileArchive file(fileName, FILE_TIMESPLINESGRAVITYFIELD_TYPE, FILE_TIMESPLINESGRAVITYFIELD_VERSION);
    file<<nameValue("GM",        GM);
    file<<nameValue("R",         R);
    file<<nameValue("degree",    splineDegree);
    file<<nameValue("timeCount", times.size());
    for(UInt i=0; i<times.size(); i++)
      file<<nameValue("time", times.at(i));
    for(UInt i=0; i<nodeCount; i++)
    {
      file<<beginGroup("node");
      file<<nameValue("cnm", cnm.at(i));
      file<<nameValue("snm", snm.at(i));
      file<<endGroup("node");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileTimeSplinesGravityfield(const FileName &fileName,
                                     Double &GM, Double &R, UInt &splineDegree,
                                     std::vector<Time> &times,
                                     std::vector<Matrix> &cnm,
                                     std::vector<Matrix> &snm)
{
  try
  {
    InFileArchive file(fileName, FILE_TIMESPLINESGRAVITYFIELD_TYPE, FILE_TIMESPLINESGRAVITYFIELD_VERSION);
    file>>nameValue("GM",        GM);
    file>>nameValue("R",         R);
    file>>nameValue("degree",    splineDegree);
    UInt timeCount;
    file>>nameValue("timeCount", timeCount);
    times.resize(timeCount);
    for(UInt i=0; i<timeCount; i++)
      file>>nameValue("time", times.at(i));
    const UInt nodeCount = timeCount+splineDegree-1;
    cnm.resize(nodeCount);
    snm.resize(nodeCount);
    for(UInt i=0; i<nodeCount; i++)
    {
      file>>beginGroup("node");
      file>>nameValue("cnm", cnm.at(i));
      file>>nameValue("snm", snm.at(i));
      file>>endGroup("node");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void InFileTimeSplinesCovariance::open(const FileName &name, UInt maxDegree, UInt minDegree)
{
  try
  {
    close();
    if(!name.empty())
    {
      UInt epochCount;
      file.open(name, FILE_TIMESPLINESCOVARIANCE_TYPE, FILE_TIMESPLINESCOVARIANCE_VERSION);
      file>>nameValue("GM",        GM_);
      file>>nameValue("R",         R_);
      file>>nameValue("minDegree", minDegree_);
      file>>nameValue("maxDegree", maxDegree_);
      file>>nameValue("degree",    splineDegree_);
      file>>nameValue("timeCount", epochCount);
      times_.resize(epochCount);
      for(UInt i=0; i<times_.size(); i++)
        file>>nameValue("time", times_.at(i));

      mustCut    = (maxDegree < maxDegree_);
      maxDegree_ = std::min(maxDegree, maxDegree_);
      minDegree_ = std::max(minDegree, minDegree_);

      // If we have a binary file, we can efficiently seek to the needed position in the file.
      if(file.canSeek())
        seekStart = file.position();
      seekSize = 0;

      indexFile = 0;
      indexData = nodeCount();;
      cov.resize(splineDegree_+1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InFileTimeSplinesCovariance::close()
{
  file.close();
  times_.clear();
  splineDegree_ = 0;
}

/***********************************************/

Matrix InFileTimeSplinesCovariance::covariance(UInt idNode)
{
  try
  {
    if(file.fileName().empty())
      throw(Exception("no file open"));

    // must restart?
    if(indexFile > idNode)
    {
      if(file.canSeek() && seekSize)
        indexFile = 0;
      else
        open(file.fileName(), maxDegree(), minDegree());
    }

    Matrix C;
    while(indexFile <= idNode)
    {
      // Seek to appropriate interval.
      if(file.canSeek() && seekSize)
      {
        file.seek(seekStart + static_cast<std::streamoff>(idNode * seekSize));
        indexFile = idNode;
      }

      file>>nameValue("node", C);
      if(mustCut)
      {
        if(C.getType() == Matrix::SYMMETRIC) // full covariance matrix?
          C = C.slice(0,0,(maxDegree_+1)*(maxDegree_+1),(maxDegree_+1)*(maxDegree_+1));
        else
          C = C.row(0,(maxDegree_+1)*(maxDegree_+1));
      }

      if(file.canSeek() && (indexFile == 0))
        seekSize = file.position() - seekStart;
      indexFile++;
    }

    return C;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix InFileTimeSplinesCovariance::covariance(const Time &time, Double factor, UInt maxDegree, UInt minDegree, Double GM, Double R)
{
  try
  {
    if(file.fileName().empty())
      throw(Exception("no file open"));
    if((time < times_.front()) || (time > times_.back()))
      throw(Exception(time.dateTimeStr()+" outside interval ["+times_.front().dateTimeStr()+", "+times_.back().dateTimeStr()+"] of <"+file.fileName().str()+">"));

    // find time interval and read missing nodes
    const UInt idx   = std::min(static_cast<UInt>(std::distance(times_.begin(), std::upper_bound(times_.begin(), times_.end(), time))), times_.size()-1)-1;
    const UInt start = (idx > indexData) ? std::max(idx, indexData+cov.size()) : idx;
    const UInt end   = (idx > indexData) ? idx+cov.size() : std::min(idx+cov.size(), indexData);
    for(UInt i=start; i<end; i++)
      cov.at(i%cov.size()) = covariance(i);
    indexData = idx;

    // spline interpolation in interval
    const Double  t     = (time-times_.at(idx)).mjd()/(times_.at(idx+1)-times_.at(idx)).mjd();
    const Vector  coeff = BasisSplines::compute(t, splineDegree_);
    Matrix        C     = factor * coeff(0) * cov.at(indexData%cov.size());
    for(UInt i=1; i<coeff.rows(); i++)
      axpy(factor * coeff(i), cov.at((indexData+i)%cov.size()), C);

    if(GM <= 0)  GM = GM_;
    if(R  <= 0)  R  = R_;
    if(maxDegree == INFINITYDEGREE) maxDegree = maxDegree_;
    minDegree = std::min(minDegree, maxDegree);

    if(GM_/GM*factor != 1.0)
      C *= std::pow(GM_/GM*factor, 2);

    // minDegree
    // ---------
    if(C.getType() == Matrix::SYMMETRIC) // full covariance matrix?
    {
      for(UInt idx=0; idx<std::pow(std::max(minDegree_, minDegree),2); idx++)
      {
        C.row(idx)    *= 0;
        C.column(idx) *= 0;
      }
    }
    else
      for(UInt idx=0; idx<std::pow(std::max(minDegree_, minDegree),2); idx++)
        C.row(idx) *= 0;

    // Quick return possible?
    // ----------------------
    if((R == R_) && (maxDegree == maxDegree_))
      return C;

    // full covariance matrix?
    // -----------------------
    if(C.getType() == Matrix::SYMMETRIC)
    {
      if(maxDegree < maxDegree_)
        C = C.slice(0, 0, (maxDegree+1)*(maxDegree+1), (maxDegree+1)*(maxDegree+1));
      else if(maxDegree>maxDegree_)
      {
        Matrix tmp = C;
        C = Matrix((maxDegree+1)*(maxDegree+1), Matrix::SYMMETRIC);
        copy(tmp, C.slice(0, 0, tmp.rows(), tmp.columns()));
      }

      if(R != R_)
      {
        Double factor = R_/R;
        UInt idx = 0;
        for(UInt n=0; n<=maxDegree; n++)
        {
          for(UInt m=0; m<2*n+1; m++)
          {
            C.row(idx)    *= factor;
            C.column(idx) *= factor;
            idx++;
          }
          factor *= R_/R;
        }
      } // (R!=R_)

      return C;
    }

    // only variances
    // --------------
    if(maxDegree < maxDegree_)
      C = C.row(0, (maxDegree+1)*(maxDegree+1));
    else if(maxDegree > maxDegree_)
    {
      Matrix tmp = C;
      C = Matrix((maxDegree+1)*(maxDegree+1), tmp.columns());
      copy(tmp, C.row(0, tmp.rows()));
    }

    if(R != R_)
    {
      Double factor = R_/R;
      UInt idx = 0;
      for(UInt n=0; n<=maxDegree; n++)
      {
        for(UInt m=0; m<2*n+1; m++)
          C.row(idx++) *= std::pow(factor, 2);
        factor *= R_/R;
      }
    } // (R!=R_)

    return C;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileTimeSplinesCovariance(const FileName &fileName,
                                    Double GM, Double R, UInt minDegree, UInt maxDegree, UInt splineDegree,
                                    const std::vector<Time> &times,
                                    const std::vector<Matrix> &sigma2)
{
  try
  {
    UInt nodeCount = times.size()+splineDegree-1;
    if(nodeCount != sigma2.size())
      throw(Exception("size of series of sigma2 does not mach to number of times."));

    OutFileArchive file(fileName, FILE_TIMESPLINESCOVARIANCE_TYPE, FILE_TIMESPLINESCOVARIANCE_VERSION);
    file<<nameValue("GM",        GM);
    file<<nameValue("R",         R);
    file<<nameValue("minDegree", minDegree);
    file<<nameValue("maxDegree", maxDegree);
    file<<nameValue("degree",    splineDegree);
    file<<nameValue("timeCount", times.size());
    for(UInt i=0; i<times.size(); i++)
      file<<nameValue("time", times.at(i));
    for(UInt i=0; i<nodeCount; i++)
      file<<nameValue("node", sigma2.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileTimeSplinesCovariance(const FileName &fileName,
                                   Double &GM, Double &R, UInt &minDegree, UInt &maxDegree, UInt &splineDegree,
                                   std::vector<Time> &times,
                                   std::vector<Matrix> &sigma2)
{
  try
  {
    InFileArchive file(fileName, FILE_TIMESPLINESCOVARIANCE_TYPE, FILE_TIMESPLINESCOVARIANCE_VERSION);
    file>>nameValue("GM",        GM);
    file>>nameValue("R",         R);
    file>>nameValue("minDegree", minDegree);
    file>>nameValue("maxDegree", maxDegree);
    file>>nameValue("degree",    splineDegree);
    UInt timeCount;
    file>>nameValue("timeCount", timeCount);
    times.resize(timeCount);
    for(UInt i=0; i<timeCount; i++)
      file>>nameValue("time", times.at(i));
    const UInt nodeCount = timeCount+splineDegree-1;
    sigma2.resize(nodeCount);
    for(UInt i=0; i<nodeCount; i++)
      file>>nameValue("node", sigma2.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
