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
    this->_maxDegree = maxDegree;
    this->_minDegree = minDegree;

    close();
    if(!name.empty())
    {
      file.open(name, FILE_TIMESPLINESGRAVITYFIELD_TYPE);
      file>>nameValue("GM", GM);
      file>>nameValue("R",  R);
      file>>nameValue("degree",  degree);
      file>>nameValue("timeCount", count);
      times.resize(count);
      for(UInt i=0; i<times.size(); i++)
        file>>nameValue("time", times.at(i));

      // If we have a binary file, we can efficiently seek to the needed position in the file.
      if(file.canSeek())
        dataStart = file.position();
      harmonics.resize(degree+1);
      for(UInt i=0; i<harmonics.size(); i++)
      {
        Matrix cnm, snm;
        file>>beginGroup("node");
        file>>nameValue("cnm", cnm);
        file>>nameValue("snm", snm);
        harmonics.at(i) = SphericalHarmonics(GM, R, cnm, snm).get(maxDegree, minDegree);
        file>>endGroup("node");
        this->_maxDegree = harmonics.at(i).maxDegree();
        if(file.canSeek() && (i==0))
          dataEpochSize = file.position() - dataStart;
      }

      index = 0;
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
}

/***********************************************/

SphericalHarmonics InFileTimeSplinesGravityfield::sphericalHarmonics(const Time &time, Double factor)
{
  try
  {
    if(file.fileName().empty())
      return SphericalHarmonics();

    if((time<times.at(0)) || (time>=times.at(times.size()-1)))
    {
      logWarning<<"No TimeSplinesGravityfield data for time "<<time.dateTimeStr()<<Log::endl;
      return SphericalHarmonics();
    }

    // Schon ueberlesen? -> dann von vorne anfangen
    if(time<times.at(index))
      open(FileName(file.fileName()), this->_maxDegree, this->_minDegree);

    // Seek to appropriate interval. The current position in the file is degree+1 intervals ahead of index
    if(file.canSeek() && (time>times.at(index+1))) // If the time is larger than the last index that has already been read in
    {
      UInt seekIndex = index;
      while(time>=times.at(seekIndex+1))
        seekIndex++;
      index = seekIndex - (degree+1);
      file.seek(dataStart + static_cast<std::streamoff>(seekIndex * dataEpochSize));
    }

    // Interpolationsintervall suchen
    while(time>=times.at(index+1))
    {
      for(UInt i=0; i+1<harmonics.size(); i++)
        harmonics.at(i) = harmonics.at(i+1);
      Matrix cnm, snm;
      file>>beginGroup("node");
      file>>nameValue("cnm", cnm);
      file>>nameValue("snm", snm);
      file>>endGroup("node");
      harmonics.at(degree) = SphericalHarmonics(GM, R, cnm, snm).get(_maxDegree, _minDegree);
      index++;
    }

    const Double       t     = (time-times.at(index)).mjd()/(times.at(index+1)-times.at(index)).mjd();
    const Vector       coeff = BasisSplines::compute(t, degree);
    SphericalHarmonics harm  = factor*coeff(0) * harmonics.at(0);
    for(UInt i=1; i<coeff.rows(); i++)
      harm += factor * coeff(i) * harmonics.at(i);
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

    OutFileArchive file(fileName, FILE_TIMESPLINESGRAVITYFIELD_TYPE);
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
    InFileArchive file(fileName, FILE_TIMESPLINESGRAVITYFIELD_TYPE);
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
      file.open(name, FILE_TIMESPLINESCOVARIANCE_TYPE);
      file>>nameValue("GM",        _GM);
      file>>nameValue("R",         _R);
      file>>nameValue("minDegree", _minDegree);
      file>>nameValue("maxDegree", _maxDegree);
      file>>nameValue("degree",    degree);

      mustCut    = (maxDegree < _maxDegree);
      _maxDegree = std::min(maxDegree, _maxDegree);
      _minDegree = std::max(minDegree, _minDegree);

      file>>nameValue("timeCount", count);
      times.resize(count);
      for(UInt i=0; i<times.size(); i++)
        file>>nameValue("time", times.at(i));

      cov.resize(degree+1);
      for(UInt i=0; i<cov.size(); i++)
      {
        file>>nameValue("node", cov.at(i));
        if(mustCut)
        {
          if(cov.at(i).getType() == Matrix::SYMMETRIC) // full covariance matrix?
            cov.at(i) = cov.at(i).slice(0,0,(_maxDegree+1)*(_maxDegree+1),(_maxDegree+1)*(_maxDegree+1));
          else
            cov.at(i) = cov.at(i).row(0,(_maxDegree+1)*(_maxDegree+1));
        }
      }

      index = 0;
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
}

/***********************************************/

Matrix InFileTimeSplinesCovariance::covariance(const Time &time, Double factor, UInt maxDegree, UInt minDegree, Double GM, Double R)
{
  try
  {
    if(file.fileName().empty())
      return Matrix();

    if((time<times.at(0)) || (time>=times.at(times.size()-1)))
    {
      logWarning<<"No TimeSplinesCovariance data for time "<<time.dateTimeStr()<<Log::endl;
      return Matrix();
    }

    // Schon ueberlesen? -> dann von vorne anfangen
    if(time<times.at(index))
      open(FileName(file.fileName()), this->_maxDegree, this->_minDegree);

    // Interpolationsintervall suchen
    while(time>=times.at(index+1))
    {
      for(UInt i=0; i+1<cov.size(); i++)
        cov.at(i) = cov.at(i+1);
      file>>nameValue("node", cov.at(degree));
      if(mustCut)
      {
        if(cov.at(degree).getType() == Matrix::SYMMETRIC) // full covariance matrix?
          cov.at(degree) = cov.at(degree).slice(0,0,(_maxDegree+1)*(_maxDegree+1),(_maxDegree+1)*(_maxDegree+1));
        else
          cov.at(degree) = cov.at(degree).row(0,(_maxDegree+1)*(_maxDegree+1));
      }
      index++;
    }

    const Double t     = (time-times.at(index)).mjd()/(times.at(index+1)-times.at(index)).mjd();
    const Vector coeff = BasisSplines::compute(t, degree);

    Matrix C = coeff(0) * cov.at(0);
    for(UInt i=1; i<coeff.rows(); i++)
      C += coeff(i) * cov.at(i);

    if(GM<=0) GM = _GM;
    if(R<=0)  R  = _R;
    if(maxDegree==INFINITYDEGREE) maxDegree = _maxDegree;
    minDegree = std::min(minDegree, maxDegree);

    if(_GM/GM*factor != 1.0)
      C *= pow(_GM/GM*factor, 2);

    // minDegree
    // ---------
    if(C.getType() == Matrix::SYMMETRIC) // full covariance matrix?
    {
      for(UInt idx=0; idx<pow(std::max(_minDegree, minDegree),2); idx++)
      {
        C.row(idx)    *= 0;
        C.column(idx) *= 0;
      }
    }
    else
      for(UInt idx=0; idx<pow(std::max(_minDegree, minDegree),2); idx++)
        C.row(idx) *= 0;

    // Quick return possible?
    // ----------------------
    if((R==_R)&&(maxDegree==_maxDegree))
      return C;

    // full covariance matrix?
    // -----------------------
    if(C.getType() == Matrix::SYMMETRIC)
    {
      if(maxDegree<_maxDegree)
        C = C.slice(0, 0, (maxDegree+1)*(maxDegree+1), (maxDegree+1)*(maxDegree+1));
      else if(maxDegree>_maxDegree)
      {
        Matrix tmp = C;
        C = Matrix((maxDegree+1)*(maxDegree+1), Matrix::SYMMETRIC);
        copy(tmp, C.slice(0, 0, tmp.rows(), tmp.columns()));
      }

      if(R!=_R)
      {
        Double factor = _R/R;
        UInt idx = 0;
        for(UInt n=0; n<=maxDegree; n++)
        {
          for(UInt m=0; m<2*n+1; m++)
          {
            C.row(idx)    *= factor;
            C.column(idx) *= factor;
            idx++;
          }
          factor *= _R/R;
        }
      } // (R!=_R)

      return C;
    }

    // only variances
    // --------------
    if(maxDegree<_maxDegree)
      C = C.row(0, (maxDegree+1)*(maxDegree+1));
    else if(maxDegree>_maxDegree)
    {
      Matrix tmp = C;
      C = Matrix((maxDegree+1)*(maxDegree+1), tmp.columns());
      copy(tmp, C.row(0, tmp.rows()));
    }

    if(R!=_R)
    {
      Double factor = _R/R;
      UInt idx = 0;
      for(UInt n=0; n<=maxDegree; n++)
      {
        for(UInt m=0; m<2*n+1; m++)
          C.row(idx++) *= pow(factor,2);
        factor *= _R/R;
      }
    } // (R!=_R)

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

    OutFileArchive file(fileName, FILE_TIMESPLINESCOVARIANCE_TYPE);
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
    InFileArchive file(fileName, FILE_TIMESPLINESCOVARIANCE_TYPE);
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
