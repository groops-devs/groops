/***********************************************/
/**
* @file fileGriddedDataTimeSeries.cpp
*
* @brief Read/write time series of gridded values.
*
* @author Torsten Mayer-Guerr
* @date 2020-07-28
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_GriddedDataTimeSeries

#include "base/import.h"
#include "base/basisSplines.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileGriddedData.h"
#include "files/fileGriddedDataTimeSeries.h"

GROOPS_REGISTER_FILEFORMAT(GriddedDataTimeSeries, FILE_GRIDDEDDATATIMESERIES_TYPE)

/***********************************************/

void InFileGriddedDataTimeSeries::open(const FileName &name)
{
  try
  {
    close();
    if(!name.empty())
    {
      UInt epochCount;

      file.open(name, FILE_GRIDDEDDATATIMESERIES_TYPE);
      file>>nameValue("splineDegree", splineDegree_);
      file>>nameValue("timeCount",    epochCount);
      file>>nameValue("dataCount",    dataCount_);
      file>>nameValue("grid",         grid_);
      times_.resize(epochCount);
      for(UInt i=0; i<times_.size(); i++)
        file>>nameValue("time", times_.at(i));

      // If we have a binary file, we can efficiently seek to the needed position in the file.
      if(file.canSeek())
        seekStart = file.position();
      seekSize = 0;

      indexFile = 0;
      indexData = nodeCount();
      data_.resize(splineDegree_+1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

void InFileGriddedDataTimeSeries::close()
{
  file.close();
  grid_ = GriddedData();
  times_.clear();
  dataCount_    = 0;
  splineDegree_ = 0;
}

/***********************************************/

Matrix InFileGriddedDataTimeSeries::data(UInt idNode) // points x data columns
{
  try
  {
    if(file.fileName().empty())
      return Matrix();

    // must restart?
    if(indexFile > idNode)
    {
      if(file.canSeek() && seekSize)
        indexFile = 0;
      else
        open(file.fileName());
    }

    Matrix data(grid_.points.size(), dataCount_);
    while(indexFile <= idNode)
    {
      // Seek to appropriate interval.
      if(file.canSeek() && seekSize)
      {
        file.seek(seekStart + static_cast<std::streamoff>(idNode * seekSize));
        indexFile = idNode;
      }

      file>>beginGroup("node");
      for(UInt i=0; i<grid_.points.size(); i++)
      {
        file>>beginGroup("points");
        for(UInt k=0; k<dataCount_; k++)
          file>>nameValue("value", data(i, k));
        file>>endGroup("points");
      }
      file>>endGroup("node");

      if(file.canSeek() && (indexFile == 0))
        seekSize = file.position() - seekStart;
      indexFile++;
    }

    return data;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix InFileGriddedDataTimeSeries::data(const Time &time)
{
  try
  {
    if((time < times_.front()) || (time > times_.back()))
      throw(Exception(time.dateTimeStr()+" outside interval ["+times_.front()%"%D_%T, "s+times_.back()%"%D_%T] of <"s+file.fileName().str()+">"));

    // find time interval and read missing nodes
    const UInt idx   = std::min(std::distance(times_.begin(), std::upper_bound(times_.begin(), times_.end(), time)),
                                static_cast<std::vector<Time>::difference_type>(times_.size()-1))-1;
    const UInt start = (idx > indexData) ? std::max(idx, indexData+data_.size()) : idx;
    const UInt end   = (idx > indexData) ? idx+data_.size() : std::min(idx+data_.size(), indexData);
    for(UInt i=start; i<end; i++)
      data_.at(i%data_.size()) = data(i);
    indexData = idx;

    // spline interploation in interval
    const Double t     = (time-times_.at(idx)).mjd()/(times_.at(idx+1)-times_.at(idx)).mjd();
    const Vector coeff = BasisSplines::compute(t, splineDegree_);
    Matrix data = coeff(0) * data_.at(indexData%data_.size());
    for(UInt i=1; i<coeff.rows(); i++)
      axpy(coeff(i), data_.at((indexData+i)%data_.size()), data);
    return data;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void writeFileGriddedDataTimeSeries(const FileName &fileName, UInt splineDegree, const std::vector<Time> &times,
                                    const GriddedData &grid, const std::vector<Matrix> &data)
{
  try
  {
    if(!grid.isValid())
      throw(Exception("GriddedData is not valid"));
    if(times.size()+splineDegree != data.size()+1)
      throw(Exception("spline degree, times.size(), and data.size() not fit"));
    const UInt dataCount = (data.size() ? data.front().columns() : 0);

    OutFileArchive file(fileName, FILE_GRIDDEDDATATIMESERIES_TYPE);
    file<<nameValue("splineDegree", splineDegree);
    file<<nameValue("timeCount",    times.size());
    file<<nameValue("dataCount",    dataCount);
    file<<nameValue("grid",         grid);
    file.comment("times");
    file.comment("=====");
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      file<<nameValue("time", times.at(idEpoch));
    for(UInt idEpoch=0; idEpoch<data.size(); idEpoch++)
    {
      file<<beginGroup("node");

      // comment
      std::string str;
      for(UInt i=0; i<dataCount; i++)
      {
        std::string str2 = "data"+i%"%i"s;
        str += str2 + std::string(26-str2.size(), ' ');
      }
      file.comment(times.at(idEpoch).dateTimeStr());
      file.comment(str);
      file.comment(std::string(str.size(), '='));

      for(UInt i=0; i<grid.points.size(); i++)
      {
        file<<beginGroup("points");
        for(UInt k=0; k<dataCount; k++)
          file<<nameValue("value", data.at(idEpoch)(i, k));
        file<<endGroup("points");
      }
      file<<endGroup("node");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/
