/***********************************************/
/**
* @file fileVariationalEquation.cpp
*
* @brief Variational equation arcs.
*
* @author Torsten Mayer-Guerr
* @date 2014-03-22
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_VariationalEquation

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileSatelliteModel.h"
#include "files/fileVariationalEquation.h"

GROOPS_REGISTER_FILEFORMAT(VariationalEquation, FILE_VARIATIONALEQUATION_TYPE)

/***********************************************/

OrbitArc VariationalEquationArc::orbitArc() const
{
  try
  {
    OrbitArc orbit;
    OrbitEpoch epoch;
    for(UInt i=0; i<times.size(); i++)
    {
      epoch.time         = times.at(i);
      epoch.position     = Vector3d(pos0(3*i+0,0), pos0(3*i+1,0), pos0(3*i+2,0));
      epoch.velocity     = Vector3d(vel0(3*i+0,0), vel0(3*i+1,0), vel0(3*i+2,0));
      orbit.push_back(epoch);
    }
    return orbit;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void FileVariationalEquation::open(const FileName &fileName)
{
  try
  {
    close();
    if(!fileName.empty())
    {
      _file.open(fileName, FILE_VARIATIONALEQUATION_TYPE);
      // check for old version
      if(_file.version() < 20150524)
        throw(Exception("version of <"+fileName.str()+"> is not supported anymore. Please recompute."));
      _file>>nameValue("satellite", _satellite);
      _file>>nameValue("arcCount",  _arcCount);
      _fileName = fileName;
      _index = 0;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void FileVariationalEquation::close()
{
  if(!_fileName.empty())
    _file.close();
  _fileName = FileName();
  _arcCount = 0;
}

/***********************************************/

VariationalEquationArc FileVariationalEquation::readArc(UInt arcNo)
{
  try
  {
    if(_fileName.empty())
      return VariationalEquationArc();

    if(arcNo>=_arcCount)
      throw(Exception("arcNo >= arcCount"));

    // arc already read -> restart
    if(arcNo<_index)
      open(FileName(_fileName));

    VariationalEquationArc arc;
    while(_index <= arcNo)
    {
      _file>>nameValue("arc", arc);
      _index++;
    }
    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void writeFileVariationalEquation(const FileName &fileName, SatelliteModelPtr satellite, VariationalEquationArc arc)
{
  writeFileVariationalEquation(fileName, satellite, std::vector<VariationalEquationArc>(1, arc));
}

/***********************************************/

void writeFileVariationalEquation(const FileName &fileName, SatelliteModelPtr satellite, std::vector<VariationalEquationArc> arcList)
{
  try
  {
    OutFileArchive file(fileName, FILE_VARIATIONALEQUATION_TYPE);
    file<<nameValue("satellite", satellite);
    file<<nameValue("arcCount", arcList.size());
    for(auto iter=arcList.begin(); iter!=arcList.end(); iter++)
      file<<nameValue("arc", *iter);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void VariationalEquationArc::save(OutArchive &oa) const
{
  oa << nameValue("times",    times);
  oa << nameValue("pos0",     pos0);
  oa << nameValue("vel0",     vel0);
  oa << nameValue("PosState", PosState);
  oa << nameValue("VelState", VelState);
  oa << nameValue("rotEarth", rotEarth);
  oa << nameValue("rotSat",   rotSat);
}

/***********************************************/

void VariationalEquationArc::load(InArchive  &ia)
{
  ia >> nameValue("times",    times);
  ia >> nameValue("pos0",     pos0);
  ia >> nameValue("vel0",     vel0);
  ia >> nameValue("PosState", PosState);
  ia >> nameValue("VelState", VelState);
  ia >> nameValue("rotEarth", rotEarth);
  ia >> nameValue("rotSat",   rotSat);
}

/***********************************************/
