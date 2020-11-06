/***********************************************/
/**
* @file archiveBinary.cpp
*
* @brief Read/write archive files in binary format.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-11
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/doodson.h"
#include "base/sphericalHarmonics.h"
#include "base/gnssType.h"
#include "archive.h"
#include "archiveBinary.h"

/***********************************************/

OutArchiveBinary::OutArchiveBinary(std::ostream &_stream, const std::string &type, UInt version) : stream(_stream)
{
  stream<<"B";
  save(type);

  std::stringstream ss;
  ss<<version;
  save(ss.str());

  // padding to 64 bit
  std::string empty((1+type.size()+ss.str().size())%8, ' ');
  if(empty.size())
    stream.write(&empty.at(0), empty.size()*sizeof(char));
}

/***********************************************/

InArchiveBinary::InArchiveBinary(std::istream &_stream) : stream(_stream), _version(0)
{
  char c;
  stream>>c;
  if(c=='b')
  {
    oldVersion = TRUE;
    return; // old version
  }
  if(c!='B')
    throw(Exception("in InArchiveBinary: magic byte expected"));
  oldVersion = FALSE;

  std::string _versionStr;
  load(typeStr);
  load(_versionStr);
  _version = versionStr2version(_versionStr);

  // padding to 64 bit
  std::string empty((1+typeStr.size()+_versionStr.size())%8, ' ');
  if(empty.size())
    stream.read(&empty.at(0), empty.size()*sizeof(char));
}

/***********************************************/

void OutArchiveBinary::save(const Double &x){stream.write(reinterpret_cast<const char*>(&x), sizeof(Double));}
void InArchiveBinary::load(Double &x)       {stream.read(reinterpret_cast<char*>(&x), sizeof(Double));}

/***********************************************/

void OutArchiveBinary::save(const Angle &x) {save(static_cast<Double>(x));}
void InArchiveBinary::load(Angle &x)        {Double w; load(w); x=Angle(w);}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const Int &x)
{
  Int64 y = static_cast<Int64>(x);
  stream.write(reinterpret_cast<const char*>(&y), sizeof(Int64));
}

/***********************************************/

void InArchiveBinary::load(Int &x)
{
  if(oldVersion)
  {
    Int32 y;
    stream.read(reinterpret_cast<char*>(&y), sizeof(Int32));
    x = static_cast<Int>(y);
  }
  else
  {
    Int64 y;
    stream.read(reinterpret_cast<char*>(&y), sizeof(Int64));
    x = static_cast<Int>(y);
  }
}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const UInt &x)
{
  UInt64 y = static_cast<UInt64>(x);
  stream.write(reinterpret_cast<const char*>(&y), sizeof(UInt64));
}

/***********************************************/

void InArchiveBinary::load(UInt &x)
{
  if(oldVersion)
  {
    UInt32 y;
    stream.read(reinterpret_cast<char*>(&y), sizeof(UInt32));
    x = static_cast<UInt>(y);
  }
  else
  {
    UInt64 y;
    stream.read(reinterpret_cast<char*>(&y), sizeof(UInt64));
    x = static_cast<UInt>(y);
  }
}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const Bool &x)
{
  UInt64 y = static_cast<UInt64>(x);
  stream.write(reinterpret_cast<const char*>(&y), sizeof(UInt64));
}

/***********************************************/

void InArchiveBinary::load(Bool &x)
{
  if(oldVersion)
    stream.read(reinterpret_cast<char*>(&x), sizeof(Bool));
  else
  {
    UInt64 y;
    stream.read(reinterpret_cast<char*>(&y), sizeof(UInt64));
    x = static_cast<Bool>(y);
  }
}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const std::string &x)
{
  save(static_cast<UInt>(x.size()));
  if(x.size())
    stream.write(&x.at(0), x.size()*sizeof(char));
}

/***********************************************/

void InArchiveBinary::load(std::string &x)
{
  UInt size;
  load(size);
  x.resize(size);
  if(x.size())
    stream.read(&x.at(0), size*sizeof(char));
}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const Time &x)
{
  save(x.mjdInt());
  save(x.mjdMod());
}

/***********************************************/

void InArchiveBinary::load(Time &x)
{
  if(oldVersion)
  {
    Int mjd, milli;
    load(mjd);
    load(milli);
    x = Time(mjd,milli/(1000.*60*60*24));
  }
  else
  {
    Int    mjdInt;
    Double mjdMod;
    load(mjdInt);
    load(mjdMod);
    x = Time(mjdInt, mjdMod);
  }
}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const Doodson &x)
{
  for(UInt i=0; i<6; i++)
    save(x.d[i]);
}

/***********************************************/

void InArchiveBinary::load(Doodson &x)
{
  if(oldVersion)
  {
    throw(Exception("Old GROOPS binary file is not supported anymore."));
  }
  else
  {
    for(UInt i=0; i<6; i++)
      load(x.d[i]);
  }
}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const GnssType &x)
{
  save(x.type);
}

/***********************************************/

void InArchiveBinary::load(GnssType &x)
{
  load(x.type);
}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const Vector &x)
{
  Matrix A = x;
  save(A);
}

/***********************************************/

void InArchiveBinary::load(Vector &x)
{
  if(oldVersion)
  {
    UInt rows;
    load(rows);
    x = Vector(rows);
    if(x.rows()!=0)
      stream.read(reinterpret_cast<char*>(x.field()), x.size()*sizeof(Double));
  }
  else
  {
    Matrix A;
    load(A);
    x = A;
  }
}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const const_MatrixSlice &x)
{
  save(static_cast<UInt>(x.getType()));
  if(x.getType()==Matrix::GENERAL)
  {
    save(x.rows()); save(x.columns());
    if(x.size()==0)
      return;
    if((!x.isRowMajorOrder()) && (x.rows()==x.ld()))
      stream.write(reinterpret_cast<const char*>(x.field()), x.size()*sizeof(Double));
    else if(!x.isRowMajorOrder())
      for(UInt s=0; s<x.columns(); s++)
        stream.write(reinterpret_cast<const char*>(&x(0,s)), x.rows()*sizeof(Double));
    else
      for(UInt s=0; s<x.columns(); s++)
        for(UInt z=0; z<x.rows(); z++)
          save(x(z,s));
  }
  else
  {
    UInt uplo = (x.isUpper()) ? 0 : 1; // compability to old version
    save(uplo); save(x.rows());
    if(x.size()==0)
      return;
    if(!x.isRowMajorOrder())
    {
      if(x.isUpper())
        for(UInt s=0; s<x.columns(); s++)
          stream.write(reinterpret_cast<const char*>(&x(0,s)), (s+1)*sizeof(Double));
      else
        for(UInt s=0; s<x.columns(); s++)
          stream.write(reinterpret_cast<const char*>(&x(s,s)), (x.rows()-s)*sizeof(Double));
    }
    else
    {
      if(x.isUpper())
        for(UInt s=0; s<x.columns(); s++)
          for(UInt z=0; z<=s; z++)
            save(x(z,s));
      else
        for(UInt s=0; s<x.columns(); s++)
          for(UInt z=s; z<x.rows(); z++)
            save(x(z,s));
    }
  }
}

/***********************************************/

void InArchiveBinary::load(Matrix &x)
{
  UInt type;
  load(type);
  if(static_cast<Matrix::Type>(type)==Matrix::GENERAL)
  {
    UInt rows, columns;
    load(rows); load(columns);
    x = Matrix(rows, columns);
    if(x.size()==0)
      return;
    stream.read(reinterpret_cast<char*>(x.field()), x.size()*sizeof(Double));
  }
  else
  {
    UInt dim, uplo;
    load(uplo); load(dim);
    x = Matrix(dim, static_cast<Matrix::Type>(type), ((uplo==1) ? Matrix::LOWER : Matrix::UPPER));
    if(x.size()==0)
      return;

    if(uplo==0)
    {
      for(UInt s=0; s<x.columns(); s++)
        stream.read(reinterpret_cast<char*>(&x(0,s)), (s+1)*sizeof(Double));
    }
    else
    {
      for(UInt s=0; s<x.columns(); s++)
        stream.read(reinterpret_cast<char*>(&x(s,s)), (x.rows()-s)*sizeof(Double));
    }
  }
}

/***********************************************/
/***********************************************/

void OutArchiveBinary::save(const SphericalHarmonics &harm)
{
  save(harm.GM());
  save(harm.R());
  save(harm.cnm());
  save(harm.snm());
  save(harm.sigma2cnm());
  save(harm.sigma2snm());
}

/***********************************************/

void InArchiveBinary::load(SphericalHarmonics &harm)
{
  Double GM, R;
  Matrix cnm, snm, sigma2cnm, sigma2snm;

  load(GM);
  load(R);
  load(cnm);
  load(snm);
  load(sigma2cnm);
  load(sigma2snm);

  harm = SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm);
}

/***********************************************/
