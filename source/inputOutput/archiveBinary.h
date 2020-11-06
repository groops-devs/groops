/***********************************************/
/**
* @file archiveBinary.h
*
* @brief Read/write archive files in binary format.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-11
*
*/
/***********************************************/

#ifndef __GROOPS_ARCHIVEBINARY__
#define __GROOPS_ARCHIVEBINARY__

#include "archive.h"

/** @addtogroup archiveGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Write archive files in binary format. */
class OutArchiveBinary : public OutArchive
{
public:
  OutArchiveBinary(std::ostream &_stream, const std::string &type, UInt version);
 ~OutArchiveBinary() {}

  OutArchiveBinary() = delete;
  OutArchiveBinary(const OutArchiveBinary &) = delete;
  OutArchiveBinary &operator=(const OutArchiveBinary &) = delete;

  ArchiveType archiveType() const {return BINARY;}

protected:
  void save(const std::string &x);
  void save(const Int         &x);
  void save(const UInt        &x);
  void save(const Double      &x);
  void save(const Bool        &x);
  void save(const Time        &x);
  void save(const Angle       &x);
  void save(const Vector      &x);
  void save(const const_MatrixSlice  &x);
  void save(const SphericalHarmonics &x);
  void save(const Doodson     &x);
  void save(const GnssType    &x);

private:
  std::ostream &stream;
};

/***** CLASS ***********************************/

/** @brief Read archive files in binary format. */
class InArchiveBinary : public InArchive
{
  std::istream &stream;
  std::string  typeStr;
  UInt         _version;
  Bool         oldVersion;

public:
  InArchiveBinary(std::istream &_stream);
 ~InArchiveBinary() {}

  ArchiveType archiveType() const {return BINARY;}
  std::string type()        const {return typeStr;}
  UInt        version()     const {return _version;}

protected:
  void load(std::string &x);
  void load(Int      &x);
  void load(UInt     &x);
  void load(Double   &x);
  void load(Bool     &x);
  void load(Time     &x);
  void load(Angle    &x);
  void load(Vector   &x);
  void load(Matrix   &x);
  void load(SphericalHarmonics &x);
  void load(Doodson  &x);
  void load(GnssType &x);
};

/***********************************************/

/// @}

#endif /* __GROOPS_ARCHIVEBINARY__ */
