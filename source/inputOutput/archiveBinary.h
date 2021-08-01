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

  ArchiveType archiveType() const override {return BINARY;}

protected:
  void save(const std::string &x) override;
  void save(const Int         &x) override;
  void save(const UInt        &x) override;
  void save(const Double      &x) override;
  void save(const Bool        &x) override;
  void save(const Time        &x) override;
  void save(const Angle       &x) override;
  void save(const Vector      &x) override;
  void save(const const_MatrixSlice  &x) override;
  void save(const SphericalHarmonics &x) override;
  void save(const Doodson     &x) override;
  void save(const GnssType    &x) override;

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

  ArchiveType archiveType() const override {return BINARY;}
  std::string type()        const override {return typeStr;}
  UInt        version()     const override {return _version;}

protected:
  void load(std::string &x) override;
  void load(Int      &x) override;
  void load(UInt     &x) override;
  void load(Double   &x) override;
  void load(Bool     &x) override;
  void load(Time     &x) override;
  void load(Angle    &x) override;
  void load(Vector   &x) override;
  void load(Matrix   &x) override;
  void load(SphericalHarmonics &x) override;
  void load(Doodson  &x) override;
  void load(GnssType &x) override;
};

/***********************************************/

/// @}

#endif /* __GROOPS_ARCHIVEBINARY__ */
