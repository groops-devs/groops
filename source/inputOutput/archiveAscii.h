/***********************************************/
/**
* @file archiveAscii.h
*
* @brief Write/read archives in ASCII format.
*
* @author Torsten Mayer-Guerr
* @date 2009-11-12
*
*/
/***********************************************/

#ifndef __GROOPS_ARCHIVEASCII__
#define __GROOPS_ARCHIVEASCII__

#include "archive.h"

/** @addtogroup archiveGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Write archive files in ASCII format. */
class OutArchiveAscii : public OutArchive
{
public:
  OutArchiveAscii(std::ostream &_stream, const std::string &type, UInt version);
 ~OutArchiveAscii(){}

  OutArchiveAscii() = delete;
  OutArchiveAscii(const OutArchiveAscii &) = delete;
  OutArchiveAscii &operator=(const OutArchiveAscii &) = delete;

  ArchiveType archiveType() const override {return ASCII;}

protected:
  void endLine() override;
  void comment(const std::string &text) override;

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
  void saveDouble(Double x, Int width=26, Int precision=18, Bool science=TRUE);

  std::ostream &stream;
  Int  i_width; // Fuer die formatierte Int Ausgabe
  Bool isNewLine;
};

/***** CLASS ***********************************/

/** @brief Read archive files in ASCII format. */
class InArchiveAscii : public InArchive
{
  std::istream &stream;
  std::string  typeStr;
  UInt        _version;

  void stripComments();
  Double readDouble(std::istream &stream);

public:
  InArchiveAscii(std::istream &_stream);
 ~InArchiveAscii() {}

  ArchiveType archiveType() const override {return ASCII;}
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

#endif /* __GROOPS__ */
