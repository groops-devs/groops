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

  ArchiveType archiveType() const {return ASCII;}

protected:
  void endLine();
  void comment(const std::string &text);

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

  ArchiveType archiveType() const {return ASCII;}
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

#endif /* __GROOPS__ */
