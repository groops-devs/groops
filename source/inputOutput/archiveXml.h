/***********************************************/
/**
* @file archiveXml.h
*
* @brief Read/write archive files in XML format.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-01
*
*/
/***********************************************/

#ifndef __GROOPS_ARCHIVEXML__
#define __GROOPS_ARCHIVEXML__

#include "archive.h"
#include "parser/xml.h"

/** @addtogroup archiveGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Write archive files in XML format. */
class OutArchiveXml : public OutArchive
{
public:
  OutArchiveXml(std::ostream &_stream, const std::string &type, UInt version);
 ~OutArchiveXml();

  OutArchiveXml() = delete;
  OutArchiveXml(const OutArchiveXml &) = delete;
  OutArchiveXml &operator=(const OutArchiveXml &) = delete;

  ArchiveType archiveType() const {return XML;}

protected:
  void startTag(const std::string &name);
  void endTag  (const std::string &name);

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
  std::ostream           &stream;
  std::stack<XmlNodePtr> stack;
};

/***** CLASS ***********************************/

/** @brief Read archive files in XML format. */
class InArchiveXml : public InArchive
{
  std::stack<XmlNodePtr> stack;
  std::string typeStr;
  UInt       _version;

public:
  InArchiveXml(std::istream &stream);
 ~InArchiveXml();

  ArchiveType archiveType() const {return XML;}
  std::string type()        const {return typeStr;}
  UInt        version()     const {return _version;}

protected:
  void startTag(const std::string &name);
  void endTag  (const std::string &name);

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

#endif /* __GROOPS_ARCHIVEXML__ */
