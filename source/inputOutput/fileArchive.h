/***********************************************/
/**
* @file fileArchive.h
*
* read and write archives files.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-01
*
*/
/***********************************************/

#ifndef __GROOPS_FILEARCHIVE__
#define __GROOPS_FILEARCHIVE__

#include "base/exception.h"
#include "inputOutput/archive.h"
#include "inputOutput/file.h"

/** @addtogroup archiveGroup */
/// @{

/***** CONSTANTS ********************************/

const UInt FILE_VERSION = 20200123;    // date of last change (ArcList, InstrumentFile restructured)
// const UInt FILE_VERSION = 20190429; // date of last change (SatelliteModel Surface hasThermalReemission)
// const UInt FILE_VERSION = 20190304; // date of last change (GnssStationInfo)
// const UInt FILE_VERSION = 20170920; // date of last change
// const UInt FILE_VERSION = 20150524; // date of last change

/***** CLASS ***********************************/

class OutFileArchive
{
  OutFile     file;
  OutArchive *archive;

public:
  OutFileArchive() : archive(nullptr)  {}
  OutFileArchive(const FileName &fileName, const std::string &type);
 ~OutFileArchive();

  OutFileArchive(const InFile &) = delete;
  OutFileArchive &operator=(const InFile &) = delete;

  void open(const FileName &fileName, const std::string &type);
  void close();

  FileName fileName() const {return file.fileName();}

  OutArchive &outArchive() {return *archive;}
  void comment(const std::string &text);

  template<typename T> inline OutFileArchive &operator<<(const T &x);
};

/***** CLASS ***********************************/

class InFileArchive
{
  InFile     file;
  InArchive *archive;

public:
  InFileArchive() : archive(nullptr) {}
  InFileArchive(const FileName &fileName, const std::string &type);
 ~InFileArchive();

  InFileArchive(const InFile &) = delete;
  InFileArchive &operator=(const InFile &) = delete;

  void open(const FileName &fileName, const std::string &type);
  void close();

  FileName    fileName() const {return file.fileName();}
  std::string type()     const;
  UInt        version()  const;

  Bool           canSeek() const;
  std::streampos position();
  void           seek(std::streampos pos);

  template<typename T> inline InFileArchive &operator>>(T &x);
  template<typename T> inline InFileArchive &operator>>(const T &x);
};

/// @}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

template<typename T>
inline OutFileArchive &OutFileArchive::operator<<(const T &x)
{
  try
  {
    if(!archive)
      throw(Exception("no file open"));
    (*archive)<<x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+file.fileName().str()+">", e)
  }
}

/***********************************************/

template<typename T>
inline InFileArchive &InFileArchive::operator>>(T &x)
{
  try
  {
    if(!archive)
      throw(Exception("no file open"));
    (*archive)>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+file.fileName().str()+">", e)
  }
}

/***********************************************/

template<typename T>
inline InFileArchive &InFileArchive::operator>>(const T &x)
{
  try
  {
    if(!archive)
      throw(Exception("no file open"));
    (*archive)>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+file.fileName().str()+">", e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
