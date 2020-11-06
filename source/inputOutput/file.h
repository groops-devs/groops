/***********************************************/
/**
* @file file.h
*
* read and write files.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-06
*
*/
/***********************************************/

#ifndef __GROOPS_FILE__
#define __GROOPS_FILE__

#include "base/exception.h"
#include "inputOutput/fileName.h"

/** @addtogroup inputOutputGroup */
/// @{

/***** CLASS ***********************************/

// Internal class
class StreamBase : virtual public std::ios
{
protected:
  std::streambuf *buffer;
  FileName        fileName_;
  Bool            canSeek_;

public:
  StreamBase();
 ~StreamBase();

  StreamBase(const StreamBase &) = delete;
  StreamBase &operator=(const StreamBase &) = delete;

  void open(const FileName &fileName, std::ios::openmode openMode);
  void close();

  FileName fileName() const {return fileName_;}
  Bool     canSeek()  const {return canSeek_;}
};

/***** CLASS ***********************************/

class OutFile : public StreamBase, public std::ostream
{
public:
  OutFile() : std::ostream(nullptr) {}
  OutFile(const FileName &fileName, std::ios::openmode openMode=std::ios::out) : std::ostream(nullptr) {open(fileName, openMode);}
 ~OutFile() {}

  OutFile(const OutFile &) = delete;
  OutFile &operator=(const OutFile &) = delete;

  void open(const FileName &fileName, std::ios::openmode openMode=std::ios::out) {StreamBase::open(fileName, openMode);}
  void close() {StreamBase::close();}

  OutFile &operator<<(std::ostream  &(*pf)(std::ostream  &)) {pf(*this); return *this;}
  OutFile &operator<<(std::ios      &(*pf)(std::ios      &)) {pf(*this); return *this;}
  OutFile &operator<<(std::ios_base &(*pf)(std::ios_base &)) {pf(*this); return *this;}
  template<typename T> inline OutFile &operator<<(const T &x);
};

/***** CLASS ***********************************/

class InFile : public StreamBase, public std::istream
{
public:
  InFile() : std::istream(nullptr) {}
  InFile(const FileName &fileName, std::ios::openmode openMode=std::ios::in) : std::istream(nullptr) {open(fileName, openMode);}
 ~InFile() {}

  InFile(const InFile &) = delete;
  InFile &operator=(const InFile &) = delete;

  void open(const FileName &fileName, std::ios::openmode openMode=std::ios::in) {StreamBase::open(fileName, openMode);}
  void close() {StreamBase::close();}

  InFile &operator>>(std::istream  &(*pf)(std::istream  &)) {pf(*this); return *this;}
  InFile &operator>>(std::ios      &(*pf)(std::ios      &)) {pf(*this); return *this;}
  InFile &operator>>(std::ios_base &(*pf)(std::ios_base &)) {pf(*this); return *this;}
  template<typename T> inline InFile &operator>>(T &x);
};

/// @}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

template<typename T> inline OutFile &OutFile::operator<<(const T &x)
{
  try
  {
    if(fileName_.empty())
      throw(Exception("no file open"));
    static_cast<std::ostream&>(*this)<<x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName_.str()+">", e)
  }
}

/***********************************************/

template<typename T>
inline InFile &InFile::operator>>(T &x)
{
  try
  {
    if(fileName_.empty())
      throw(Exception("no file open"));
    static_cast<std::istream&>(*this)>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName_.str()+">", e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
