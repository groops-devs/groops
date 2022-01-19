/***********************************************/
/**
* @file file.cpp
*
* read and write files.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-06
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/constants.h"
#include "base/string.h"
#include "external/compress.h"
#include "file.h"

/***** CLASS ***********************************/

#ifdef GROOPS_DISABLE_Z
#else

namespace zlib
{
#include <zlib.h>
}

// This is a reimplementation of gzstream (Deepak Bandyopadhyay, Lutz Kettner)
class StreambufGZ : public std::streambuf
{
private:
  static constexpr int bufferSize = 47+256;    // size of data buff
  // totals 512 bytes under g++ for igzstream at the end.

  zlib::gzFile       file;               // file handle for compressed file
  char               buffer[bufferSize]; // data buffer
  Bool               opened;             // open/close state of stream
  std::ios::openmode mode;               // I/O mode

  StreambufGZ::int_type flush_buffer();

public:
  StreambufGZ();
 ~StreambufGZ() {close();}

  bool is_open() const {return opened;}

  StreambufGZ *open(const FileName &fileName, std::ios::openmode openMode);
  StreambufGZ *close();

  virtual StreambufGZ::int_type underflow() override;
  virtual StreambufGZ::int_type overflow(StreambufGZ::int_type c) override;
  virtual StreambufGZ::int_type sync() override;
};

/***********************************************/

// Separate the writing of the buffer from overflow() and sync() operation.
StreambufGZ::int_type StreambufGZ::flush_buffer()
{
  auto w = pptr() - pbase();
  if(zlib::gzwrite(file, pbase(), w) != w)
    return traits_type::eof();
  pbump(-w);
  return w;
}

/***********************************************/

StreambufGZ::StreambufGZ() : opened(FALSE)
{
  setp(buffer, buffer+(bufferSize-1));
  setg(buffer+4, buffer+4, buffer+4); // beginning of putback area, read position, end position
  // ASSERT: both input & output capabilities will not be used together
}

/***********************************************/

StreambufGZ *StreambufGZ::open(const FileName &fileName, std::ios::openmode openMode)
{
  if(is_open())
    return nullptr;
  mode = openMode;
  // no append nor read/write mode
  if((mode & std::ios::ate) || (mode & std::ios::app) || ((mode & std::ios::in) && (mode & std::ios::out)))
    throw(Exception("openMode combination not allowed for .gz files"));
  file = zlib::gzopen(fileName.c_str(), (mode & std::ios::out) ? "wb" : "rb");
  if(!file)
    return nullptr;
  opened = TRUE;
  return this;
}

/***********************************************/

StreambufGZ *StreambufGZ::close()
{
  if(is_open())
  {
    sync();
    opened = FALSE;
    if(zlib::gzclose(file) == Z_OK)
      return this;
  }
  return nullptr;
}

/***********************************************/

// used for input buffer only
StreambufGZ::int_type StreambufGZ::underflow()
{
  if(gptr() && (gptr() < egptr()))
    return traits_type::to_int_type(*gptr());

  if(!(mode & std::ios::in) || !opened)
    return traits_type::eof();

  // Josuttis' implementation of inbuf
  auto n_putback = gptr() - eback();
  if(n_putback > 4)
    n_putback = 4;
  memcpy(buffer+(4-n_putback), gptr()-n_putback, n_putback);

  auto num = zlib::gzread(file, buffer+4, bufferSize-4);
  if(num <= 0) // ERROR or EOF
    return traits_type::eof();

  // reset buffer pointers
  setg(buffer+(4-n_putback), buffer+4, buffer+4+num); // beginning of putback area, read position, end of buffer

  return traits_type::to_int_type(*gptr()); // return next character
}

/***********************************************/

// used for output buffer only
StreambufGZ::int_type StreambufGZ::overflow(StreambufGZ::int_type c)
{
  if(!(mode & std::ios::out) || !opened)
    return traits_type::eof();
  if(c != traits_type::eof())
  {
    *pptr() = c;
    pbump(1);
  }
  if(flush_buffer() == traits_type::eof())
    return traits_type::eof();
  return c;
}

/***********************************************/

StreambufGZ::int_type StreambufGZ::sync()
{
  if(pptr() && (pptr() > pbase()) && (flush_buffer() == traits_type::eof()))
    return -1;
  return 0;
}

#endif // LIB_Z

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

StreamBase::StreamBase() : buffer(nullptr), canSeek_(FALSE) {}
StreamBase::~StreamBase() {close();}

/***********************************************/

void StreamBase::open(const FileName &fileName, std::ios::openmode openMode)
{
  try
  {
    close();
    if(fileName.empty())
      return;
    this->fileName_ = fileName;
    this->canSeek_  = TRUE;

    // determine format from extension
    std::string fileFormat = String::upperCase(fileName.packExtension());

    // determine format from magic bytes if possible
    if(openMode == std::ios::in)
    {
      std::ifstream file(fileName.c_str(), std::ios::binary);
      unsigned char magic[2] = {0};
      file.read(reinterpret_cast<char*>(magic), sizeof(magic));
      const unsigned char magicCompress[2] = {0x1f, 0x9d};
      const unsigned char magicZlib[2]     = {0x1f, 0x8b};
      if(std::memcmp(magic, magicZlib, sizeof(magic)) == 0)
        fileFormat = "GZ";
      else if(std::memcmp(magic, magicCompress, sizeof(magic)) == 0)
        fileFormat = "Z";
    }

    if(fileFormat == "GZ")
    {
#ifdef GROOPS_DISABLE_Z
      throw(Exception("compiled without Z library"));
#else
      buffer = new StreambufGZ();
      std::ios::init(buffer);
      if(!static_cast<StreambufGZ*>(buffer)->open(fileName, openMode))
        clear(rdstate() | std::ios::badbit);
      canSeek_ = FALSE;
#endif
    }
    else if(fileFormat == "Z")
    {
      if(openMode != std::ios::in)
        throw(Exception("old .Z compress implemented for input only"));

      buffer = new std::stringbuf(decompress_file(fileName.c_str()));
      std::ios::init(buffer);
    }
    else
    {
      buffer = new std::filebuf();
      std::ios::init(buffer);
      if(!static_cast<std::filebuf*>(buffer)->open(fileName.str(), openMode))
        clear(rdstate() | std::ios::badbit);
    }

    if(!good())
      throw(Exception("error by opening file"));
    exceptions(std::ios::badbit);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/

void StreamBase::close()
{
  try
  {
    if(buffer)
    {
      delete buffer;
      buffer = nullptr;
      std::ios::init(nullptr);
    }
    fileName_ = FileName();
    canSeek_  = FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName_.str()+">", e)
  }
}

/***********************************************/
