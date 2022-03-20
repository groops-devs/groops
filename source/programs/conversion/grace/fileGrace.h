/***********************************************/
/**
* @file fileGrace.h
*
* @brief GRACE L1A/L1B files (binary/ascii).
*
* @author Torsten Mayer-Guerr
* @date 2017-12-17
*
*/
/***********************************************/

#ifndef __GROOPS_FILEGRACE__
#define __GROOPS_FILEGRACE__

#include "base/vector3d.h"
#include "inputOutput/file.h"

/***** CLASS ***********************************/

/** @brief GRACE L1A/L1B files (binary/ascii).
* @ingroup programsConversionGroup
*/
class FileInGrace
{
  InFile file;
  Bool   isBinary;

  // change the byte order for n Bytes (e.g. Intel <-> SUN systems)
  void changeByteOrder(char *bytes, UInt n) const;

public:
  class Int8
  {
  public:
    Int value;
    operator Byte() {return value;}
  };

  class Flag
  {
  public:
    void *value;
    std::size_t size;
    Flag(void *x, std::size_t s) : value(x), size(s) {}
  };

  template<typename T> static Flag flag(T &x) {return Flag(&x, sizeof(x));}

  FileInGrace(const FileName &fileName, UInt &numberOfRecords);

  FileInGrace &operator>>(Double &x);
  FileInGrace &operator>>(Float  &x);
  FileInGrace &operator>>(Int32  &x);
  FileInGrace &operator>>(Int64  &x);
  FileInGrace &operator>>(UInt16 &x);
  FileInGrace &operator>>(UInt32 &x);
  FileInGrace &operator>>(UInt64 &x);
  FileInGrace &operator>>(Char   &x);
  FileInGrace &operator>>(Int8   &x);
  FileInGrace &operator>>(const Flag &x);
  FileInGrace &operator>>(Vector3d &vec);
};

/**********************************************************************/

#endif /* __GROOPS__ */

