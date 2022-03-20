/***********************************************/
/**
* @file fileGrace.cpp
*
* @brief GRACE L1A/L1B files (isBinary/ascii).
*
* @author Torsten Mayer-Guerr
* @date 2017-12-17
*
*/
/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "inputOutput/logging.h"
#include "inputOutput/file.h"
#include "fileGrace.h"

/***********************************************/

FileInGrace::FileInGrace(const FileName &fileName, UInt &numberOfRecords)
{
  try
  {
    numberOfRecords = 0;

    try
    {
      file.open(fileName, std::ios::in|std::ios::binary);
      file.exceptions(std::ios::badbit | std::ios::failbit);
    }
    catch(std::exception &/*e*/)
    {
      logError<<"error by opening file <"<<fileName<<">: continue..."<<Log::endl;
      return;
    }

    std::string line;
    while(std::getline(file, line))
    {
      // GRACE
      if(line.find("NUMBER OF DATA RECORDS") == 0)
        numberOfRecords = String::toInt(line.substr(31, 10));
      if(line.find("FILE FORMAT") == 0)
        isBinary = (String::toInt(line.substr(31, 2)) == 0);
      if(line.find("END OF HEADER") == 0)
        break;
      // GRACE-FO
      if(line.find("    num_records") == 0)
      {
        numberOfRecords = String::toInt(line.substr(16, 10));
        isBinary = FALSE;
      }
      if(line.find("# End of YAML header") == 0)
        break;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void FileInGrace::changeByteOrder(char *bytes, UInt n) const
{
  for(UInt i=0; i<n/2; i++)
    std::swap(bytes[i], bytes[n-i-1]);
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(Double &x)
{
  try
  {
    if(isBinary)
    {
      file.read(reinterpret_cast<char*>(&x), sizeof(x));
      changeByteOrder(reinterpret_cast<char*>(&x), sizeof(x));
    }
    else
      file>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(Float &x)
{
  try
  {
    if(isBinary)
    {
      file.read(reinterpret_cast<char*>(&x), sizeof(x));
      changeByteOrder(reinterpret_cast<char*>(&x), sizeof(x));
    }
    else
      file>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(Int32 &x)
{
  try
  {
    if(isBinary)
    {
      file.read(reinterpret_cast<char*>(&x), sizeof(x));
      changeByteOrder(reinterpret_cast<char*>(&x), sizeof(x));
    }
    else
      file>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(Int64 &x)
{
  try
  {
    if(isBinary)
    {
      file.read(reinterpret_cast<char*>(&x), sizeof(x));
      changeByteOrder(reinterpret_cast<char*>(&x), sizeof(x));
    }
    else
      file>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(UInt16 &x)
{
  try
  {
    if(isBinary)
    {
      file.read(reinterpret_cast<char*>(&x), sizeof(x));
      changeByteOrder(reinterpret_cast<char*>(&x), sizeof(x));
    }
    else
      file>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/***********************************************/

FileInGrace &FileInGrace::operator>>(UInt32 &x)
{
  try
  {
    if(isBinary)
    {
      file.read(reinterpret_cast<char*>(&x), sizeof(x));
      changeByteOrder(reinterpret_cast<char*>(&x), sizeof(x));
    }
    else
      file>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(UInt64 &x)
{
  try
  {
    if(isBinary)
    {
      file.read(reinterpret_cast<char*>(&x), sizeof(x));
      changeByteOrder(reinterpret_cast<char*>(&x), sizeof(x));
    }
    else
      file>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(Char &x)
{
  try
  {
    if(isBinary)
    {
      file.read(reinterpret_cast<char*>(&x), sizeof(x));
      changeByteOrder(reinterpret_cast<char*>(&x), sizeof(x));
    }
    else
      file>>x;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(Int8 &x)
{
  try
  {
    if(isBinary)
    {
      Byte b;
      file.read(&b, 1);
      x.value = b;
    }
    else
      file>>x.value;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(const Flag &x)
{
  try
  {
    if(isBinary)
    {
      file.read(reinterpret_cast<char*>(x.value), x.size);
      changeByteOrder(reinterpret_cast<char*>(x.value), x.size);
    }
    else
    {
      UInt64 v = 0;
      Byte   c;
      for(UInt i=0; i<8*x.size; i++)
      {
        file>>c;
        v = (v<<1) + c-'0';
      }
      switch(x.size)
      {
        case 1 : *reinterpret_cast<Byte  *>(x.value) = v; break;
        case 2 : *reinterpret_cast<UInt16*>(x.value) = v; break;
        case 4 : *reinterpret_cast<UInt32*>(x.value) = v; break;
        case 8 : *reinterpret_cast<UInt64*>(x.value) = v; break;
        default:
          throw(Exception("something strange"));
      };
    }
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

FileInGrace &FileInGrace::operator>>(Vector3d &vec)
{
  try
  {
    *this>>vec.x()>>vec.y()>>vec.z();
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
