/***********************************************/
/**
* @file fileArchive.cpp
*
* @brief Read and write archives files.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-01
*
*/
/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "inputOutput/logging.h"
#include "inputOutput/archive.h"
#include "inputOutput/archiveXml.h"
#include "inputOutput/archiveBinary.h"
#include "inputOutput/archiveAscii.h"
#include "inputOutput/file.h"
#include "inputOutput/fileArchive.h"

/***** FUNCTIONS *******************************/

OutFileArchive::OutFileArchive(const FileName &fileName, const std::string &type, UInt version) : archive(nullptr)
{
  open(fileName, type, version);
}

/***********************************************/

OutFileArchive::~OutFileArchive()
{
  close();
}

/***********************************************/

void OutFileArchive::open(const FileName &fileName, const std::string &type, UInt version)
{
  try
  {
    close();
    if(fileName.empty())
      return;

    // determine format from extension
    const std::string extension = String::upperCase(fileName.typeExtension());
    if(extension == "XML")
    {
      file.open(fileName);
      file.exceptions(std::ios::badbit | std::ios::failbit);
      archive = new OutArchiveXml(file, type, version);
    }
    else if(extension == "DAT")
    {
      file.open(fileName, std::ios::binary | std::ios::out);
      file.exceptions(std::ios::badbit | std::ios::failbit);
      archive = new OutArchiveBinary(file, type, version);
    }
    else
    {
      file.open(fileName);
      file.exceptions(std::ios::badbit | std::ios::failbit);
      archive = new OutArchiveAscii(file, type, version);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/

void OutFileArchive::close()
{
  if(archive)
  {
    delete archive;
    archive = nullptr;
  }
  file.close();
}

/***********************************************/

void OutFileArchive::comment(const std::string &text)
{
  if(archive)
    archive->comment(text);
}

/***********************************************/
/***********************************************/

InFileArchive::InFileArchive(const FileName &fileName, const std::string &type, UInt version) : archive(nullptr)
{
  open(fileName, type, version);
}

/***********************************************/

InFileArchive::~InFileArchive()
{
  close();
}

/***********************************************/

void InFileArchive::open(const FileName &fileName, const std::string &typeStr, UInt version)
{
  try
  {
    close();
    if(fileName.empty())
      return;

    // determine format from extension
    const std::string extension = String::upperCase(fileName.typeExtension());
    if(extension=="XML")
    {
      file.open(fileName);
      file.exceptions(std::ios::badbit | std::ios::failbit);
      char c;
      file>>c;
      file.putback(c);
      if(c!='<')
        throw(Exception("Seems not to be a valid XML file."));
      archive = new InArchiveXml(file);
    }
    else if(extension=="DAT")
    {
      file.open(fileName, std::ios::binary | std::ios::in);
      file.exceptions(std::ios::badbit | std::ios::failbit);
      char c;
      file>>c;
      file.putback(c);
      if((c!='b') && (c!='B'))
        throw(Exception("Seems not to be a groops binary file."));
      if(c=='b')
        logWarning<<"File <"<<fileName<<"> is a rather old GROOPS binary file: will not be supported in future"<<Log::endl;
      archive = new InArchiveBinary(file);
    }
    else
    {
      file.open(fileName);
      file.exceptions(std::ios::badbit | std::ios::failbit);
      archive = new InArchiveAscii(file);
    }

    // check type
    if((!typeStr.empty()) && (!type().empty()) && (type()!=typeStr))
      throw(Exception("file type is '"+type()+"' but must be '"+typeStr+"'"));

    // check version
    if(this->version() > version)
      logWarning<<"File <"<<fileName<<"> is created with a newer version ("<<this->version()<<") of GROOPS ("<<version<<"). "
                <<"This may causes problems. You should update your GROOPS software."<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/

void InFileArchive::close()
{
  if(archive)
  {
    delete archive;
    archive = nullptr;
  }
  file.close();
}

/***********************************************/

std::string InFileArchive::type() const
{
  if(!archive)
    throw(Exception("In InFileArchive::type: no file open"));
  return archive->type();
}

/***********************************************/

UInt InFileArchive::version() const
{
  if(!archive)
    throw(Exception("InFileArchive::version: no file open"));
  return archive->version();
}

/***********************************************/

Bool InFileArchive::canSeek() const
{
  if(!archive || (archive->archiveType() != InArchive::BINARY))
    return FALSE;
  return file.canSeek();
}

/***********************************************/

std::streampos InFileArchive::position()
{
  if(!archive)
    throw(Exception("InFileArchive::position: no file open"));
  return file.tellg();
}

/***********************************************/

void InFileArchive::seek(std::streampos pos)
{
  if(!archive)
    throw(Exception("InFileArchive::seek: no file open"));
  file.seekg(pos);
}

/***********************************************/
