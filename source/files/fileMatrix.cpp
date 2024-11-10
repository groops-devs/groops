/***********************************************/
/**
* @file fileMatrix.cpp
*
* @brief Read/write Matrix.
*
* @author Torsten Mayer-Guerr
* @date 2013-02-07
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_Matrix

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"

GROOPS_REGISTER_FILEFORMAT(Matrix, FILE_MATRIX_TYPE)

/***********************************************/

void writeFileMatrix(const FileName &fileName, const const_MatrixSlice &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_MATRIX_TYPE, FILE_MATRIX_VERSION);
    file<<nameValue("matrix", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileMatrix(const FileName &fileName, Matrix &x)
{
  try
  {
    InFileArchive file(fileName, ""/*arbitrary type*/, std::max(FILE_MATRIX_VERSION, FILE_INSTRUMENT_VERSION));
    if(file.type().empty() || (file.type() == FILE_MATRIX_TYPE))
      file>>nameValue("matrix", x);
    else if(file.type() == FILE_INSTRUMENT_TYPE)
      x = InstrumentFile::read(fileName).matrix();
    else
      throw(Exception("file type is '"+file.type()+"' but must be '"+FILE_MATRIX_TYPE+"' or '"+FILE_INSTRUMENT_TYPE+"'"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
