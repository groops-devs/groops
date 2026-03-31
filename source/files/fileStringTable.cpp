/***********************************************/
/**
* @file fileStringTable.cpp
*
* @brief Read/write table of strings.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-03
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_StringList
#define DOCSTRING_FILEFORMAT_StringTable

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "fileStringTable.h"

GROOPS_REGISTER_FILEFORMAT(StringList,  "stringList")
GROOPS_REGISTER_FILEFORMAT(StringTable, "stringTable")

/***********************************************/

static Bool stripComments(std::istream &stream)
{
  try
  {
    char c;
    if(!(stream>>c))
      return FALSE;
    stream.putback(c);
    if(c != '#')
      return TRUE;
    std::string line;
    std::getline(stream, line); // skip rest of line
    return stripComments(stream);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileStringTable(const FileName &fileName, std::vector<std::vector<std::string>> &x)
{
  try
  {
    InFile file(fileName);

    while(stripComments(file))
    {
      std::string line;
      std::getline(file, line);
      std::stringstream ss(line);

      std::vector<std::string> dataLine;
      while(stripComments(ss))
      {
        std::string data;
        ss>>std::quoted(data);
        dataLine.push_back(data);
      }

      if(dataLine.size())
        x.push_back(dataLine);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/

void writeFileStringTable(const FileName &fileName, const std::vector<std::vector<std::string>> &x)
{
  try
  {
    OutFile file(fileName);
    for(const auto &line : x)
    {
      for(UInt i=0; i<line.size(); i++)
      {
        if((!line.at(i).empty()) && (line.at(i).find_first_of(" \t\n\"#") == std::string::npos))
          file<<line.at(i);              // string without special characters
        else
          file<<std::quoted(line.at(i)); // string in quotes
        if(i < line.size()-1)            // no space at the end of a line
          file<<" ";
      }
      file<<std::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/
