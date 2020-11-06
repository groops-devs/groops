/***********************************************/
/**
* @file fileFormatRegister.h
*
* @brief Desciption of file format in groops.
*
* @author Torsten Mayer-Guerr
* @date 2020-06-15
*
*/
/***********************************************/

#ifndef __GROOPS_FILEFORMATGREGISTER__
#define __GROOPS_FILEFORMATGREGISTER__

#include "base/import.h"

/** @addtogroup filesGroup */
/// @{

/***** DEFINE **********************************/

/** @brief Macro for file format registration in documenation.
* Call: GROOPS_REGISTER_FILEFORMAT(Title, "type").
* A string named "docstring{Title}" must exist. */
#define GROOPS_REGISTER_FILEFORMAT(_title, _type)\
class _FileFormat##_title : public FileFormat\
{\
public:\
  std::string title() const override {return #_title;}\
  std::string type()  const override {return _type;}\
  std::string documentation() const override {return docstring##_title;}\
};\
static _FileFormat##_title _fileFormat##_title;\
// end macro

/***** CLASS ***********************************/

/** @brief Class registration in documentation.
* Use macro GROOPS_REGISTER_FILEFORMAT. */
class FileFormat
{
public:
  FileFormat() {list(this);}
  virtual ~FileFormat() {}

  virtual std::string title() const = 0;
  virtual std::string type()  const = 0;
  virtual std::string documentation() const = 0;

  static std::vector<FileFormat*> list(FileFormat *fileFormat=nullptr)
  {
    static std::vector<FileFormat*> list_;
    if(fileFormat != nullptr)
      list_.push_back(fileFormat);
    return list_;
  }

  static void sort(std::vector<FileFormat*> &list)
  {
    std::sort(list.begin(), list.end(), [](FileFormat *a, FileFormat *b) {return a->type() < b->type();});
  }
};

/***********************************************/

/// @}

#endif /* __GROOPS__ */

