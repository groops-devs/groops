/***********************************************/
/**
* @file program.cpp
*
* @brief Interface for applications in groops.
*
* @author Torsten Mayer-Guerr
* @date 2008-07-27
*
*/
/***********************************************/

// define GROOPS_TAGS to autogenerate tagStrings
#define GROOPS_STRINGIFY(x) ,#x
#define GROOPS_TAGS(tag_, ...)\
enum Tags {tag_, __VA_ARGS__};\
const char *tagStrings[] = {#tag_ _GROOPS_FOR_EACH(GROOPS_STRINGIFY, __VA_ARGS__)};

#include "base/import.h"
#include "config/configRegister.h"
#include "program.h"

/***********************************************/

std::vector<Program::Program*> Program::Program::programList(Program *program)
{
  static std::vector<Program*> list;
  if(program != nullptr)
    list.push_back(program);
  return list;
}

/***********************************************/

// sort programlist according to first tag, followed by alphabetical order of names
void Program::Program::sortList(std::vector<Program*> &list)
{
  std::sort(list.begin(), list.end(), [](Program *a, Program *b) {return (a->tags().at(0) != b->tags().at(0)) ? (a->tags().at(0) < b->tags().at(0)) : (a->name() < b->name());});
}

/***********************************************/

std::vector<Program::RenamedProgram::Renamed>
  Program::RenamedProgram::renamedList(const Renamed &renamed)
{
  static std::vector<Renamed> list;
  if(!renamed.oldName.empty())
    list.push_back(renamed);
  return list;
}

/***********************************************/
