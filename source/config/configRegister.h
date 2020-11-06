/***********************************************/
/**
* @file configRegister.h
*
* @brief Interface for classes in groops.
*
* @author Torsten Mayer-Guerr
* @date 2018-05-20
*
*/
/***********************************************/

#ifndef __GROOPS_CONFIGREGISTER__
#define __GROOPS_CONFIGREGISTER__

#include "base/import.h"
#include "config/config.h"
#include "config/generateDocumentation.h"

/** @addtogroup configGroup */
/// @{

/***** DEFINE **********************************/

// workaround for Microsoft compiler
// https://renenyffenegger.ch/notes/development/languages/C-C-plus-plus/preprocessor/macros/__VA_ARGS__/index
#define _GROOPS_PASS_ON(...) __VA_ARGS__

/// calls func for every following argument, e.g. _GROOPS_FOR_EACH(test, a, b, c) -> test(a) test(b) test(c)
#define _GROOPS_FOR_EACH(_func, ...) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GRROPS_FUNC)(_GROOPS_func, _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_NARG)(__VA_ARGS__)))(_func, __VA_ARGS__))

// Macro magic: Overloading Macro on Number of Arguments
// https://stackoverflow.com/questions/11761703/overloading-macro-on-number-of-arguments/11763277#11763277
#define _GROOPS_ARG_N(_01, _02, _03, _04, _05, _06, _07, _08, _09, _10,\
                      _11, _12, _13, _14, _15, _16, _17, _18, _19, _20,\
                      _21, _22, _23, _24, _25, _26, _27, _28, _29, _30,\
                      _31, _32, _33, _34, _35, _36, _37, _38, _39, _40,\
                      _41, _42, _43, _44, _45, _46, _47, _48, _49, N, ...) N

#define _GROOPS_SEQUENCE 49, 48, 47, 46, 45, 44, 43, 42, 41, 40,\
                         39, 38, 37, 36, 35, 34, 33, 32, 31, 30,\
                         29, 28, 27, 26, 25, 24, 23, 22, 21, 20,\
                         19, 18, 17, 16, 15, 14, 13, 12, 11, 10,\
                         09, 08, 07, 06, 05, 04, 03, 02, 01, 00

#define _GROOPS_func00(_func)
#define _GROOPS_func01(_func, x)      _func(x)
#define _GROOPS_func02(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func01)(_func, __VA_ARGS__))
#define _GROOPS_func03(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func02)(_func, __VA_ARGS__))
#define _GROOPS_func04(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func03)(_func, __VA_ARGS__))
#define _GROOPS_func05(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func04)(_func, __VA_ARGS__))
#define _GROOPS_func06(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func05)(_func, __VA_ARGS__))
#define _GROOPS_func07(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func06)(_func, __VA_ARGS__))
#define _GROOPS_func08(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func07)(_func, __VA_ARGS__))
#define _GROOPS_func09(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func08)(_func, __VA_ARGS__))
#define _GROOPS_func10(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func09)(_func, __VA_ARGS__))
#define _GROOPS_func11(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func10)(_func, __VA_ARGS__))
#define _GROOPS_func12(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func11)(_func, __VA_ARGS__))
#define _GROOPS_func13(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func12)(_func, __VA_ARGS__))
#define _GROOPS_func14(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func13)(_func, __VA_ARGS__))
#define _GROOPS_func15(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func14)(_func, __VA_ARGS__))
#define _GROOPS_func16(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func15)(_func, __VA_ARGS__))
#define _GROOPS_func17(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func16)(_func, __VA_ARGS__))
#define _GROOPS_func18(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func17)(_func, __VA_ARGS__))
#define _GROOPS_func19(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func18)(_func, __VA_ARGS__))
#define _GROOPS_func20(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func19)(_func, __VA_ARGS__))
#define _GROOPS_func21(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func20)(_func, __VA_ARGS__))
#define _GROOPS_func22(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func21)(_func, __VA_ARGS__))
#define _GROOPS_func23(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func22)(_func, __VA_ARGS__))
#define _GROOPS_func24(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func23)(_func, __VA_ARGS__))
#define _GROOPS_func25(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func24)(_func, __VA_ARGS__))
#define _GROOPS_func26(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func25)(_func, __VA_ARGS__))
#define _GROOPS_func27(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func26)(_func, __VA_ARGS__))
#define _GROOPS_func28(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func27)(_func, __VA_ARGS__))
#define _GROOPS_func29(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func28)(_func, __VA_ARGS__))
#define _GROOPS_func30(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func29)(_func, __VA_ARGS__))
#define _GROOPS_func31(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func30)(_func, __VA_ARGS__))
#define _GROOPS_func32(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func31)(_func, __VA_ARGS__))
#define _GROOPS_func33(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func32)(_func, __VA_ARGS__))
#define _GROOPS_func34(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func33)(_func, __VA_ARGS__))
#define _GROOPS_func35(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func34)(_func, __VA_ARGS__))
#define _GROOPS_func36(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func35)(_func, __VA_ARGS__))
#define _GROOPS_func37(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func36)(_func, __VA_ARGS__))
#define _GROOPS_func38(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func37)(_func, __VA_ARGS__))
#define _GROOPS_func39(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func38)(_func, __VA_ARGS__))
#define _GROOPS_func40(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func39)(_func, __VA_ARGS__))
#define _GROOPS_func41(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func40)(_func, __VA_ARGS__))
#define _GROOPS_func42(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func41)(_func, __VA_ARGS__))
#define _GROOPS_func43(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func42)(_func, __VA_ARGS__))
#define _GROOPS_func44(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func43)(_func, __VA_ARGS__))
#define _GROOPS_func45(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func44)(_func, __VA_ARGS__))
#define _GROOPS_func46(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func45)(_func, __VA_ARGS__))
#define _GROOPS_func47(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func46)(_func, __VA_ARGS__))
#define _GROOPS_func48(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func47)(_func, __VA_ARGS__))
#define _GROOPS_func49(_func, x, ...) _func(x) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_func48)(_func, __VA_ARGS__))

#define _GRROPS_FUNC_(name, n) name##n
#define _GRROPS_FUNC(name, n) _GRROPS_FUNC_(name, n)
#define _GROOPS_NARG(...)  _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_NARG_)(__VA_ARGS__,_GROOPS_SEQUENCE))
#define _GROOPS_NARG_(...) _GROOPS_PASS_ON(_GROOPS_PASS_ON(_GROOPS_ARG_N)(__VA_ARGS__))

/***********************************************/

/** @brief Macro for class registration in schema and documenation.
* Call: GROOPS_REGISTER_CLASS_WITHOUT_SUBS(ClassName, "schemaType").
* A string named "docstring{ClassName}" must exist.
* The class must provide a static function "create(Config &config, const std::string &name)".
* @ingroup config */
#define GROOPS_REGISTER_CLASS_WITHOUT_SUBS(_className, _typeName)\
_GROOPS_REGISTER_CLASS(_className, _typeName)\
inline void _Class##_className::generateDocumentation(Documentation &documentation) const\
{\
  Config config;\
  _className::create(config, typeName());\
  documentation.writeText(docstring##_className);\
  documentation.writeConfigTable(config);\
}\
// end macro

/***********************************************/

/** @brief Macro for class registration in schema and documenation.
* Call: GROOPS_REGISTER_CLASS(ClassName, "schemaType", subClasses...)
* A string named "docstring{ClassName}" must exist for the class and subclasses.
* The class must provide a static function "create(Config &config, const std::string &name)".
* The sub classes must have a constructor of type {ClassName}(Config &config).
* @ingroup config */
#define GROOPS_REGISTER_CLASS(_className, _typeName, ...)\
_GROOPS_REGISTER_CLASS(_className, _typeName)\
inline void _Class##_className::generateDocumentation(Documentation &documentation) const\
{\
  documentation.writeText(docstring##_className);\
  Config config;\
  _GROOPS_FOR_EACH(_GROOPS_CLASS_DOCUMENTATION, __VA_ARGS__)\
}\
// end macro

/***** DEFINE **********************************/

/** @brief Register an old name of a renamed class.
* This macro can be used to keep old config files working.
* @ingroup config */
#define GROOPS_RENAMED_CLASS(_oldname, _newname, _time) \
class _RenamedSchemaClass##_oldname : public RenamedSchemaClass\
{\
public:\
  _RenamedSchemaClass##_oldname(const Renamed &renamed) : RenamedSchemaClass(renamed) {}\
};\
static _RenamedSchemaClass##_oldname renamedSchemaClass##_oldname(RenamedSchemaClass::Renamed(#_oldname, #_newname, _time));\
// end macro

/***********************************************/

/** @brief Generate the readConfig function for a class.
* Call: GROOPS_READCONFIG_CLASS(ClassName, "schemaType").
* @ingroup config */
#define GROOPS_READCONFIG_CLASS(_className, _typeName)\
template<> Bool readConfig(Config &config, const std::string &name, _className##Ptr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)\
{\
  try\
  {\
    if(isCreateSchema(config))\
    {\
      if(mustSet == Config::DEFAULT)\
        throw(Exception("In readConfig("+config.currentNodeName()+"."+name+", type="+#_typeName+"): Config::DEFAULT is not allowed. Please change to Config::OPTIONAL"));\
      config.xselement(name, _typeName, mustSet, Config::ONCE, defaultValue, annotation); \
      return FALSE;\
    }\
\
    if(!hasName(config, name, mustSet))\
      return FALSE;\
    var = _className::create(config, name);\
    return TRUE;\
  }\
  catch(std::exception &e)\
  {\
    GROOPS_RETHROW(e)\
  }\
}\
// end macro

/***********************************************/

/** @brief Generate the readConfig function for an unbounded class.
* Call: GROOPS_READCONFIG_UNBOUNDED_CLASS(ClassName, "schemaType").
* @ingroup config */
#define GROOPS_READCONFIG_UNBOUNDED_CLASS(_className, _typeName)\
template<> Bool readConfig(Config &config, const std::string &name, _className##Ptr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)\
{\
  try\
  {\
    if(isCreateSchema(config))\
    {\
      config.xselement(name, _typeName, mustSet, Config::UNBOUNDED, defaultValue, annotation); \
      return FALSE;\
    }\
\
    Bool found = hasName(config, name, mustSet);\
    if(found || (mustSet == Config::DEFAULT))\
      var = _className##Ptr(new _className(config, name));\
    return found;\
  }\
  catch(std::exception &e)\
  {\
    GROOPS_RETHROW(e)\
  }\
}\
// end macro

/***********************************************/

// Internal macro
#define _GROOPS_REGISTER_CLASS(_className, _typeName)\
class _Class##_className : public SchemaClass\
{\
public:\
  std::string typeName() const {return _typeName;}\
  void registerConfigSchema(Config &config) const {_className::create(config, typeName());}\
  void generateDocumentation(Documentation &documentation) const;\
};\
static _Class##_className _class##_className;\
// end macro


/***********************************************/

// Internal macro
#define _GROOPS_CLASS_DOCUMENTATION(_className)\
{\
  _className tmp(config);\
  documentation.writeText(docstring##_className);\
  documentation.writeConfigTable(config);\
}\
// end macro

/***** CLASS ***********************************/

/** @brief Class registration in schema and documentation.
* Use macros GROOPS_REGISTER_CLASS and GROOPS_REGISTER_CLASS_WITHOUT_SUBS.
* @ingroup config */
class SchemaClass
{
public:
  SchemaClass() {classList(this);}
  virtual ~SchemaClass() {}

  virtual std::string typeName() const = 0;
  virtual void registerConfigSchema(Config &config) const = 0;
  virtual void generateDocumentation(Documentation &documentation) const = 0;

  static std::vector<SchemaClass*> classList(SchemaClass *schemaClass=nullptr)
  {
    static std::vector<SchemaClass*> list;
    if(schemaClass != nullptr)
      list.push_back(schemaClass);
    return list;
  }

  static void sort(std::vector<SchemaClass*> &list)
  {
    std::sort(list.begin(), list.end(), [](SchemaClass *a, SchemaClass *b) {return a->typeName() < b->typeName();});
  }
};

/***** CLASS ***********************************/

/** @brief Interface for renamed classes.
* @ingroup config */
class RenamedSchemaClass
{
public:
  class Renamed
  {
  public:
    std::string oldName;
    std::string newName;
    Time        time;
    Renamed(const std::string &oldName_, const std::string &newName_, const Time &time_) : oldName(oldName_), newName(newName_), time(time_) {}
  };

  RenamedSchemaClass(const Renamed &renamed) {renamedList(renamed);}
  virtual ~RenamedSchemaClass() {}

  static std::vector<Renamed> renamedList(const Renamed &renamed = Renamed("", "", Time()))
  {
    static std::vector<Renamed> list;
    if(!renamed.oldName.empty())
      list.push_back(renamed);
    return list;
  }
};

/***********************************************/

/// @}

#endif /* __GROOPS__ */

