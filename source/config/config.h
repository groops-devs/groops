/***********************************************/
/**
* @file config.h
*
* @brief Read a configuration file or writes the configuration options in a XSD-Schema.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2008-07-22
*/
/***********************************************/

#ifndef __GROOPS_CONFIG__
#define __GROOPS_CONFIG__

#include "base/import.h"
#include "parser/xml.h"
#include "parser/expressionParser.h"
#include "inputOutput/fileName.h"
#include "parallel/parallel.h"

/** @addtogroup configGroup */
/// @{

/***** TYPES ***********************************/

class Loop;
typedef std::shared_ptr<Loop> LoopPtr;

/***** CLASS ***********************************/

/** @brief Read a configuration file or writes the configuration options in a XSD-Schema. */
class Config
{
public:
  enum Appearance      {MUSTSET, DEFAULT, OPTIONAL};
  enum AppearanceCount {ONCE, UNBOUNDED};
  enum ComplexType     {CHOICE, SEQUENCE, COMPLEXTYPE};

  /** @brief Creates a Configuration from a XML file.
  * Elements from @a Config can be read with @a readConfig. */
  Config(FileName &fileName, const std::map<std::string, std::string> &commandlineGlobals = std::map<std::string, std::string>());

  Config(Config &&)                 = default; //!< Move is allowed
  Config(const Config &)            = delete;  //!< Copy constructor not allowed
  Config &operator=(const Config &) = delete;  //!< assignment not allowed

  /** @brief Gives the parsed content of the node with @a name. */
  Bool getConfigValue(const std::string &name, const std::string &type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation, std::string &text);

  /** @brief Gives the parsed content of the node with @a name. */
  Bool getConfigValue(const std::string &name, const std::string &type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation, Double &v);

  /** @brief Removes the node with @a name from @a config and stores it in @a conf. */
  Bool getConfig(const std::string &name, Config::Appearance mustSet, Config &conf);

  /** @brief Removes all nodes with @a name from @a config and stores it in @a conf.
  * Number of elements is unknown as loops and conditions are not evaluated.  */
  Bool getUnboundedConfig(const std::string &name, Config &conf);

  /** @brief Reads the first variable(s) @p var.
  * The nodes are not removed from the config.
  * Should be used after @a readConfigLater. */
  template<typename T> void read(T &var, VariableList &variableList) const;

  /** @brief Get the current global variable list. */
  VariableList &getVarList() {return varList;}

  // --- Schema mode ----

  /** @brief Default Constructor.
  * Schema mode. Only used by generateDocumentation. */
  Config();

  /** @brief Writes Configuration options to a XSD-Schema file.
  * Every program and class will be called. A call of @a readConfig do
  * not read any variables but identify the configuration options.
  * These options will be saved in a XSD-Schema file. */
  static void writeSchema(const std::string &fileName);

  /** @brief Register an element to a XSD-Schema. */
  void xselement(const std::string &name, const std::string &type, Config::Appearance mustSet, Config::AppearanceCount count, const std::string &defaultValue, const std::string &annotation);

  /** @brief XML node with elements for documentationTable. */
  XmlNodePtr table();

  /** @brief nested name of current xmlNode. */
  std::string currentNodeName() const;

protected:
  // Internal class
  class StackNode
  {
  public:
    std::string  name;
    ComplexType  type;
    XmlNodePtr   xmlNode;        // complex node with children
    XmlNodePtr   xmlLastChild;   // last processed child
    LoopPtr      loopPtr;
    Bool         loopNext;       // move to next iteration before
    VariableList loopVarListOld; // varList without loop variables

    StackNode(XmlNodePtr _xmlNode, ComplexType _type, const std::string &_name) : name(_name), type(_type), xmlNode(_xmlNode) {}
  };

  std::stack<StackNode>  stack;
  XmlNodePtr             global;
  VariableList           varList;
  Bool                   createSchema;

  // stack management
  void         push(XmlNodePtr xmlNode, ComplexType type, const std::string &name=std::string());
  void         pop();

  // normal mode
  XmlNodePtr   resolveLink(XmlNodePtr xmlNode) const;
  XmlNodePtr   getChild(const std::string &name, Bool remove=TRUE);
  XmlNodePtr   getChildWithLoopCheck(const std::string &name, Bool remove=TRUE);
  void         notEmptyWarning();
  std::string  copy(Config &config, const VariableList &variableList) const;

  // schema mode
  XmlNodePtr   createSchemaNode(const std::string &name);
  void         xssimpleType(const std::string &name, const std::string &baseType);
  void         xscomplexElement(const std::string &name, Config::ComplexType type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);
  void         addAppInfo(const std::string &text);
public:
  void         setNodeUnbounded(); // need access from readConfig(std::vector<T> &var, ...)

  friend class ProgramConfig;
  friend Bool isCreateSchema(Config &config);
  friend void renameDeprecatedConfig(Config &config, const std::string &oldName, const std::string &newName, const Time &time);
  friend void renameDeprecatedChoice(Config &config, std::string &type, const std::string &oldName, const std::string &newName, const Time &time);
  friend Bool hasName(Config &config, const std::string &name, Config::Appearance mustSet);
  friend Bool readConfigSequence(Config &config, const std::string &name, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);
  friend void endSequence(Config &config);
  friend Bool readConfigChoice(Config &config, const std::string &name, std::string &choice, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);
  friend Bool readConfigChoiceElement(Config &config, const std::string &name, const std::string &choice, const std::string &annotation);
  friend void endChoice(Config &config);
};

/***** CLASS ***********************************/

/** @brief Config elements of a program.
* It can be read with @a readConfig().
* The config is evaluated and the is program executed with @a run().
* @ingroup config */
class ProgramConfig : public Config
{
public:
  void run(VariableList &variableList, Parallel::CommunicatorPtr comm) const;
};

/***** FUNCTIONS ***********************************/

/** @brief Check if functions are called in schema mode.
* No code is allowed except the readConfig functions in schema mode.
* @ingroup config */
Bool isCreateSchema(Config &config);

/** @brief Has Config an element with @a name?
* If @a mustSet is MUSTSET and name is not found, an exception is thrown.
* @ingroup config */
Bool hasName(Config &config, const std::string &name, Config::Appearance mustSet=Config::OPTIONAL);

/** @brief Indicate a renamed config element.
* Must be called before the @a readConfig functions.
* @ingroup config */
void renameDeprecatedConfig(Config &config, const std::string &oldName, const std::string &newName, const Time &time);

/** @brief Indicate a renamed choice element.
* Must be called after @a readConfigChoice and before the @a readConfigChoiceElement functions.
* @ingroup config */
void renameDeprecatedChoice(Config &config, std::string &type, const std::string &oldName, const std::string &newName, const Time &time);

/** @brief Starts a group of readConfig elements.
* If @a mustSet is MUSTSET and name is not found, an exception is thrown.
* @ingroup config */
Bool readConfigSequence(Config &config, const std::string &name, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/** @brief Ends a group of readConfig elements.
* Must be called only if the @a readConfigSequence return TRUE before.
* @ingroup config */
void endSequence(Config &config);

/** @brief Starts a choice.
* If result is TRUE the selected choice is given in @a choice.
* If @a mustSet is MUSTSET and name is not found, an exception is thrown.
* @ingroup config */
Bool readConfigChoice(Config &config, const std::string &name, std::string &choice, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/** @brief Check the selected choice.
* Returns TRUE if @a name is the selected @a choice.
* In schema mode all possible @a readConfigChoiceElement must be called.
* @ingroup config */
Bool readConfigChoiceElement(Config &config, const std::string &name, const std::string &choice, const std::string &annotation="");

/** @brief Ends a group of choice elements.
* Must be called only if the @a readConfigChoice return TRUE before.
* @ingroup config */
void endChoice(Config &config);

/** @brief Reads a variable from @a config.
* If @a name is in the @a config, the value is set in @a var and TRUE is returned.
* If @a name is not in the @a config or the content is empty, FALSE is returned. In this case the value of @a var depends on @a mustSet:
* @a mustSet = Config::MUSTSET: an exception is thrown.
* @a mustSet = Config::OPTIONAL: @a var is left untouched and @a defaultValue is not used (@a defaultValue might be set as hint for groopsGui).
* @a mustSet = Config::DEFAULT: @a var is set to @a defaultValue (or in case an empty class is created).
* @ingroup config */
template<typename T> Bool readConfig(Config &config, const std::string &name, T &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/** @brief Removes the config of a variable from @a config and stores it in configNew.
* It can be later evluated with configNew.read(). @a var is untouched and only used to determine the type.
* @ingroup config */
template<typename T> Bool readConfigLater(Config &config, const std::string &name, T &var, Config &configNew, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/** @brief Removes the config of a variable from @a config and stores it in configNew.
* It can be later evluated with configNew.read(). @a var is untouched and only used to determine the type.
* @ingroup config */
template<typename T> Bool readConfigLater(Config &config, const std::string &name, std::vector<T> &var, std::vector<Config> &configNew, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/
/***** INLINES   *******************************/
/***********************************************/

template<typename T> inline void Config::read(T &var, VariableList &variableList) const
{
  try
  {
    Config config;
    const std::string name = copy(config, variableList);
    readConfig(config, name, var, MUSTSET, "", "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<typename T> inline Bool readConfig(Config &config, const std::string &name, std::vector<T> &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(isCreateSchema(config))
  {
    if(mustSet == Config::DEFAULT)
      throw(Exception("In readConfig("+config.currentNodeName()+"."+name+"): Config::DEFAULT is not allowed for std::vector. Please change to Config::OPTIONAL"));
    T tmp;
    readConfig(config, name, tmp, mustSet, defaultValue, annotation);
    config.setNodeUnbounded();
    return FALSE;
  }

  if(hasName(config, name, mustSet))
    do
    {
      var.resize(var.size()+1);
      readConfig(config, name, var.back(), mustSet, defaultValue, annotation);
    }
    while(hasName(config, name, Config::OPTIONAL));

  return (var.size()!=0);
}

/***********************************************/

template<typename T> inline Bool readConfigLater(Config &config, const std::string &name, T &var, Config &configNew, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(isCreateSchema(config))
  {
    readConfig(config, name, var, mustSet, defaultValue, annotation);
    return FALSE;
  }

  // get type and unbounded
  Config configXsd;
  readConfig(configXsd, name, var, mustSet, defaultValue, annotation);
  XmlNodePtr xmlNode = configXsd.table()->getNextChild();

  if(xmlNode->getAttribute("maxOccurs"))
    return config.getUnboundedConfig(name, configNew);
  return config.getConfig(name, mustSet, configNew);
}

/***********************************************/

template<typename T> inline Bool readConfigLater(Config &config, const std::string &name, std::vector<T> &var, std::vector<Config> &configNew, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(isCreateSchema(config))
  {
    readConfig(config, name, var, mustSet, defaultValue, annotation);
    return FALSE;
  }

  if(hasName(config, name, mustSet))
    do
    {
      configNew.resize(configNew.size()+1);
      config.getConfig(name, mustSet, configNew.back());
    }
    while(hasName(config, name, Config::OPTIONAL));

  return (configNew.size()!=0);
}

/***********************************************/

#endif
