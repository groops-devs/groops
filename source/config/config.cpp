/***********************************************/
/**
* @file config.cpp
*
* @brief Read a configuration file or writes the configuration options in a XSD-Schema.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2008-07-22
*/
/***********************************************/

#include "base/import.h"
#include "config/configRegister.h"
#include "base/doodson.h"
#include "base/gnssType.h"
#include "inputOutput/file.h"
#include "parser/xml.h"
#include "parser/stringParser.h"
#include "parser/expressionParser.h"
#include "parallel/parallel.h"
#include "classes/condition/condition.h"
#include "classes/loop/loop.h"
#include "programs/program.h"
#include "config.h"

/***********************************************/
/*** Stack management **************************/
/***********************************************/

void Config::push(XmlNodePtr xmlNode, ComplexType type, const std::string &name)
{
  std::string name_ = (name.empty()) ? Config::currentNodeName()+"."+xmlNode->getName() : name;
  stack.emplace(xmlNode, type, name_);
}

/***********************************************/

void Config::pop()
{
  if(stack.top().loopPtr) // finish old loop
    varList = stack.top().loopVarListOld; // restore old varList
  stack.pop();
}

/***********************************************/

std::string Config::currentNodeName() const
{
  return stack.top().name;
}

/***********************************************/
/*** Normal mode (readConfig) ******************/
/***********************************************/

Config::Config(FileName &fileName, const std::map<std::string, std::string> &commandlineGlobals)
{
  try
  {
    createSchema = FALSE;

    InFile stream(fileName);
    XmlNodePtr root = XmlNode::read(stream);

    // global: replace typename with label
    // -----------------------------------
    global = root->getChild("global");
    if(!global)
      global = XmlNode::create("global");
    XmlNodePtr globalFile = XmlNode::create("global");
    while(global->hasChildren())
    {
      XmlNodePtr xmlNode = global->getNextChild();
      XmlAttrPtr label = xmlNode->getAttribute("label");
      if(label)
        xmlNode->setName(label->getText());
      globalFile->addChild(xmlNode);
    }

    // add global elements from command line
    // -------------------------------------
    for(auto it=commandlineGlobals.begin(); it!=commandlineGlobals.end(); it++)
    {
      // replace or add new global variable
      XmlNodePtr xmlNode = globalFile->findChild(it->first);
      if(xmlNode)
      {
        xmlNode->setText(it->second);
        xmlNode->getAttribute("link"); // remove link
      }
      else
        writeXml(globalFile, it->first, it->second);
    }

    // Additional global input files
    // -----------------------------
    std::vector<FileName> fileNameGlobal(root->getChildCount("inputfileGlobal"));
    for(UInt i=0; i<fileNameGlobal.size(); i++)
      fileNameGlobal.at(i) = root->getChild("inputfileGlobal")->getText();
    for(UInt i=fileNameGlobal.size(); i-->0;)
    {
      Config config(fileNameGlobal.at(i));
      while(config.global->hasChildren())
      {
        XmlNodePtr xmlNode = config.global->getNextChild();
        if(!globalFile->findChild(xmlNode->getName()))
          globalFile->prependChild(xmlNode);
      }
    }

    // search global variables
    // -----------------------
    while(globalFile->hasChildren())
    {
      XmlNodePtr xmlNode = globalFile->getNextChild();
      global->addChild(xmlNode);
      if(!xmlNode->hasChildren()) // not complex type?
      {
        XmlAttrPtr link = xmlNode->findAttribute("link");
        addVariable(xmlNode->getName(), ((link) ? "{"+link->getText()+"}" : xmlNode->getText()), varList);
      }
    }

    global = resolveLink(global);
    root = resolveLink(root);
    push(root, SEQUENCE, "groops");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

XmlNodePtr Config::resolveLink(XmlNodePtr xmlNode) const
{
  try
  {
    if(!xmlNode)
      return xmlNode;

    XmlAttrPtr link = xmlNode->findAttribute("link");
    if(link)
    {
      XmlNodePtr xmlNodeGlobal = global->findChild(link->getText());
      if(!xmlNodeGlobal || xmlNodeGlobal->getAttribute("resolving"))
        throw(Exception(std::string("cannot resolve link: ")+xmlNode->getName()+" -> "+link->getText()+"'"));

      XmlNodePtr xmlNodeNew = xmlNodeGlobal->clone();

      // add temporary resolving attribute to prevent link loops
      writeAttribute(xmlNodeGlobal, "resolving", 1);

      xmlNodeNew = resolveLink(xmlNodeNew);
      xmlNodeNew->setName(xmlNode->getName());
      while(xmlNode->hasAttribute())
        xmlNodeNew->addAttribute(xmlNode->getNextAttribute());

      // remove attributes as link was resolved successfully
      xmlNodeGlobal->getAttribute("resolving");
      xmlNodeNew->getAttribute("link");

      xmlNode = xmlNodeNew;
    }

    // test children
    for(XmlNodePtr &child : xmlNode->getChildren())
      child = resolveLink(child);

    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

XmlNodePtr Config::getChild(const std::string &name, Bool remove)
{
  try
  {
    XmlNodePtr xmlNode = getChildWithLoopCheck(name, remove);
    if(!xmlNode)
      return xmlNode;

    XmlAttrPtr attr = xmlNode->getAttribute("condition");
    if(!attr)
      return xmlNode;

    // expand condition
    XmlNodePtr xmlNodeCondition = global->findChild(attr->getText());
    if(!xmlNodeCondition)
      throw(Exception(std::string("cannot resolve condition link: '")+xmlNode->getName()+"' -> '"+attr->getText()+"'"));
    XmlNodePtr tmp = XmlNode::create("tmp");
    tmp->addChild(xmlNodeCondition->clone()); // make copy
    push(tmp, SEQUENCE, currentNodeName());
    ConditionPtr condition;
    readConfig(*this, xmlNodeCondition->getName(), condition, Config::MUSTSET, "", "");
    pop();

    // check condition
    if(!condition->condition(varList))
    {
      if(!remove)
        getChildWithLoopCheck(name, TRUE/*remove*/); // node disabled -> can always be removed
      return getChild(name, remove);    // node disabled -> try next element
    }

    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

XmlNodePtr Config::getChildWithLoopCheck(const std::string &name, Bool remove)
{
  try
  {
    // finish old iteration
    // --------------------
    if(stack.top().loopPtr && stack.top().loopNext && !stack.top().loopPtr->iteration(varList))
    {
      varList     = stack.top().loopVarListOld; // restore old varList
      stack.top().xmlNode->getChild(stack.top().xmlLastChild->getName());
      stack.top().loopPtr = LoopPtr(nullptr);
    }

    // if not loop: get new child
    // --------------------------
    if(!stack.top().loopPtr)
    {
      // get new child
      XmlNodePtr xmlNode = stack.top().xmlNode->findChild(name);
      stack.top().xmlLastChild = xmlNode;
      if(xmlNode == nullptr)
        return xmlNode;

      XmlAttrPtr attr = xmlNode->getAttribute("loop");
      if(attr == nullptr)
      {
        if(remove)
          stack.top().xmlNode->getChild(name);
        return xmlNode;
      }

      // expand loop
      XmlNodePtr xmlNodeLoop = global->findChild(attr->getText());
      if(!xmlNodeLoop)
        throw(Exception(std::string("cannot resolve loop link: ")+xmlNode->getName()+" -> "+attr->getText()+"'"));
      XmlNodePtr tmp = XmlNode::create("tmp");
      tmp->addChild(xmlNodeLoop->clone()); // make copy
      push(tmp, SEQUENCE, currentNodeName());
      LoopPtr loopPtr;
      readConfig(*this, xmlNodeLoop->getName(), loopPtr, Config::MUSTSET, "", "");
      pop();

      // init loop
      stack.top().loopVarListOld = varList;
      stack.top().loopPtr        = loopPtr;
      if(!stack.top().loopPtr->iteration(varList)) // empty loop?
      {
        stack.top().loopNext = TRUE;
        stack.top().xmlNode->getChild(name); // remove child
        return getChild(name, remove);       // node disabled -> try next element
      }
    }

    // Now we are in a loop
    // --------------------
    if(stack.top().xmlLastChild->getName() != name)
      throw(Exception("loop error"));
    stack.top().loopNext = remove;
    return stack.top().xmlLastChild->clone();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Config::getConfigValue(const std::string &name, const std::string &type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation, std::string &text)
{
  if(createSchema)
  {
    if((mustSet == Config::DEFAULT) && defaultValue.empty())
      throw(Exception("In readConfig("+currentNodeName()+"."+name+", type="+type+"): Config::DEFAULT with empty defaultValue"));
    xselement(name, type, mustSet, ONCE, defaultValue, annotation);
    return FALSE;
  }

  XmlNodePtr child = getChild(name);
  if(child != nullptr)
  {
    if(child->hasChildren())
      throw(Exception(currentNodeName()+"."+child->getName()+" with unexpected children"));
    text = child->getText();
  }

  if(text.empty() && (mustSet == MUSTSET))
    throw(Exception("config element '"+currentNodeName()+"' must contain '"+name+"'"));
  if(text.empty() && (mustSet == DEFAULT))
    text = defaultValue;

  Bool resolved;
  text = StringParser::parse(name, text, varList, resolved);
  return !text.empty();
}

/***********************************************/

Bool Config::getConfigValue(const std::string &name, const std::string &type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation, Double &v)
{
  std::string text;
  try
  {
    if(!getConfigValue(name, type, mustSet, defaultValue, annotation, text))
      return FALSE;
    v = Expression::parse(text)->evaluate(varList);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("'"+name+"' = ' "+text+"'", e)
  }
}

/***********************************************/

Bool Config::getConfigValue(const std::string &name, const std::string &type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation, Config &conf)
{
  try
  {
    if(createSchema)
    {
      xselement(name, type, mustSet, ONCE, defaultValue, annotation);
      return FALSE;
    }

    XmlNodePtr xmlNode = XmlNode::create(name);
    XmlNodePtr child   = getChild(name);
    if(child)
      xmlNode->addChild(child);
    else if(mustSet == MUSTSET)
      throw(Exception("config element '"+currentNodeName()+"' must contain '"+name+"'"));

    // make copy
    conf.push(xmlNode, Config::SEQUENCE, currentNodeName());
    conf.createSchema = FALSE;
    conf.global       = global;
    conf.varList      = varList;

    return child != nullptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("'"+name+"'", e)
  }
}

/***********************************************/

Bool Config::getUnboundedConfigValues(const std::string &name, const std::string &type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation, Config &conf)
{
  try
  {
    if(createSchema)
    {
      xselement(name, type, mustSet, UNBOUNDED, defaultValue, annotation);
      return FALSE;
    }

    // finish old iteration
    // --------------------
    if(stack.top().loopPtr && stack.top().loopNext && !stack.top().loopPtr->iteration(varList))
    {
      varList     = stack.top().loopVarListOld; // restore old varList
      stack.top().xmlNode->getChild(stack.top().xmlLastChild->getName());
      stack.top().loopPtr = LoopPtr(nullptr);
    }

    if(stack.top().loopPtr)
      throw(Exception("Unexpected loop attribute"));

    XmlNodePtr xmlNode = XmlNode::create(name);
    for(;;)
    {
      XmlNodePtr child = stack.top().xmlNode->getChild(name);
      if(!child)
        break;
      xmlNode->addChild(child);
    }

    // make copy
    conf.push(xmlNode, Config::SEQUENCE, currentNodeName());
    conf.createSchema = FALSE;
    conf.global       = global;
    conf.varList      = varList;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("'"+name+"'", e)
  }
}

/***********************************************/

void Config::notEmptyWarning()
{
  try
  {
    if(!Parallel::isMaster())
      return;
    if(stack.top().xmlNode->hasChildren())
    {
      logWarning<<"*** Warning: unknown variables in '"<<currentNodeName()<<"' (ignored):"<<Log::endl;
      while(stack.top().xmlNode->hasChildren())
        logWarning<<"  name = "<<stack.top().xmlNode->getNextChild()->getName()<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Config::copy(Config &config, const VariableList &variableList)
{
  try
  {
    config.push(stack.top().xmlNode->clone(), Config::SEQUENCE, currentNodeName());
    config.createSchema = FALSE;
    config.global       = global;
    config.varList      = varList;
    config.varList     += variableList;
    return stack.top().xmlNode->getName();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/*** Schema mode *******************************/
/***********************************************/

// Only used by generateDocumentation
Config::Config()
{
  try
  {
    createSchema = TRUE;
    XmlNodePtr rootNode = XmlNode::create("xs:schema");
    writeAttribute(rootNode, "xmlns:xs",             "http://www.w3.org/2001/XMLSchema");
    writeAttribute(rootNode, "elementFormDefault",   "qualified");
    writeAttribute(rootNode, "attributeFormDefault", "unqualified");
    push(rootNode, SEQUENCE, "xs:schema");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Config::writeSchema(const std::string &fileName)
{
  try
  {
    Config config;
    XmlNodePtr rootNode = config.stack.top().xmlNode;
    config.stack.top().type = Config::COMPLEXTYPE;

    // all types are string types as they can contain expressions
    config.xssimpleType("int",        "xs:string"); //"xs:int");
    config.xssimpleType("uint",       "xs:string"); //"xs:nonNegativeInteger");
    config.xssimpleType("double",     "xs:string"); //"xs:double");
    config.xssimpleType("angle",      "xs:string"); //"xs:double");
    config.xssimpleType("boolean",    "xs:string"); //"xs:boolean");
    config.xssimpleType("time",       "xs:string"); //"xs:double");
    config.xssimpleType("doodson",    "xs:string");
    config.xssimpleType("string",     "xs:string");
    config.xssimpleType("filename",   "xs:string");
    config.xssimpleType("expression", "xs:string");
    config.xssimpleType("gnssType",   "xs:string");

    std::vector<SchemaClass*> classList = SchemaClass::classList();
    SchemaClass::sort(classList);
    for(UInt i=0; i<classList.size(); i++)
    {
      const std::string name = config.stack.top().name;
      classList.at(i)->registerConfigSchema(config);
      if(name != config.stack.top().name)
        throw(Exception("In class "+classList.at(i)->typeName()+": Missing endSequence() or endChoice()?"));
    }

    XmlNodePtr xmlNodeTypes = rootNode->clone(); // type for the global section

    // list of programs
    std::string choice;
    readConfigChoice(config, "programType", choice, Config::MUSTSET, "", "");
    const auto renamedList = Program::RenamedProgram::renamedList();
    for(auto &renamed : renamedList)
      renameDeprecatedChoice(config, choice, renamed.oldName, renamed.newName, renamed.time);
    std::vector<Program::Program*> programList = Program::Program::programList();
    Program::Program::sortList(programList);
    for(auto program : programList)
    {
      readConfigChoiceElement(config, program->name(), choice, program->description());
      const std::string name = config.stack.top().name;
      program->run(config);
      if(name != config.stack.top().name)
        throw(Exception("In program "+program->name()+": Missing endSequence() or endChoice()?"));
      for(UInt idx : program->tags())
        config.addAppInfo("tag: "s+Program::tagStrings[idx]);
    }
    endChoice(config); // end program

    // type definition finished, now the elements
    config.stack.top().type = SEQUENCE;
    config.xscomplexElement("groops", Config::SEQUENCE, Config::MUSTSET, "", "GROOPS (Gravity Recovery Object Oriented Programming System)");
    renameDeprecatedConfig(config, "programme", "program", date2time(2020, 6, 3));

    config.xselement("inputfileGlobal", "filename", Config::OPTIONAL, Config::UNBOUNDED, "", "global settings");

    // global section
    readConfigSequence(config, "global", Config::MUSTSET, "", "global settings");
    const auto renamedClassList = RenamedSchemaClass::renamedList();
    for(auto &renamed : renamedClassList)
      renameDeprecatedConfig(config, renamed.oldName, renamed.newName, renamed.time);
    while(xmlNodeTypes->hasChildren())
    {
      std::string name;
      readAttribute(xmlNodeTypes->getNextChild(), "name", name, TRUE/*mustSet*/);
      config.xselement(name, name, Config::OPTIONAL, Config::UNBOUNDED, "", "");
    }
    endSequence(config); // end global

    // program
    config.xselement("program", "programType", Config::OPTIONAL, Config::UNBOUNDED, "", "");
    endSequence(config); // end <groops>
    OutFile file(fileName);
    XmlNode::write(file, rootNode);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

XmlNodePtr Config::createSchemaNode(const std::string &name)
{
  try
  {
    XmlNodePtr xmlNode = createXmlNode(stack.top().xmlNode, name);
    stack.top().xmlLastChild = xmlNode;
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Config::xssimpleType(const std::string &name, const std::string &base)
{
  try
  {
    XmlNodePtr xmlNode = createSchemaNode("xs:simpleType");
    writeAttribute(xmlNode, "name", name);
    writeAttribute(createXmlNode(xmlNode, "xs:restriction"), "base", base);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Config::xselement(const std::string &name, const std::string &type, Config::Appearance mustSet, Config::AppearanceCount count, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!createSchema)
      throw(Exception("modus of config is not createSchema"));

    XmlNodePtr xmlNode = createSchemaNode("xs:element");
    writeAttribute(xmlNode, "name", name);
    writeAttribute(xmlNode, "type", type);
    if(mustSet != MUSTSET)
      writeAttribute(xmlNode, "minOccurs", "0");
    if(count == UNBOUNDED)
      writeAttribute(xmlNode, "maxOccurs", "unbounded");
    if(!defaultValue.empty())
      writeAttribute(xmlNode, "default", defaultValue);
    if(!annotation.empty())
      writeXml(createXmlNode(xmlNode, "xs:annotation"), "xs:documentation", annotation);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Config::xscomplexElement(const std::string &name, Config::ComplexType type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    XmlNodePtr xmlNode = createSchemaNode((stack.top().type == COMPLEXTYPE) ? "xs:complexType" : "xs:element");
    writeAttribute(xmlNode, "name", name);
    if(!annotation.empty())
      writeXml(createXmlNode(xmlNode, "xs:annotation"), "xs:documentation", annotation);

    if(stack.top().type != Config::COMPLEXTYPE) // is not root level?
    {
      if(mustSet == OPTIONAL)
        writeAttribute(xmlNode, "minOccurs", "0");
      if(!defaultValue.empty())
        writeAttribute(xmlNode, "default", defaultValue);
      xmlNode = createXmlNode(xmlNode, "xs:complexType"); // <xs:element> -> <xs:complexType>
    }

    if(type==SEQUENCE)
      push(createXmlNode(xmlNode, "xs:sequence"), SEQUENCE, currentNodeName()+"."+name);
    else
      push(createXmlNode(xmlNode, "xs:choice"), CHOICE, currentNodeName()+"."+name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Config::addAppInfo(const std::string &text)
{
  try
  {
    StackNode top = stack.top();
    stack.pop();
    if(!stack.empty())
    {
      XmlNodePtr xmlNode = stack.top().xmlLastChild;
      if(!xmlNode)
        throw(Exception(top.name+ " " +text));
      XmlNodePtr xmlAnnotation = xmlNode->findChild("xs:annotation");
      if(!xmlAnnotation)
      {
        xmlAnnotation = XmlNode::create("xs:annotation");
        xmlNode->prependChild(xmlAnnotation);
      }
      writeXml(xmlAnnotation, "xs:appinfo", text);
    }
    stack.push(top);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Config::setNodeUnbounded()
{
  try
  {
    if(!stack.top().xmlLastChild)
      throw(Exception("schemaNode = nullptr"));
    writeAttribute(stack.top().xmlLastChild, "maxOccurs", "unbounded");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

XmlNodePtr Config::table()
{
  try
  {
    return stack.top().xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void ProgramConfig::run(VariableList &variableList)
{
  try
  {
    Config config;
    const std::string name = copy(config, variableList);

    std::string type;
    while(readConfigChoice(config, name, type, OPTIONAL, "", ""))
    {
      for(auto &renamed : Program::RenamedProgram::renamedList())
        renameDeprecatedChoice(config, type, renamed.oldName, renamed.newName, renamed.time);

      for(auto &program : Program::Program::programList())
        if(readConfigChoiceElement(config, program->name(), type, ""))
        {
          std::string comment;
          StackNode top = stack.top();
          stack.pop(); // coment is given in <program> not in <choiceElement>
          XmlAttrPtr attr = stack.top().xmlNode->getAttribute("comment");
          if(attr)
            comment = attr->getText();
          stack.push(top);

          if(comment.empty())
            logStatus<<"--- "<<program->name()<<" ---"<<Log::endl;
          else
          {
            Bool resolved;
            comment = StringParser::parse("comment", comment, getVarList(), resolved);
            logStatus<<"--- "<<program->name()<<" ("<<comment<<") ---"<<Log::endl;
          }
          Parallel::barrier();
          program->run(config);
          break;
        }

      endChoice(config);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

LoopPtr LoopConfig::read(VariableList &variableList)
{
  try
  {
    Config config;
    const std::string name = copy(config, variableList);
    LoopPtr loop;
    readConfig(config, name, loop, MUSTSET, "", "");
    return loop;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/*** Functions *********************************/
/***********************************************/

Bool isCreateSchema(Config &config)
{
  return config.createSchema;
}

/***********************************************/

void renameDeprecatedConfig(Config &config, const std::string &oldName, const std::string &newName, const Time &time)
{
  try
  {
    if(isCreateSchema(config))
    {
      config.addAppInfo("rename: "+oldName+" = "+newName);
      return;
    }

    for(;;)
    {
      XmlNodePtr xmlChild = config.stack.top().xmlNode->findChild(oldName);
      if(!xmlChild)
        break;
      if(Parallel::isMaster())
        logWarning<<"In '"<<config.currentNodeName()<<"':"<<" config element '"<<oldName<<"' has new name '"<<newName<<"' since "<<time.dateStr()<<Log::endl;
      xmlChild->setName(newName);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void renameDeprecatedChoice(Config &config, std::string &type, const std::string &oldName, const std::string &newName, const Time &time)
{
  try
  {
    if(isCreateSchema(config))
    {
      config.addAppInfo("rename: "+oldName+" = "+newName);
      return;
    }

    if(type == oldName)
    {
      if(Parallel::isMaster())
        logWarning<<"'"<<config.currentNodeName()<<"':' has new choice name '"<<newName<<"' since "<<time.dateStr()<<Log::endl;
      type = newName;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool hasName(Config &config, const std::string &name, Config::Appearance mustSet)
{
  try
  {
    if(isCreateSchema(config))
      return TRUE;

    const Bool found = (config.getChild(name, FALSE/*remove*/) != nullptr);
    if((mustSet == Config::MUSTSET) && !found)
      throw(Exception("config '"+config.currentNodeName()+"' must contain '"+name+"'"));
    return found;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool readConfigSequence(Config &config, const std::string &name, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(isCreateSchema(config))
  {
    config.xscomplexElement(name, Config::SEQUENCE, mustSet, defaultValue, annotation);
    return TRUE;
  }

  XmlNodePtr child = config.getChild(name);
  if(child==nullptr)
  {
    if(mustSet == Config::MUSTSET)
      throw(Exception("config '"+config.currentNodeName()+"' must contain '"+name+"'"));
    return FALSE;
  }
  config.push(child, Config::SEQUENCE);
  return TRUE;
}

/***********************************************/

void endSequence(Config &config)
{
  if(isCreateSchema(config))
  {
    Bool hasChildren = config.stack.top().xmlNode->hasChildren();
    config.pop();
    if(!hasChildren) // without children -> remove complexType
      config.stack.top().xmlLastChild->getChild("xs:complexType");
    return;
  }

  config.notEmptyWarning();
  config.pop();
}

/***********************************************/

Bool readConfigChoice(Config &config, const std::string &name, std::string &choice, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(isCreateSchema(config))
  {
    config.xscomplexElement(name, Config::CHOICE, mustSet, defaultValue, annotation);
    return TRUE;
  }

  XmlNodePtr xmlNode = config.getChild(name);
  if(!xmlNode)
  {
    if(mustSet == Config::MUSTSET)
      throw(Exception("config '"+config.currentNodeName()+"' must contain '"+name+"'"));
    choice = "";
    return FALSE;
  }
  XmlNodePtr child = xmlNode->getNextChild();
  if(!child)
    throw(Exception("config choice element '"+name+"' in '"+config.currentNodeName()+"' has no child"));
  choice = child->getName();

  config.push(xmlNode, Config::CHOICE);
  config.push(child,   Config::CHOICE, config.currentNodeName()+"("+choice+")");
  return TRUE;
}

/***********************************************/

Bool readConfigChoiceElement(Config &config, const std::string &name, const std::string &choice, const std::string &annotation)
{
  if(isCreateSchema(config))
  {
    if(config.stack.top().type == Config::SEQUENCE) // is last choice element a SEQUENCE?
      endSequence(config);
    config.xscomplexElement(name, Config::SEQUENCE, Config::MUSTSET, "", annotation);
    return TRUE;
  }

  return (name == choice);
}

/***********************************************/

void endChoice(Config &config)
{
  if(isCreateSchema(config))
  {
    if(config.stack.top().type == Config::SEQUENCE) // is last choice element a SEQUENCE?
      endSequence(config);
    config.pop();
    return;
  }

  config.notEmptyWarning();
  config.pop();
  config.pop();
}

/***********************************************/
/*** Read Simple Elements **********************/
/***********************************************/

// read Int
template<> Bool readConfig(Config &config, const std::string &name, Int &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  Double v;
  Bool   found = config.getConfigValue(name, "int", mustSet, defaultValue, annotation, v);
  if(found)
    var = static_cast<Int>(round(v));
  return found;
}

/***********************************************/

// read UInt
template<> Bool readConfig(Config &config, const std::string &name, UInt &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  Double v;
  Bool   found = config.getConfigValue(name, "uint", mustSet, defaultValue, annotation, v);
  if(found)
    var = static_cast<UInt>(round(v));
  return found;
}

/***********************************************/

// read Double
template<> Bool readConfig(Config &config, const std::string &name, Double &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  Double v;
  Bool   found = config.getConfigValue(name, "double", mustSet, defaultValue, annotation, v);
  if(found)
    var = v;
  return found;
}

/***********************************************/

// read std::string
template<> Bool readConfig(Config &config, const std::string &name, std::string &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string text;
  Bool found = config.getConfigValue(name, "string", mustSet, defaultValue, annotation, text);
  if(found)
    var = text;
  return found;
}

/***********************************************/

// read Bool
template<> Bool readConfig(Config &config, const std::string &name, Bool &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  Double v;
  Bool   found = config.getConfigValue(name, "boolean", mustSet, defaultValue, annotation, v);
  if(found)
    var = static_cast<Bool>(v);
  return found;
}

/***********************************************/

// read Angle
template<> Bool readConfig(Config &config, const std::string &name, Angle &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  Double v;
  Bool   found = config.getConfigValue(name, "angle", mustSet, defaultValue, annotation, v);
  if(found)
    var = Angle(v*DEG2RAD);
  return found;
}

/***********************************************/

// read Time
template<> Bool readConfig(Config &config, const std::string &name, Time &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  Double v;
  Bool   found = config.getConfigValue(name, "time", mustSet, defaultValue, annotation, v);
  if(found)
    var = mjd2time(v);
  return found;
}

/***********************************************/

// read Doodson
template<> Bool readConfig(Config &config, const std::string &name, Doodson &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string text;
  Bool found = config.getConfigValue(name, "doodson", mustSet, defaultValue, annotation, text);
  if(found)
    var = Doodson(text);
  return found;
}

/***********************************************/

// read FileName
template<> Bool readConfig(Config &config, const std::string &name, FileName &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string text;
  Bool found = config.getConfigValue(name, "filename", mustSet, defaultValue, annotation, text);
  if(found)
    var = FileName(text);
  return found;
}

/***********************************************/

// read Expression
template<> Bool readConfig(Config &config, const std::string &name, ExpressionVariablePtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string text;
  Bool found = config.getConfigValue(name, "expression", mustSet, defaultValue, annotation, text);
  if(found)
    var = ExpressionVariablePtr(new ExpressionVariable(name, text));
  return found;
}

/***********************************************/

// read GnssType
template<> Bool readConfig(Config &config, const std::string &name, GnssType &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string text;
  Bool found = config.getConfigValue(name, "gnssType", mustSet, defaultValue, annotation, text);
  if(found)
    var = GnssType(text);
  return found;
}

/***********************************************/

// read Program
template<> Bool readConfig(Config &config, const std::string &name, ProgramConfig &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  return config.getUnboundedConfigValues(name, "programType", mustSet, defaultValue, annotation, var);
}

/***********************************************/

// read Loop
template<> Bool readConfig(Config &config, const std::string &name, LoopConfig &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  return config.getConfigValue(name, "loopType", mustSet, defaultValue, annotation, var);
}

/***********************************************/
/***********************************************/
