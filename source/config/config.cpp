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
  VariableList varList;
  std::map<std::string, XmlNodePtr> links;
  if(stack.size())
  {
    varList = stack.top().varList;
    links   = stack.top().links;
  }
  std::string name_ = (name.empty()) ? Config::currentNodeName()+"."+xmlNode->getName() : name;
  stack.emplace(xmlNode, type, name_);
  stack.top().varList = varList;
  stack.top().links   = links;
}

/***********************************************/

void Config::pop()
{
  if(stack.top().loopPtr) // finish old loop
    throw(Exception("exit node with unfinished loop"));
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

    push(root, SEQUENCE, "groops");

    // local variables before global?
    // ------------------------------
    XmlNodePtr xmlNode;
    while((xmlNode = root->findNextChild()))
    {
      XmlAttrPtr label = xmlNode->getAttribute("label");
      if(!label)                           // is not local variable?
        break;
      xmlNode = root->getNextChild();      // remove
      xmlNode->setName(label->getText());  // replace typename with label
      stack.top().links[xmlNode->getName()] = xmlNode;
    }

    // global: replace typename with label
    // -----------------------------------
    XmlNodePtr global = root->getChild("global");
    if(global)
      for(XmlNodePtr &child : global->getChildren())
      {
        XmlAttrPtr label = child->getAttribute("label");
        if(!label)
          throw(Exception("elements in global must have attribute 'label'"));
        child->setName(label->getText());
        stack.top().links[child->getName()] = child;
      }

    // add global elements from command line
    // -------------------------------------
    for(const auto &command : commandlineGlobals)
    {
      XmlNodePtr xmlNode = XmlNode::create(command.first);
      xmlNode->setText(command.second);
      stack.top().links[command.first] = xmlNode;
    }

    // set global variables
    // --------------------
    for(auto &pair : stack.top().links)
      if(!pair.second->hasChildren())  // not complex type?
      {
        XmlAttrPtr link = pair.second->findAttribute("link");
        stack.top().varList.setVariable(pair.second->getName(), ((link) ? "{"+link->getText()+"}" : pair.second->getText()));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Config::hasName(const std::string &name)
{
  for(;;)
  {
    if(stack.top().xmlLastChild)
    {
      if(stack.top().xmlLastChild->getName() != name)
        throw(Exception("loop error"));
      return TRUE;
    }

    if(!stack.top().xmlNode->findChild(name))
      return FALSE;

    // local variables before?
    // -----------------------
    XmlNodePtr xmlNode;
    while((xmlNode = stack.top().xmlNode->getNextChild()))
    {
      XmlAttrPtr label = xmlNode->getAttribute("label");
      if(label) // xmlNode is local variable
      {
        xmlNode->setName(label->getText());  // replace typename with label
        stack.top().links[xmlNode->getName()] = xmlNode;
        if(!xmlNode->hasChildren())          // not complex type? -> variable
        {
          XmlAttrPtr link = xmlNode->findAttribute("link");
          stack.top().varList.setVariable(xmlNode->getName(), ((link) ? "{"+link->getText()+"}" : xmlNode->getText()));
        }
        continue;
      }
      if(xmlNode->getName() == name)
        break;
      // unknown element: move to end
      stack.top().xmlNode->addChild(xmlNode);
    }

    if(xmlNode == nullptr)
      return FALSE;

    // -----------------------
    std::function<void(XmlNodePtr&)> resolveLink = [&](XmlNodePtr &xmlNode)
    {

      XmlAttrPtr link = xmlNode->getAttribute("link");
      if(link)
      {
        XmlNodePtr xmlNodeLink = stack.top().links[link->getText()];
        if(!xmlNodeLink || xmlNodeLink->getAttribute("resolving"))
          throw(Exception(std::string("cannot resolve link: ")+xmlNode->getName()+" -> "+link->getText()+"'"));


        // add temporary resolving attribute to prevent link loops
        writeAttribute(xmlNodeLink, "resolving", 1);

        XmlNodePtr xmlNodeNew = xmlNodeLink->clone();
        resolveLink(xmlNodeNew);
        xmlNodeNew->setName(xmlNode->getName());
        while(xmlNode->hasAttribute())
          xmlNodeNew->addAttribute(xmlNode->getNextAttribute());
        xmlNode = xmlNodeNew;

        // remove attributes as link was resolved successfully
        xmlNodeLink->getAttribute("resolving");
      }
    };
    // -----------------------

    resolveLink(xmlNode);

    // loop?
    // -----
    if(stack.top().loopPtr)
    {
      if(!stack.top().loopPtr->iteration(stack.top().varList)) // finish loop
      {
        stack.top().varList = stack.top().loopVarListOld; // restore old varList
        stack.top().loopPtr = nullptr;
        continue;
      }
      stack.top().xmlNode->prependChild(xmlNode->clone());  // restore for next loop
    }
    else
    {
      // expand loop
      XmlAttrPtr attr = xmlNode->getAttribute("loop");
      if(attr)
      {
        XmlNodePtr xmlNodeLoop = stack.top().links[attr->getText()];
        if(!xmlNodeLoop)
          throw(Exception(std::string("cannot resolve loop link: ")+xmlNode->getName()+" -> '"+attr->getText()+"'"));
        XmlNodePtr tmp = XmlNode::create("tmp");
        tmp->addChild(xmlNodeLoop->clone());  // make copy
        push(tmp, SEQUENCE, currentNodeName());
        LoopPtr loopPtr;
        readConfig(*this, xmlNodeLoop->getName(), loopPtr, Config::MUSTSET, "", "");
        pop();
        stack.top().loopPtr = loopPtr;
        stack.top().loopVarListOld = stack.top().varList;
        stack.top().xmlNode->prependChild(xmlNode);  // restore for next loop
        continue;  // start loop
      }
    }

    // check condition
    XmlAttrPtr attr = xmlNode->getAttribute("condition");
    if(attr)
    {
      // expand condition
      XmlNodePtr xmlNodeCondition = stack.top().links[attr->getText()];
      if(!xmlNodeCondition)
        throw(Exception(std::string("cannot resolve condition link: '")+xmlNode->getName()+"' -> '"+attr->getText()+"'"));
      attr = nullptr;
      XmlNodePtr tmp = XmlNode::create("tmp");
      tmp->addChild(xmlNodeCondition->clone()); // make copy
      push(tmp, SEQUENCE, currentNodeName());
      ConditionPtr condition;
      readConfig(*this, xmlNodeCondition->getName(), condition, Config::MUSTSET, "", "");
      pop();

      // check condition
      if(!condition->condition(stack.top().varList))
        continue; // node disabled -> try next element
    }

    stack.top().xmlLastChild = xmlNode;
    return TRUE;
  } // for(;;)
}

/***********************************************/

XmlNodePtr Config::getChild(const std::string &name)
{
  try
  {
    hasName(name); // get current node in xmlLastChild
    XmlNodePtr xmlNode;
    std::swap(stack.top().xmlLastChild, xmlNode);
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Config::getConfigText(const std::string &name, const std::string &type, Config::Appearance mustSet,
                           const std::string &defaultValue, const std::string &annotation, Bool parse, std::string &text)
{
  if(createSchema)
  {
    if((mustSet == Config::DEFAULT) && defaultValue.empty())
      throw(Exception("In readConfig("+currentNodeName()+"."+name+", type="+type+"): Config::DEFAULT with empty defaultValue"));
    xselement(name, type, mustSet, ONCE, defaultValue, annotation);
    return FALSE;
  }

  XmlNodePtr child = getChild(name);
  if(child)
  {
    if(child->hasChildren())
      throw(Exception(currentNodeName()+"."+child->getName()+" with unexpected children"));
    text = child->getText();
  }

  if(text.empty() && (mustSet == MUSTSET))
    throw(Exception("config element '"+currentNodeName()+"' must contain '"+name+"'"));
  if(text.empty() && (mustSet == DEFAULT))
    text = defaultValue;

  if(parse)
  {
    Bool resolved = TRUE;
    text = StringParser::parse(name, text, stack.top().varList, resolved);
    if(!resolved)
      throw(Exception("In readConfig("+currentNodeName()+"."+name+", type="+type+")='"+text+"': unresolved variables"));
  }
  return !text.empty();
}

/***********************************************/

Bool Config::getConfigValue(const std::string &name, const std::string &type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation, Double &v)
{
  std::string text;
  try
  {
    if(!getConfigText(name, type, mustSet, defaultValue, annotation, TRUE, text))
      return FALSE;
    v = ExpressionVariable::parse(text, stack.top().varList);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("'"+name+"' = '"+text+"'", e)
  }
}

/***********************************************/

Bool Config::getConfig(const std::string &name, Config::Appearance mustSet, Config &conf)
{
  try
  {
    XmlNodePtr xmlNode = XmlNode::create(name);
    XmlNodePtr child   = getChild(name);
    if(child)
      xmlNode->addChild(child);
    else if(mustSet == MUSTSET)
      throw(Exception("config element '"+currentNodeName()+"' must contain '"+name+"'"));

    // make copy
    conf.createSchema = FALSE;
    conf.push(xmlNode, Config::SEQUENCE, currentNodeName());
    conf.stack.top().links   = stack.top().links;
    conf.stack.top().varList = stack.top().varList;

    return child != nullptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("'"+name+"'", e)
  }
}

/***********************************************/

Bool Config::getUnboundedConfig(const std::string &name, Config &conf)
{
  try
  {
    if(stack.top().xmlLastChild || stack.top().loopPtr)
      throw(Exception("last element not yet finished processed"));

    XmlNodePtr xmlNode = XmlNode::create(name);
    XmlNodePtr child;
    while((child = stack.top().xmlNode->findNextChild()))
    {
      if((child->getName() == name) || child->findAttribute("label")) // child is element or local variable
        xmlNode->addChild(stack.top().xmlNode->getNextChild());
      else if(!stack.top().xmlNode->findChild(name))
        break;
    }

    // make copy
    conf.createSchema = FALSE;
    conf.push(xmlNode, Config::SEQUENCE, currentNodeName());
    conf.stack.top().links   = stack.top().links;
    conf.stack.top().varList = stack.top().varList;

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
    if(stack.top().xmlNode->hasChildren())
    {
      logWarningOnce<<"*** Warning: unknown variables in '"<<currentNodeName()<<"' (ignored):"<<Log::endl;
      while(stack.top().xmlNode->hasChildren())
        logWarningOnce<<"  name = "<<stack.top().xmlNode->getNextChild()->getName()<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Config::copy(Config &config, const VariableList &variableList) const
{
  try
  {
    config.createSchema = FALSE;
    config.push(stack.top().xmlNode->clone(), Config::SEQUENCE, currentNodeName());
    config.stack.top().links    = stack.top().links;
    config.stack.top().varList  = stack.top().varList;
    config.stack.top().varList += variableList;
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
      program->run(config, Parallel::selfCommunicator());
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

void ProgramConfig::run(VariableList &variableList, Parallel::CommunicatorPtr comm) const
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
          StackNode top = config.stack.top();
          config.stack.pop(); // coment is given in <program> not in <choiceElement>
          XmlAttrPtr attr = config.stack.top().xmlNode->getAttribute("comment");
          if(attr)
            comment = attr->getText();
          config.stack.push(top);

          Parallel::barrier(comm);
          if(comment.empty())
            logStatus<<"--- "<<program->name()<<" ---"<<Log::endl;
          else
          {
            try
            {
              Bool resolved = TRUE;
              comment = StringParser::parse("comment", comment, config.getVarList(), resolved);
            }
            catch(std::exception &/*e*/) {}
            logStatus<<"--- "<<program->name()<<" ("<<comment<<") ---"<<Log::endl;
          }
          program->run(config, comm);
          Parallel::barrier(comm);
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
      logWarningOnce<<"In '"<<config.currentNodeName()<<"':"<<" config element '"<<oldName<<"' has new name '"<<newName<<"' since "<<time.dateStr()<<Log::endl;
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
      logWarningOnce<<"'"<<config.currentNodeName()<<"':' has new choice name '"<<newName<<"' since "<<time.dateStr()<<Log::endl;
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

    const Bool found = config.hasName(name);
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
  Bool found = config.getConfigText(name, "string", mustSet, defaultValue, annotation, TRUE, text);
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
  Bool found = config.getConfigText(name, "doodson", mustSet, defaultValue, annotation, TRUE, text);
  if(found)
    var = Doodson(text);
  return found;
}

/***********************************************/

// read FileName
template<> Bool readConfig(Config &config, const std::string &name, FileName &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string text;
  Bool found = config.getConfigText(name, "filename", mustSet, defaultValue, annotation, FALSE, text);
  if(found)
  {
    Bool resolved = TRUE;
    std::string textResolved = StringParser::parse(name, text, VariableList(), resolved);
    if(resolved)
      var = FileName(textResolved);
    else
      var = FileName(text, config.getVarList());
  }
  return found;
}

/***********************************************/

// read Expression
template<> Bool readConfig(Config &config, const std::string &name, ExpressionVariablePtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string text;
  Bool found = config.getConfigText(name, "expression", mustSet, defaultValue, annotation, FALSE, text);
  if(found)
    var = std::make_shared<ExpressionVariable>(name, text, config.getVarList());
  return found;
}

/***********************************************/

// read GnssType
template<> Bool readConfig(Config &config, const std::string &name, GnssType &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string text;
  Bool found = config.getConfigText(name, "gnssType", mustSet, defaultValue, annotation, TRUE, text);
  if(found)
    var = GnssType(text);
  return found;
}

/***********************************************/

// read Program
template<> Bool readConfig(Config &config, const std::string &name, ProgramConfig &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(isCreateSchema(config))
  {
    config.xselement(name, "programType", mustSet, Config::UNBOUNDED, defaultValue, annotation);
    return FALSE;
  }

  return config.getUnboundedConfig(name, var);
}

/***********************************************/
/***********************************************/
