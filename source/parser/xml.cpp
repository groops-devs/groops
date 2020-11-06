/***********************************************/
/**
* @file xml.cpp
*
* @brief Nodes of a XML tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2004-11-14
*
*/
/***********************************************/

#include <expat.h>
#include "base/importStd.h"
#include "base/string.h"
#include "base/angle.h"
#include "base/time.h"
#include "xml.h"

/***********************************************/

XmlNodePtr XmlNode::clone() const
{
  XmlNodePtr ptr = create(name_);
  ptr->text_ = text_;
  for(auto attr : attribute)
    ptr->addAttribute(std::make_shared<XmlAttribute>(*attr));
  for(auto child : children)
    ptr->addChild(child->clone());
  return ptr;

}

/***********************************************/

UInt XmlNode::getChildCount(const std::string &name)
{
  return std::count_if(children.begin(), children.end(), [&name](auto child) {return child->getName() == name;});
}

/***********************************************/

XmlNodePtr XmlNode::getChild(const std::string &name)
{
  auto iter = std::find_if(children.begin(), children.end(), [&name](auto child) {return child->getName() == name;});
  if(iter == children.end())
    return XmlNodePtr();
  XmlNodePtr ptr = *iter;
  children.erase(iter);
  return ptr;
}

/***********************************************/

XmlNodePtr XmlNode::findChild(const std::string &name)
{
  auto iter = std::find_if(children.begin(), children.end(), [&name](auto child) {return child->getName() == name;});
  return (iter != children.end()) ? *iter : XmlNodePtr();
}

/***********************************************/

XmlNodePtr XmlNode::getNextChild()
{
  XmlNodePtr ptr;
  if(!children.empty())
  {
    ptr = children.front();
    children.pop_front();
  }
  return ptr;
}

/***********************************************/
/***********************************************/

XmlAttrPtr XmlNode::getAttribute(const std::string &name)
{
  auto iter = std::find_if(attribute.begin(), attribute.end(), [&name](auto attr) {return attr->getName() == name;});
  if(iter == attribute.end())
    return XmlAttrPtr();
  XmlAttrPtr ptr = *iter;
  attribute.erase(iter);
  return ptr;
}

/***********************************************/

XmlAttrPtr XmlNode::findAttribute(const std::string &name)
{
  auto iter = std::find_if(attribute.begin(), attribute.end(), [&name](auto attr) {return attr->getName() == name;});
  return (iter != attribute.end()) ? *iter : XmlAttrPtr();
}


/***********************************************/

XmlAttrPtr XmlNode::getNextAttribute()
{
  XmlAttrPtr ptr;
  if(!attribute.empty())
  {
    ptr = attribute.front();
    attribute.pop_front();
  }
  return ptr;
}

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

class XmlReadFile
{
public:
  XmlNodePtr    root;
  std::stack<XmlNodePtr> stack; // XML-Knoten, die noch nicht vollstaendig bearbeitet wurden
};

/***********************************************/

// StartHandler
static void XmlNodeStartElement(XmlReadFile *file, const XML_Char *name, const XML_Char **atts)
{
  XmlNodePtr ptr = XmlNode::create(name);
  if(file->stack.empty())
    file->root = ptr;

  // Attributes
  Bool disabled = FALSE;
  for(UInt i=0; atts[i]!=nullptr; i+=2)
  {
    if(atts[i]==std::string("disabled"))
    {
      if(atts[i+1] == std::string("1"))
        disabled = TRUE;
    }
    else
    {
      XmlAttrPtr attr(new XmlAttribute());
      attr->name = atts[i];
      attr->text = atts[i+1];
      ptr->addAttribute(attr);
    }
  }

  if((!file->stack.empty()) && (!disabled))
    file->stack.top()->addChild(ptr);
  file->stack.push(ptr);

}

/***********************************************/

// TextHandler
static void XmlNodeTextElement(XmlReadFile *file, const XML_Char *s, int len)
{
  file->stack.top()->addText(std::string(s, len));
}

/***********************************************/

// EndHandler
static void XmlNodeEndElement(XmlReadFile *file, const XML_Char */*name*/)
{
  XmlNodePtr ptr = file->stack.top();
  file->stack.top() = XmlNodePtr();
  file->stack.pop();

  // Remove white spaces at the beginning of text.
  std::size_t pos = ptr->getText().find_first_not_of(" \n\t");
  if(pos!=std::string::npos)
    ptr->setText(ptr->getText().substr(pos));
  else
    ptr->setText(std::string());
}

/***********************************************/

XmlNodePtr XmlNode::read(std::istream &stream)
{
  try
  {
    stream.exceptions(std::ios::badbit|std::ios::failbit);

    // Init Parser
    XML_Parser expatParser = XML_ParserCreate(nullptr);

    // Set Handler
    std::shared_ptr<XmlReadFile> file(new XmlReadFile);
    XML_SetUserData(expatParser, file.get());
    XML_SetElementHandler(expatParser, reinterpret_cast<XML_StartElementHandler>(XmlNodeStartElement), reinterpret_cast<XML_EndElementHandler>(XmlNodeEndElement));
    XML_SetCharacterDataHandler(expatParser, reinterpret_cast<XML_CharacterDataHandler>(XmlNodeTextElement));

    static const UInt BUFFERSIZE = 1024*1024;
    Char *buffer = new Char[BUFFERSIZE];

    do
    {
      stream.get(buffer, BUFFERSIZE, '\0');
      auto bufferCount = stream.gcount();
      XML_Status status = XML_Parse(expatParser, buffer, bufferCount, stream.eof());

      if(status!=XML_STATUS_OK)
      {
        std::stringstream ss;
        ss<<"Parser Error in row "<<XML_GetCurrentLineNumber(expatParser);
        ss<<", column "<<XML_GetCurrentColumnNumber(expatParser);
        ss<<": "<<XML_ErrorString(XML_GetErrorCode(expatParser));
        throw(Exception(ss.str()));
      }
    }
    while(!stream.eof());

    delete[] buffer;
    XML_ParserFree(expatParser);

    return file->root;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void XmlNode::write(std::ostream &stream, UInt depth)
{
  // start tag
  stream<<std::string(depth, '\t')<<"<"<<getName();

  // Attribute
  for(auto attr : attribute)
    stream<<" "<<attr->getName()<<"=\""<<sanitizeXML(attr->getText())<<"\"";

  // short form?
  if(getText().empty() && children.empty())
  {
    stream<<"/>"<<std::endl;;
    return;
  }

  stream<<">"<<sanitizeXML(getText());

  // children
  if(!children.empty())
  {
    stream<<std::endl;
    for(auto child : children)
      child->write(stream, depth+1);
    stream<<std::string(depth, '\t');
  }

  // end tag
  stream<<"</"<<getName()<<">"<<std::endl;
}

/***********************************************/

void XmlNode::write(std::ostream &stream, const XmlNodePtr &root)
{
  try
  {
    stream.exceptions(std::ios::badbit|std::ios::failbit);
    stream<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>"<<std::endl;
    root->write(stream);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string XmlNode::sanitizeXML(const std::string &str)
{
  return String::replaceAll(str, {{"&", "&amp;"}, {">", "&gt;"}, {"<", "&lt;"}, {"'", "&apos;"}, {"\"", "&quot;"}});
}

/***********************************************/
/***** FUNCTIONS *******************************/
/***********************************************/

XmlNodePtr getNextChild(const XmlNodePtr &ptr, std::string &name)
{
  XmlNodePtr child = ptr->getNextChild();
  if(child)
    name = child->getName();
  return child;
}

/***********************************************/

UInt childCount(const XmlNodePtr &parent, const std::string &name, Bool mustGreaterZero)
{
  UInt count = parent->getChildCount(name);
  if((count==0) && mustGreaterZero)
    throw(Exception("'"+parent->getName()+"' must contain '"+name+"'"));
  return count;
}

/***********************************************/
