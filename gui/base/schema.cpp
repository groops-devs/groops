/***********************************************/
/**
* @file schema.cpp
*
* @brief XSD (XML schema) tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#include <QTextStream>
#include "base/importGroops.h"
#include "base/xml.h"
#include "schema.h"

/***********************************************/

XsdElementPtr XsdComplex::getXsdElement(QString name) const
{
  try
  {
    XsdElementPtr xsdElement(nullptr);
    while(!xsdElement)
    {
      auto iter = std::find_if(element.begin(), element.end(), [&](XsdElementPtr elem){ return elem->name == name; });
      if(iter != element.end())
      {
        xsdElement = *iter;
        break;
      }

      auto iterRenames = std::find_if(renames.begin(), renames.end(), [&](auto pair){ return pair.first == name; });
      if(iterRenames == renames.end())
        break;
      name = iterRenames->second;
    }

    return xsdElement;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Schema::Schema(const QString &fileName)
{
  try
  {
    XmlNodePtr xmlNode = XmlNode::readFile(fileName);

    // read ComplexTypes
    UInt count = childCount(xmlNode, QString("xs:complexType"));
    for(UInt i=0; i<count; i++)
    {
      XmlNodePtr xmlNode2 = getChild(xmlNode, QString("xs:complexType"));
      XsdElementPtr element(new XsdElement());
      complexType.push_back(element);

      // read name
      XmlAttrPtr attr = xmlNode2->getAttribute("name");
      if(attr==nullptr)
        throw(Exception("xs:complexType without name"));
      element->name = attr->getText();
      element->complex = readComplex(xmlNode2);
    }

    // read elements
    rootElement = readElement(getChild(xmlNode,QString("xs:element")));
    setComplexPointer(rootElement); // set recursively pointer to the complex types
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XsdElementPtr Schema::readElement(XmlNodePtr xmlNode)
{
  try
  {
    if(xmlNode==nullptr)
      throw(Exception("xs:element is null"));

    XsdElementPtr element(new XsdElement());
    element->complex = nullptr;

    // read name
    XmlAttrPtr attr = xmlNode->getAttribute("name");
    if(attr==nullptr)
      throw(Exception("xs:element without name"));
    element->name = attr->getText();

    // annotation, tags, and renames
    std::map<QString, QString> renames;
    XmlNodePtr child = getChild(xmlNode, "xs:annotation");
    if(child!=nullptr)
    {
      XmlNodePtr doc = getChild(child, "xs:documentation");
      if(doc!=nullptr)
        element->annotation = doc->getText();

      XmlNodePtr info = getChild(child, "xs:appinfo");
      while(info!=nullptr)
      {
        if(info->getText().startsWith("tag:"))
          element->tags.append(info->getText().split(" ", QString::SkipEmptyParts)[1]);
        else if(info->getText().startsWith("rename:")) // format: "rename: oldName = newName"
        {
          QStringList splits = info->getText().split(" ", QString::SkipEmptyParts);
          renames[splits[1]] = splits[3];
        }
        info = getChild(child, "xs:appinfo");
      }
    }

    // must set?
    element->optional = false;
    attr = xmlNode->getAttribute("minOccurs");
    if(attr!=nullptr)
    {
      UInt minOccur; attr->getValue(minOccur);
      if(minOccur==0)
        element->optional = true;
    }

    // unbounded
    element->unbounded = false;
    attr = xmlNode->getAttribute("maxOccurs");
    if((attr!=nullptr)&&(attr->getText()=="unbounded"))
      element->unbounded = true;

    // default value
    attr = xmlNode->getAttribute("default");
    if(attr!=nullptr)
      element->defaultValue = attr->getText();

    // type
    attr = xmlNode->getAttribute("type");
    if(attr!=nullptr)
      element->type = attr->getText();

    XmlNodePtr xmlNode2 = getChild(xmlNode, "xs:complexType");
    if(xmlNode2 != nullptr)
    {
      element->complex = readComplex(xmlNode2);
      // create unique identifier
      QTextStream ss(&element->type);
      ss<<element->complex->type<<"Start-";
      for(UInt i=0; i<element->complex->element.size(); i++)
        ss<<element->complex->element.at(i)->type<<"-";
      ss<<element->complex->type<<"End";
      element->complex->renames.insert(renames.begin(), renames.end());
    }

    return element;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XsdComplexPtr Schema::readComplex(XmlNodePtr xmlNode)
{
  try
  {
    XsdComplexPtr complex(new XsdComplex());
    complex->isReady = false;

    // read renames
    XmlNodePtr xmlNode2 = getChild(xmlNode, "xs:annotation");
    if(xmlNode2!=nullptr)
    {
      XmlNodePtr info = getChild(xmlNode2, "xs:appinfo");
      while(info!=nullptr)
      {
        if(info->getText().startsWith("rename:")) // format: "rename: oldName = newName"
        {
          QStringList splits = info->getText().split(" ", QString::SkipEmptyParts);
          complex->renames[splits[1]] = splits[3];
        }
        info = getChild(xmlNode2, "xs:appinfo");
      }
    }

    complex->type = "sequence";
    xmlNode2 = getChild(xmlNode, "xs:sequence");
    if(xmlNode2 == nullptr)
    {
      xmlNode2 = getChild(xmlNode, "xs:choice");
      complex->type = "choice";
    }
    if (xmlNode2 == nullptr)
      throw(Exception("xs:complex must contain xs:sequence or xs:choice"));

     // read children
    UInt count = childCount(xmlNode2, "xs:element", true);
    complex->element.resize(count);
    for(UInt i=0; i<count; i++)
      complex->element.at(i) = readElement(getChild(xmlNode2, "xs:element"));

    complexList.push_back(complex);
    return complex;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Schema::setComplexPointer(XsdElementPtr element)
{
  try
  {
    if(!element->complex)
      for(UInt i=0; i<complexType.size(); i++)
        if(element->type == complexType.at(i)->name)
        {
          element->complex = complexType.at(i)->complex;
          break;
        }

    if((element->complex) && (!element->complex->isReady))
    {
      element->complex->isReady = true;
      for(UInt i=0; i<element->complex->element.size(); i++)
        setComplexPointer(element->complex->element.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

std::vector<XsdElementPtr> Schema::programList() const
{
  std::vector<XsdElementPtr> programList;

  for(auto element : complexType)
    if(element->name == "programmeType" || element->name == "programType")
      for(auto program : element->complex->element)
        programList.push_back(program);

  return programList;
}

/***********************************************/

QString Schema::programType() const
{

  for(auto element : complexType)
    if(element->name == "programmeType" || element->name == "programType")
      return element->name;

  return QString();
}

/***********************************************/

bool Schema::validateSchema(QString fileName)
{
  try
  {
    if(fileName.isEmpty())
      return false;
    Schema schema(fileName);
  }
  catch(...)
  {
    return false;
  }
  return true;
}

/***********************************************/
/***********************************************/

