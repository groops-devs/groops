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
#include <QFile>
#include "base/importGroops.h"
#include "base/xml.h"
#include "schema.h"

/***********************************************/

XsdElementPtr XsdComplex::getXsdElement(QString name) const
{
  try
  {
    auto iter = std::find_if(elements.begin(), elements.end(), [&](XsdElementPtr element) {return element->names.contains(name);});
    if(iter == elements.end())
      return nullptr;
    return *iter;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool Schema::readFile(QString fileName)
{
  try
  {
    complexType.clear();
    rootElement = nullptr;
    if(fileName.isEmpty())
      return false;

    QFile file(fileName);
    if(!file.open(QFile::ReadOnly | QFile::Text))
      return false;
    QString errorMessage;
    XmlNodePtr xmlNode = XmlNode::read(&file, errorMessage);
    if(!xmlNode)
      return false;

    // read ComplexTypes
    UInt count = xmlNode->getChildCount(QString("xs:complexType"));
    for(UInt i=0; i<count; i++)
      complexType.push_back(readElement(xmlNode->getChild(QString("xs:complexType")), {}, true));
    // read elements
    rootElement = readElement(xmlNode->getChild(QString("xs:element")), {}, false);

    setComplexPointer(rootElement); // set recursively pointer to the complex types
    return true;
  }
  catch(std::exception &e)
  {
    complexType.clear();
    rootElement = nullptr;
    return false;
  }
}

/***********************************************/

XsdElementPtr Schema::readElement(XmlNodePtr xmlNode, const std::map<QString, QString> &renames, bool isComplexType)
{
  try
  {
    if(!xmlNode)
      throw(Exception("xs:element is null"));

    XsdElementPtr element(new XsdElement());
    element->complex = nullptr;

    // read name
    XmlAttrPtr attr = xmlNode->getAttribute("name");
    if(!attr)
      throw(Exception("xs:element without name"));
    element->names = QStringList(attr->text);
    // renames
    auto findRename = [&](const QString name){return std::find_if(renames.begin(), renames.end(), [&](auto pair){return pair.second == name;});};
    auto iterRenames = findRename(element->names.back());
    while(iterRenames != renames.end())
    {
      element->names.push_back(iterRenames->first);
      iterRenames = findRename(element->names.back());
    }

    // annotation, tags, and renames
    std::map<QString, QString> childRenames;
    XmlNodePtr child = xmlNode->getChild("xs:annotation");
    if(child)
    {
      XmlNodePtr doc = child->getChild("xs:documentation");
      if(doc)
        element->annotation = doc->getText();

      for(XmlNodePtr info = child->getChild("xs:appinfo"); info; info = child->getChild("xs:appinfo"))
      {
        if(info->getText().startsWith("rename:")) // format: "rename: oldName = newName"
        {
          QStringList splits = info->getText().split(" ", Qt::SkipEmptyParts);
          childRenames[splits[1]] = splits[3];
        }
        else if(info->getText().startsWith("tag:"))
          element->tags.append(info->getText().split(" ", Qt::SkipEmptyParts)[1]);
      }
    }

    // attributes
    QString minOccurs, maxOccurs;
    readAttribute (xmlNode, "type",      element->type);
    readAttribute (xmlNode, "default",   element->defaultValue);
    readAttribute (xmlNode, "minOccurs", minOccurs);
    readAttribute (xmlNode, "maxOccurs", maxOccurs);
    element->optional  = (minOccurs == "0");
    element->unbounded = (maxOccurs == "unbounded");

    XmlNodePtr xmlNode2;
    if(isComplexType)
      xmlNode2 = xmlNode;
    else
      xmlNode2 = xmlNode->getChild("xs:complexType");
    if(xmlNode2)
    {
      element->complex = readComplex(xmlNode2, childRenames);
      if(element->type.isEmpty()) // create unique identifier
      {
        QTextStream ss(&element->type);
        ss<<element->complex->type<<"Start-";
        for(auto child : element->complex->elements)
          ss<<child->type<<"-";
        ss<<element->complex->type<<"End";
      }
    }

    return element;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XsdComplexPtr Schema::readComplex(XmlNodePtr xmlNode, const std::map<QString, QString> &renames)
{
  try
  {
    XsdComplexPtr complex(new XsdComplex());
    complex->isReady = false;

    complex->type = "sequence";
    XmlNodePtr xmlNode2 = xmlNode->getChild("xs:sequence");
    if(!xmlNode2)
    {
      xmlNode2 = xmlNode->getChild("xs:choice");
      complex->type = "choice";
    }
    if(!xmlNode2)
      throw(Exception("xs:complex must contain xs:sequence or xs:choice"));

     // read children
    UInt count = xmlNode2->getChildCount("xs:element");
    if(!count)
      throw(Exception("xs:sequence/xs::choice must contain xs:element"));
    complex->elements.resize(count);
    for(UInt i=0; i<count; i++)
      complex->elements.at(i) = readElement(xmlNode2->getChild("xs:element"), renames, false);

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
    {
      auto iter = std::find_if(complexType.begin(), complexType.end(), [&](const auto &child) {return child->names.front() == element->type;});
      if(iter != complexType.end())
        element->complex = (*iter)->complex;
    }

    if((element->complex) && (!element->complex->isReady))
    {
      element->complex->isReady = true;
      for(auto &child : element->complex->elements)
        setComplexPointer(child);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

