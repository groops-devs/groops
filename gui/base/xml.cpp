/***********************************************/
/**
* @file xml.cpp
*
* @brief Nodes of an XML tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2004-11-14
*/
/***********************************************/

#include <QDebug>
#include <QString>
#include <QTextStream>
#include <QXmlStreamReader>
#include "base/importGroops.h"
#include "base/xml.h"

/***********************************************/

UInt XmlNode::getChildCount(const QStringList &names)
{
  return std::count_if(children.begin(), children.end(), [&names](auto child) {return names.contains(child->getName());});
}

/***********************************************/

XmlNodePtr XmlNode::getChild(const QStringList &names)
{
  auto iter = std::find_if(children.begin(), children.end(), [&names](auto child) {return names.contains(child->getName());});
  if(iter == children.end())
    return XmlNodePtr();
  XmlNodePtr ptr = *iter;
  children.erase(iter);
  return ptr;
}

/***********************************************/

XmlNodePtr XmlNode::getNextChild()
{
  XmlNodePtr ptr(nullptr);
  if(!children.empty())
  {
    ptr = children.front();
    children.pop_front();
  }
  return ptr;
}

/***********************************************/

XmlNodePtr XmlNode::peekNextChild()
{
  if(!children.empty())
    return children.front();
  return nullptr;
}

/***********************************************/

XmlAttrPtr XmlNode::getAttribute(const QString &name)
{
  auto iter = std::find_if(attribute.begin(), attribute.end(), [&name](auto attr) {return attr->name == name;});
  if(iter == attribute.end())
    return XmlAttrPtr();
  XmlAttrPtr ptr = *iter;
  attribute.erase(iter);
  return ptr;
}

/***********************************************/
/***********************************************/

XmlNodePtr XmlNode::parse(QXmlStreamReader &reader, QString &errorMessage)
{
  try
  {
    // header
    auto type = reader.readNext();
    if(type != QXmlStreamReader::StartDocument)
    {
      errorMessage = "XML header expected";
      return nullptr;
    }

    XmlNodePtr             root;
    std::stack<XmlNodePtr> stack;
    for(;;)
    {
      auto type = reader.readNext();
      if(type == QXmlStreamReader::EndDocument)
        break;

      if(type == QXmlStreamReader::Invalid)
      {
        errorMessage = QString("XML parser error (row=%1, col=%2):\n%3").arg(reader.lineNumber()).arg(reader.columnNumber()).arg(reader.errorString());
        return nullptr;
      }

      if(type == QXmlStreamReader::Comment)
      {
        XmlNodePtr xmlNode = XmlNode::create("COMMENT");
        xmlNode->text += reader.text().toString();
        if(!stack.empty())
          stack.top()->addChild(xmlNode);
      }

      if(type == QXmlStreamReader::StartElement)
      {
        XmlNodePtr xmlNode = XmlNode::create(reader.qualifiedName().toString());
        if(stack.empty())
          root = xmlNode;

        for(auto &attribute : reader.attributes())
        {
          XmlAttrPtr attr(new XmlAttribute());
          attr->name = attribute.name().toString();
          attr->text = attribute.value().toString();
          xmlNode->addAttribute(attr);
        }

        if(!stack.empty())
          stack.top()->addChild(xmlNode);
        stack.push(xmlNode);
      } // if(StartElement)

      if(type == QXmlStreamReader::Characters)
      {
        if(!reader.isWhitespace())
          stack.top()->text += reader.text().toString();
      }

      if(type == QXmlStreamReader::EndElement)
      {
        XmlNodePtr xmlNode = stack.top();
        if(xmlNode->name != reader.qualifiedName().toString())
        {
          errorMessage = QString("XML parser error: (row=%1, col=%2):\nunexpected end element").arg(reader.lineNumber()).arg(reader.columnNumber());
          return nullptr;
        }
        stack.pop();
      }  // if(EndElement)
    } // for(;;)

    return root;
  }
  catch(std::exception &e)
  {
    errorMessage = QString("XML parser error: critical error");
    return nullptr;
  }
}

/***********************************************/

XmlNodePtr XmlNode::read(const QByteArray &data, QString &errorMessage)
{
  QXmlStreamReader reader(data);
  return parse(reader, errorMessage);
}

/***********************************************/

XmlNodePtr XmlNode::read(QIODevice *device, QString &errorMessage)
{
  QXmlStreamReader reader(device);
  return parse(reader, errorMessage);
}

/***********************************************/
/***********************************************/

void XmlNode::writeFile(const QString &fileName, const XmlNodePtr &root)
{
  try
  {
    QFile file(fileName);
    if(file.open(QFile::WriteOnly | QFile::Truncate))
    {
      QTextStream stream(&file);
      write(stream, root, false);
    }
    else
      throw(Exception("cannot open file"));
  }
  catch(std::exception &e)
  {
    throw(Exception(_GROOPS_ERRORLINE+" filename='"+fileName.toStdString()+"'\n"+e.what()));
  }
}

/***********************************************/

void XmlNode::write(QTextStream &stream, const XmlNodePtr &root, bool writeCommentsAsElements)
{
  try
  {
    stream<<"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    root->write(stream, 0, writeCommentsAsElements);
  }
  catch(std::exception &e)
  {
    throw(Exception(std::string("In XmlNode::write:\n")+e.what()));
  }
}

/***********************************************/

void XmlNode::write(QTextStream &stream, UInt depth, bool writeCommentsAsElements)
{
  auto sanitizeXML = [](const QString &in)
  {
    if(in.isEmpty())
      return QString();
    std::vector<QString> chars = {"&", ">", "<", "\"", "'"};
    std::vector<QString> replacements = {"&amp;", "&gt;", "&lt;", "&quot;", "&apos;"};
    QString out(in);
    for(UInt k=0; k<chars.size(); k++)
      out = out.replace(chars.at(k), replacements.at(k));
    return out;
  };

  // special node -> translate to XML comment
  if(!writeCommentsAsElements && (getName() == "COMMENT"))
  {
    stream<<QString(depth, '\t')<<"<!--"<<sanitizeXML(getText())<<"-->\n";
    return;
  }

  // start tag
  stream<<QString(depth, '\t')<<"<"<<getName();

  // attribute
  for(auto attr : attribute)
    stream<<" "<<attr->name<<"=\""<<sanitizeXML(attr->text)<<"\"";

  // short form?
  if(getText().isEmpty() && children.empty())
  {
    stream<<"/>\n";
    return;
  }
  stream<<">"<<sanitizeXML(getText());

  // children
  if(!children.empty())
  {
    stream<<"\n";
    for(auto child : children)
      child->write(stream, depth+1, writeCommentsAsElements);
    stream<<QString(depth, '\t');
  }

  // end tag
  stream<<"</"<<getName()<<">\n";
}

/***********************************************/
/***********************************************/

XmlAttrPtr readAttribute(const XmlNodePtr &node, const QString &name, QString &var)
{
  try
  {
    XmlAttrPtr attr = node->getAttribute(name);
    if(attr)
      var = attr->text;
    return attr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlAttrPtr writeAttribute(const XmlNodePtr &node, const QString &name, const QString &var)
{
  XmlAttrPtr attr = XmlAttrPtr(new XmlAttribute(name, var));
  node->addAttribute(attr);
  return attr;
}

/***********************************************/
