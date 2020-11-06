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

#include <QString>
#include <QTextStream>
#include <QFile>
#include <QDomNode>
#include "base/importGroops.h"
#include "base/xml.h"

/***********************************************/

XmlNodePtr XmlNode::clone() const
{
  XmlNodePtr ptr = create(name);
  ptr->text  = text;
  for(std::list<XmlAttrPtr>::const_iterator im=attribute.begin(); im!=attribute.end(); im++)
    ptr->addAttribute(XmlAttrPtr(new XmlAttribute(**im)));
  for(std::list<XmlNodePtr>::const_iterator im=children.begin(); im!=children.end(); im++)
    ptr->addChild((*im)->clone());
  return ptr;

}

/***********************************************/

UInt XmlNode::getChildCount(const QStringList &names)
{
  UInt count = 0;
  for(std::list<XmlNodePtr>::iterator im=children.begin(); im!=children.end(); im++)
    if(names.contains((*im)->getName()))
      count++;
  return count;
}

/***********************************************/

XmlNodePtr XmlNode::getChild(const QStringList &names)
{
  XmlNodePtr ptr(nullptr);

  for(std::list<XmlNodePtr>::iterator im=children.begin(); im!=children.end(); im++)
    if(names.contains((*im)->getName()))
    {
      ptr = *im;
      *im = XmlNodePtr(nullptr);
      children.erase(im);
      return ptr;
    }

  return ptr;
}

/***********************************************/

XmlNodePtr XmlNode::findChild(const QString &name)
{
  for(std::list<XmlNodePtr>::iterator im=children.begin(); im!=children.end(); im++)
    if((*im)->getName()==name)
      return *im;

  return XmlNodePtr(nullptr);
}

/***********************************************/

XmlNodePtr XmlNode::getNextChild()
{
  XmlNodePtr ptr(nullptr);
  if(!children.empty())
  {
    ptr = children.front();
    children.front() = XmlNodePtr(nullptr);
    children.pop_front();
  }
  return ptr;
}

/***********************************************/
/***********************************************/

XmlAttrPtr XmlNode::getAttribute(const QString &name)
{
  for(std::list<XmlAttrPtr>::iterator im=attribute.begin(); im!=attribute.end(); im++)
    if((*im)->getName()==name)
    {
      XmlAttrPtr ptr = *im;
      *im = XmlAttrPtr(nullptr);
      attribute.erase(im);
      return ptr;
    }

  return XmlAttrPtr(nullptr);
}

/***********************************************/

XmlAttrPtr XmlNode::getNextAttribute()
{
  XmlAttrPtr ptr(nullptr);
  if(!children.empty())
  {
    ptr = attribute.front();
    attribute.front() = XmlAttrPtr(nullptr);
    attribute.pop_front();
  }
  return ptr;
}

/***********************************************/

XmlNodePtr XmlNode::readFile(const QString &name)
{
  try
  {
    QFile file(name);
    if(!file.open(QFile::ReadOnly | QFile::Text))
      throw(Exception("cannot open file"));

    int          errorLine, errorColumn;
    QString      errorStr;
    QDomDocument doc;
    if (!doc.setContent(&file, false, &errorStr, &errorLine, &errorColumn))
    {
      std::stringstream ss;
      ss<<"parser error: (row="<<errorLine<<", col="<<errorColumn<<errorStr.toStdString();
      throw(Exception(ss.str()));
    }
    return create(doc.documentElement());
  }
  catch(std::exception &e)
  {
    throw(Exception(_GROOPS_ERRORLINE+" filename='"+name.toStdString()+"'\n"+e.what()));
  }
}

/***********************************************/

XmlNodePtr XmlNode::create(const QDomElement &dom)
{
  XmlNodePtr node(new XmlNode(dom.tagName()));

  if(dom.hasAttributes())
  {
    QDomNamedNodeMap att = dom.attributes();
    for(int i=0; i<att.length(); i++)
      node->addAttribute(XmlAttrPtr(new XmlAttribute(att.item(i).toAttr().name(), att.item(i).toAttr().value())));
  }

  QDomElement dom2 = dom;
  QDomElement domChild = dom2.firstChildElement();
  while(!domChild.isNull())
  {
    node->addChild(create(domChild));
    QDomElement domChild2 = domChild;
    domChild = domChild.nextSiblingElement();
    dom2.removeChild(domChild2);
  }

  node->text = dom2.text();

  return node;
}

/***********************************************/

void XmlNode::writeFile(const QString &name, const XmlNodePtr &root)
{
  try
  {
    QFile file(name);
    if(file.open(QFile::WriteOnly | QFile::Truncate))
    {
      QTextStream stream(&file);
      write(stream, root);
    }
    else
      throw(Exception("cannot open file"));
  }
  catch(std::exception &e)
  {
    throw(Exception(_GROOPS_ERRORLINE+" filename='"+name.toStdString()+"'\n"+e.what()));
  }
}

/***********************************************/

void XmlNode::write(QTextStream &stream, const XmlNodePtr &root)
{
  try
  {
    //stream.exceptions(std::ios::badbit|std::ios::failbit);
    stream<<"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    root->write(stream);
  }
  catch(std::exception &e)
  {
    throw(Exception(std::string("In XmlNode::write:\n")+e.what()));
  }
}

/***********************************************/

void XmlNode::write(QTextStream &stream, UInt depth)
{
  // start tag
  for(UInt i=0; i<depth; i++)
    stream<<"\t";
  stream<<"<"<<getName();

  // attribute
  for(std::list<XmlAttrPtr>::iterator im=attribute.begin(); im!=attribute.end(); im++)
    stream<<" "<<(*im)->getName()<<"=\""<<sanitizeXML((*im)->getText())<<"\"";

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
    for(std::list<XmlNodePtr>::iterator im=children.begin(); im!=children.end(); im++)
      (*im)->write(stream, depth+1);
    for(UInt i=0; i<depth; i++)
      stream<<"\t";
  }

  // end tag
  stream<<"</"<<getName()<<">\n";
}

/***********************************************/

QString XmlNode::sanitizeXML(const QString& in)
{
  if(in.isEmpty())
    return QString();

  QString out(in);

  std::vector<QString> chars = {"&", ">", "<", "\"", "'"};
  std::vector<QString> replacements = {"&amp;", "&gt;", "&lt;", "&quot;", "&apos;"};

  for(UInt k = 0; k<chars.size(); k++)
      out = out.replace(chars.at(k), replacements.at(k));

  return out;
}

/***********************************************/
/***** FUNCTIONS *******************************/
/***********************************************/

XmlNodePtr getNextChild(const XmlNodePtr &ptr, QString &name)
{
  XmlNodePtr child = ptr->getNextChild();
  name = child->getName();
  return child;
}

/***********************************************/

UInt childCount(const XmlNodePtr &parent, const QString &name, Bool mustGreaterZero)
{
  UInt count = parent->getChildCount(name);
  if((count==0) && mustGreaterZero)
    throw(Exception(("'"+parent->getName()+"' must contain '"+name+"'").toStdString()));
  return count;
}


/***********************************************/
