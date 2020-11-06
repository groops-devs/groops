/***********************************************/
/**
* @file xml.h
*
* @brief Nodes of an XML tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2004-11-14
*/
/***********************************************/

#ifndef __GROOPSGUI__XML__
#define __GROOPSGUI__XML__

#include <QString>
#include <QTextStream>
#include <sstream>
#include "base/importGroops.h"

/***** TYPES ***********************************/

class   QDomElement;
class   XmlAttribute;
class   XmlNode;
typedef std::shared_ptr<XmlNode>      XmlNodePtr;
typedef std::shared_ptr<XmlAttribute> XmlAttrPtr;

/***** CLASS ***********************************/

/**
* @brief Attribues of an XML node.
* @see XmlNode
*/
class XmlAttribute
{
public:
/**
* @brief Name of the attribute.
*/
QString name;

/**
* @brief Content/text of the attribute.
*/
QString text;

XmlAttribute() {}
XmlAttribute(QString _name, QString _text) : name(_name), text(_text) {}

/**
* @brief Returns name of the attribute.
*/
const QString &getName() const {return name;}

/**
* @brief Returns content/text of the attribute.
*/
const QString &getText() const {return text;}

/**
* @brief Return text/content of attribute interpreted as a value.
*/
template<typename T>
void getValue(T &var) const;

/**
* @brief Set value as text/content of attribute.
*/
template<typename T>
void setValue(const T &var);
};

/***** CLASS ***********************************/

/**
* @brief Node of an XML tree.
*
* A node can only have one(!) text part that comes before any children.
* Nodes are removed from the node tree while reading, so each node can only be read once.
* Memory is managed via std::shared_prt, so it does not have to be freed manually.
*/
class XmlNode
{
  QString           name;
  QString           text;
  std::list<XmlNodePtr> children;
  std::list<XmlAttrPtr> attribute;

  void write(QTextStream &stream, UInt depth=0);
  XmlNode(const XmlNode &node);
  XmlNode &operator = (const XmlNode &node);

public:
  XmlNode(const QString &name) {this->name = name;}

/**
* @brief Create a (deep) copy of this node and all of its children.
*/
XmlNodePtr clone() const;

/**
* @brief Returns the name of the node.
*/
const QString &getName() const {return name;}

/**
* @brief Set the name of the node.
*/
void setName(const QString &name) {this->name = name;}

//@{
/**
* @brief Returns the text/content of the node.
*/
const QString &getText() const {return text;}

/**
* @brief Sets the text/content of the node.
*/
void setText(const QString &text) {this->text = text;}

/**
* @brief Append @p text to the text/content of a node.
*/
void addText(const QString &text) {this->text += text;}

/**
* @brief Return text/content of attribute interpreted as a value.
*/
template<typename T>
void getValue(T &var) const;

/**
* @brief Set value as text/content of attribute.
*/
template<typename T>
void setValue(const T &var);
//@}

//@{
/**
* @brief Returns TRUE of node has children, FALSE otherwise.
*/
Bool hasChildren() const {return !children.empty();}

/**
* @brief Returns number of children with name @p name.
*/
UInt getChildCount(const QString &name) { return getChildCount(QStringList({name})); }
UInt getChildCount(const QStringList &names);

/**
* @brief Returns the first child with name @p name.
*
* The child is removed from the node tree.
* Returns a nullptr if child does not exist.
*/
XmlNodePtr getChild(const QString &name) { return getChild(QStringList({name})); }
XmlNodePtr getChild(const QStringList &names);

/**
* @brief Returns the first child with name @p name.
*
* The child is not (!) removed from the node tree.
* Returns a nullptr if child does not exist.
*/
XmlNodePtr findChild(const QString &name);

/**
* @brief Adds a child to the node.
*
* A node must not be added to the node tree multiple times.
* (Clone with @a clone().)
*/
void addChild(const XmlNodePtr &child) {children.push_back(child);}

/**
* @brief Returns the next child.
*
* The child is removed from the node tree.
* Returns a nullptr if child does not exist.
*/
XmlNodePtr getNextChild();
//@}

//@{
/**
* @brief Returns true of the node has an attribute, FALSE otherwise.
*/
Bool hasAttribute() const {return !attribute.empty();}

/**
* @brief Returns the attribute with name @p name.
*
* Returns a nullptr if attribute does not exist.
*/
XmlAttrPtr getAttribute(const QString &name);

/**
* @brief Adds a new attribute to the node.
*/
void addAttribute(const XmlAttrPtr &attr) { if(attr) attribute.push_back(attr);}

/**
* @brief Returns the next attribute.
*
* The attribute is removed from the node.
* Returns a nullptr if attribute does not exist.
*/
XmlAttrPtr getNextAttribute();
//@}

/**
* @brief Creates a new node with name @p name.
*/
static XmlNodePtr create(const QString &name);

static XmlNodePtr create(const QDomElement &dom);

/**
* @brief Reads an XML tree from a @p stream.
* @param stream Input stream.
* @return Root of the XML tree.
*/
static XmlNodePtr read(std::istream &stream);

/**
* @brief Reads an XML tree from a file.
* @param name File name.
* @return Root of the XML tree.
*/
static XmlNodePtr readFile(const QString &name);


/**
* Writes an XML tree into a @p stream.
* @param stream Output stream.
* @param root Root of the XML tree.
*/
static void write(QTextStream &stream, const XmlNodePtr &root);

/**
* Writes an XML tree into a file.
* @param name File name.
* @param root Root of the XML tree.
*/
static void writeFile(const QString &name, const XmlNodePtr &root);

/**
* Replace characters with entity references.
* @param in QString.
* @return out QString with character entitiy references.
*/
static QString sanitizeXML(const QString& in);

};

/***** FUNCTIONS ***********************************/

XmlNodePtr createXmlNode(const XmlNodePtr &parent, const QString &name);

/**
* @brief Returns number of children of node @p parent with name @p name.
*
* If @p mustGreaterZero==true an exception is thrown if there are no children.
*/
UInt childCount(const XmlNodePtr &parent, const QString &name, Bool mustGreaterZero=false);

/**
* @brief Returns the next child of node @p parent.
*
* @param[out] name Name of the child.
*/
XmlNodePtr getNextChild(const XmlNodePtr &parent, QString &name);

/**
* @brief Returns the next child of node @p parent with name @p name.
*
* If @p mustSet==true an exception is thrown if there are no children with @p name.
*/
XmlNodePtr getChild(const XmlNodePtr &parent, const QString &name, Bool mustSet=false);

template<class T>
XmlNodePtr readXml(const XmlNodePtr &parent, const QString &name, T &var, Bool mustSet=false);

template<class T>
XmlNodePtr writeXml(const XmlNodePtr &parent, const QString &name, const T &var);

//@{
template<class T>
XmlAttrPtr readAttribute(const XmlNodePtr &node, const QString &name, T &var, Bool mustSet=false);

template<class T>
XmlAttrPtr writeAttribute(const XmlNodePtr &node, const QString &name, const T &var);
//@}


/***********************************************/
/***** INLINES *********************************/
/***********************************************/

inline XmlNodePtr XmlNode::create(const QString &name)
{
  return XmlNodePtr(new XmlNode(name));
}

/***********************************************/

template<typename T>
inline void XmlAttribute::getValue(T &var) const
{
  try
  {
    QTextStream ss(const_cast<QString*>(&text));
    ss>>var;
/*    if(!(ss.good()||ss.eof()))
      throw Exception("stream error");*/
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

template<>
inline void XmlAttribute::getValue(QString &var) const
{
  var = text;
}

/***********************************************/

template<typename T>
inline void XmlAttribute::setValue(const T &var)
{
  std::stringstream ss;
  ss.setf(std::ios::scientific,std::ios::floatfield);
  ss.precision(14);
  ss<<var;
  text = QString::fromStdString(ss.str());
  if(!(ss.good()||ss.eof()))
    throw Exception("XmlAttribute.setValue:\n stream error");
}

/***********************************************/

template<>
inline void XmlAttribute::setValue(const QString &var)
{
  text = var;
}

/***********************************************/

template<typename T>
inline void XmlNode::getValue(T &var) const
{
  std::stringstream ss(text.toStdString());
  ss>>var;
  if(!(ss.good()||ss.eof()))
    throw Exception("XmlNode.getValue:\n stream error: '"+text.toStdString()+"'");
}

template<>
inline void XmlNode::getValue(QString &var) const
{
  var = getText();
}

/***********************************************/

template<typename T>
inline void XmlNode::setValue(const T &var)
{
  std::stringstream ss;
  ss.setf(std::ios::scientific,std::ios::floatfield);
  ss.precision(14);
  ss<<var;
  text = QString::fromStdString(ss.str());
  if(!(ss.good()||ss.eof()))
    throw Exception("XmlNode.setValue:\n stream error");
}

template<>
inline void XmlNode::setValue(const QString &var)
{
  text = var;
}

/***********************************************/

template<class T>
inline XmlNodePtr readXml(const XmlNodePtr &parent, const QString &name, T &var, Bool mustSet)
{
  XmlNodePtr child = getChild(parent, name, mustSet);
  if(child!=nullptr)
    child->getValue(var);
  return child;
}

/***********************************************/

template<class T>
inline XmlNodePtr writeXml(const XmlNodePtr &parent, const QString &name, const T &var)
{
  XmlNodePtr child = createXmlNode(parent, name);
  child->setValue(var);
  return child;
}

/***********************************************/

inline XmlNodePtr getChild(const XmlNodePtr &parent, const QString &name, Bool mustSet)
{
  try
  {
    XmlNodePtr child = parent->getChild(name);
    if((child==nullptr) && mustSet)
      throw(Exception("'"+parent->getName().toStdString()+"' must contain '"+name.toStdString()+"'"));
    return child;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

inline XmlNodePtr createXmlNode(const XmlNodePtr &parent, const QString &name)
{
  XmlNodePtr child = XmlNode::create(name);
  parent->addChild(child);
  return child;
}

/***********************************************/

template<class T>
inline XmlAttrPtr readAttribute(const XmlNodePtr &node, const QString &name, T &var, Bool mustSet)
{
  try
  {
    XmlAttrPtr attr = node->getAttribute(name);
    if((attr==nullptr) && mustSet)
      throw(Exception(("'"+node->getName()+"' muss contain '"+name+"'").toStdString()));
    if(attr!=nullptr)
      attr->getValue(var);
    return attr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

template<class T>
inline XmlAttrPtr writeAttribute(const XmlNodePtr &node, const QString &name, const T &var)
{
  XmlAttrPtr attr = XmlAttrPtr(new XmlAttribute());
  attr->name = name;
  attr->setValue(var);
  node->addAttribute(attr);
  return attr;
}

/***********************************************/
/***********************************************/

#endif /* __GROOPSGUI__XML__ */
