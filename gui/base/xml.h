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
#include <QStringList>
#include <QFile>
#include <sstream>
#include "base/importGroops.h"

/***** TYPES ***********************************/

class   QXmlStreamReader;
class   QTextStream;
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
  /** @brief Name of the attribute. */
  QString name;

  /** @brief Content/text of the attribute. */
  QString text;

  XmlAttribute() {}
  XmlAttribute(QString _name, QString _text) : name(_name), text(_text) {}
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

  static XmlNodePtr parse(QXmlStreamReader &reader, QString &errorMessage);
  void write(QTextStream &stream, UInt depth, bool writeCommentsAsElements);

public:
  XmlNode(const QString &name) {this->name = name;}
  XmlNode(const XmlNode &node) = delete;
  XmlNode &operator=(const XmlNode &node) = delete;

  /** @brief Returns the name of the node. */
  const QString &getName() const {return name;}

  /** @brief Returns the text/content of the node. */
  const QString &getText() const {return text;}

  /** @brief Sets the text/content of the node. */
  void setText(const QString &text) {this->text = text;}

  /** @brief Returns TRUE of node has children, FALSE otherwise. */
  bool hasChildren() const {return !children.empty();}

  /** @brief Returns number of children with name @p name.  */
  UInt getChildCount(const QStringList &names);
  UInt getChildCount(const QString &name) {return getChildCount(QStringList(name));}

  /** @brief Returns the first child with name @p name.
  * The child is removed from the node tree.
  * Returns a nullptr if child does not exist.  */
  XmlNodePtr getChild(const QStringList &names);
  XmlNodePtr getChild(const QString &name) {return getChild(QStringList({name}));}

  /** @brief Adds a child to the node.
  * A node must not be added to the node tree multiple times. */
  void addChild(const XmlNodePtr &child) {children.push_back(child);}

  /** @brief Returns the next child.
  * The child is removed from the node tree.
  * Returns a nullptr if child does not exist. */
  XmlNodePtr getNextChild();

  /** @brief Returns the next child.
  * Returns a nullptr if child does not exist. */
  XmlNodePtr peekNextChild();

  /** @brief Returns the attribute with name @p name.
  * Returns a nullptr if attribute does not exist. */
  XmlAttrPtr getAttribute(const QString &name);

  /** @brief Adds a new attribute to the node. */
  void addAttribute(const XmlAttrPtr &attr) {if(attr) attribute.push_back(attr);}

  /** @brief Creates a new node with name @p name. */
  static XmlNodePtr create(const QString &name) {return std::make_shared<XmlNode>(name);}

  /** @brief Reads an XML tree from stream.
  * @param name File name.
  * @return Root of the XML tree. */
  static XmlNodePtr read(const QByteArray &data, QString &errorMessage);

  /** @brief Reads an XML tree from stream.
  * @param name File name.
  * @return Root of the XML tree. */
  static XmlNodePtr read(QIODevice *device, QString &errorMessage);

  /** Writes an XML tree into a @p stream.
  * @param stream Output stream.
  * @param root Root of the XML tree. */
  static void write(QTextStream &stream, const XmlNodePtr &root, bool writeCommentsAsElements=false);

  /** Writes an XML tree into a file.
  * @param name File name.
  * @param root Root of the XML tree. */
  static void writeFile(const QString &name, const XmlNodePtr &root);
};

/***** FUNCTIONS ***********************************/

XmlAttrPtr readAttribute (const XmlNodePtr &node, const QString &name, QString &var);
XmlAttrPtr writeAttribute(const XmlNodePtr &node, const QString &name, const QString &var);

/***********************************************/

#endif /* __GROOPSGUI__XML__ */
