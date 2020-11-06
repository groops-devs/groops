/***********************************************/
/**
* @file xml.h
*
* @brief Nodes of a XML tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2004-11-14
*
*/
/***********************************************/

#ifndef __GROOPS_XML__
#define __GROOPS_XML__

#include "base/importStd.h"
#include "base/angle.h"
#include "base/time.h"

/** @addtogroup parserGroup */
/// @{

/***** TYPES ***********************************/

class   XmlAttribute;
class   XmlNode;
typedef std::shared_ptr<XmlNode>      XmlNodePtr;
typedef std::shared_ptr<XmlAttribute> XmlAttrPtr;

/***** CLASS ***********************************/

/** @brief Attribute of a XML-node.
* @ingroup parserGroup
* @see XmlNode */
class XmlAttribute
{
public:
  /// Name of Attribute.
  std::string name;

  /// Content of Attribute.
  std::string text;

  /// Name of Attribute.
  const std::string &getName() const {return name;}

  /// Content of Attribute.
  const std::string &getText() const {return text;}

  /** @brief Interpret content as type of @a var.
  * @param[out] var is filled with the content of attribute. */
  template<typename T> void getValue(T &var) const;

  /** @brief Set content of attribute. */
  template<typename T> void setValue(const T &var);
};

/***** CLASS ***********************************/

/** @brief Node of a XML-tree.
* @ingroup parserGroup
* For simple XML files. This means only on block of text is allowed for each node.
*
* Der Baum wird beim auslesen direkt abgebaut,
* so dass nur einmal auslesen moeglich ist.
* Der Speicher wird mit std::shared_ptr verwaltet,
* so dass man keinen Speicher freigeben muss. */
class XmlNode
{
  std::string           name_;
  std::string           text_;
  std::list<XmlNodePtr> children;
  std::list<XmlAttrPtr> attribute;

  void write(std::ostream &stream, UInt depth=0);

public:
  /// Constructor.
  explicit XmlNode(const std::string &name) : name_(name) {}

  XmlNode(const XmlNode &node) = delete;
  XmlNode &operator=(const XmlNode &node) = delete;

  /** @brief Deep copy.
  * Creates a copy of the node and of all children. */
  XmlNodePtr clone() const;

  /** @brief Name of the node. */
  const std::string &getName() const {return name_;}

  /** @brief Set name of the node. */
  void setName(const std::string &name) {name_ = name;}

  /** @brief Content of the node. */
  const std::string &getText() const {return text_;}

  /** @brief Set content of the node. */
  void setText(const std::string &text) {text_ = text;}

  /** @brief Append @a text to the content of the node. */
  void addText(const std::string &text) {text_ += text;}

  /** @brief Interpret content as type of @a var.
  * @param[out] var is filled with the content of node. */
  template<typename T> void getValue(T &var) const;

  /** @brief Set content of node. */
  template<typename T> void setValue(const T &var);

  /** @brief Has the node children nodes? */
  Bool hasChildren() const {return !children.empty();}

  std::list<XmlNodePtr> &getChildren() {return children;}

  /** @brief Number of children with @a name. */
  UInt getChildCount(const std::string &name);

  /** @brief Returns the the first child with @a name.
  * The child is removed from tree. If child does not exist, a NULL pointer is returned. */
  XmlNodePtr getChild(const std::string &name);

  /** @brief Returns the the first child with @a name.
  * The child is NOT removed from tree. If child does not exist, a NULL pointer is returned. */
  XmlNodePtr findChild(const std::string &name);

  /** @brief Append a new child.
  * It is not allowed to have the same node multiple times in the tree.
  * (Create a copy with @a clone() before). */
  void addChild(const XmlNodePtr &child) {children.push_back(child);}

  /** @brief Insert a new child at begin.
  * It is not allowed to have the same node multiple times in the tree.
  * (Create a copy with @a clone() before). */
  void prependChild(const XmlNodePtr &child) {children.push_front(child);}

  /** @brief Returns the next child.
  * The child is removed from tree. If child not exists, a NULL pointer is returned. */
  XmlNodePtr getNextChild();

  /** @brief Has the node attributes? */
  Bool hasAttribute() const {return !attribute.empty();}

  /** @brief Returns the the first attribute with @a name.
  * The attribute is removed from the node. If attribute does not exist, a NULL pointer is returned. */
  XmlAttrPtr getAttribute(const std::string &name);

  /** @brief Returns the the first attribute with @a name.
  * The attribute is NOT removed from the node. If attribute does not exist, a NULL pointer is returned. */
  XmlAttrPtr findAttribute(const std::string &name);

  /** @brief Append a new attribute.
  * It is not allowed to have the same attribute multiple times in the tree. */
  void addAttribute(const XmlAttrPtr &attr) {attribute.push_back(attr);}

  /** @brief Returns the next attribute.
  * The attribute is removed from the node. If attribute does not exist, a NULL pointer is returned. */
  XmlAttrPtr getNextAttribute();

  /** @brief XmlNode factory. */
  static XmlNodePtr create(const std::string &name);

  /** @brief Read XML tree from stream.
  * @param stream Input stream.
  * @return Root node of the XML tree. */
  static XmlNodePtr read(std::istream &stream);

  /** @brief Write XML tree to stream.
  * @param stream output stream.
  * @param root Root node of the XML tree. */
  static void write(std::ostream &stream, const XmlNodePtr &root);

  /** @brief Replace characters with entity references
  * @param in string to be parsed
  * @return out string with character entitiy references
  */
  static std::string sanitizeXML(const std::string& in);
};


/***** FUNCTIONS ***********************************/

/** @brief Anzahl der Kinder mit Namen @a name.
* If @a mustSet==TRUE, an Exception is thrown.
* @relates XmlNode */
UInt childCount(const XmlNodePtr &parent, const std::string &name, Bool mustGreaterZero=FALSE);

/** @brief Returns the next child.
* In @a name the name of the child is returned.
* If child not exists, a NULL pointer is returned.
* @relates XmlNode */
XmlNodePtr getNextChild(const XmlNodePtr &parent, std::string &name);

/** @brief Returns the the first child with @a name.
* The child is removed from tree. If child not exists and @a mustSet==TRUE,
* an Exception is thrown otherwise a NULL pointer is returned.
* @relates XmlNode */
XmlNodePtr getChild(const XmlNodePtr &parent, const std::string &name, Bool mustSet=FALSE);

/** @brief Reads a variable from XML.
* Locking for child with @a name and interpret the content as type of @a var.
* The child is removed from tree. If child not exists the @a var is untouched
* and if additionally @a mustSet==TRUE, an Exception is thrown.
* @relates XmlNode */
template<class T>
XmlNodePtr readXml(const XmlNodePtr &parent, const std::string &name, T &var, Bool mustSet=FALSE);

/** @brief Write a variable to a XML tree.
* A XML child with @a name is created and is appended to @a parent.
* The value of @a var is written as content of the child.
* @relates XmlNode */
template<class T>
XmlNodePtr writeXml(const XmlNodePtr &parent, const std::string &name, const T &var);

/** @brief Reads a variable from XmlAttribute.
* Locking for a attribute with @a name and interpret the content as type of @a var.
* The attribute is removed from node. If attribute not exists the @a var is untouched
* and if additionally @a mustSet==TRUE, an Exception is thrown.
* @relates XmlNode */
template<class T>
XmlAttrPtr readAttribute(const XmlNodePtr &node, const std::string &name, T &var, Bool mustSet=FALSE);

/** @brief Write a variable as attribute to a XML node.
* A XmlAttribute with @a name is created and is appended to @a node.
* The value of @a var is written as content of the attribute.
* @relates XmlNode */
template<class T>
XmlAttrPtr writeAttribute(const XmlNodePtr &node, const std::string &name, const T &var);

/** @brief XmlNode factory.
* A XML child with @a name is created and is appended to @a parent.
* @relates XmlNode */
XmlNodePtr createXmlNode(const XmlNodePtr &parent, const std::string &name);

/// @}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

template<typename T> void XmlAttribute::getValue(T &var) const
{
  try
  {
    std::stringstream ss(text);
    ss>>var;
    if(!(ss.good()||ss.eof()))
      throw Exception("stream error");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

template<> inline void XmlAttribute::getValue(std::string &var) const // Spezialization
{
  var = text;
}

/***********************************************/

inline XmlNodePtr XmlNode::create(const std::string &name)
{
  return XmlNodePtr(new XmlNode(name));
}

/***********************************************/

template<typename T> void XmlAttribute::setValue(const T &var)
{
  try
  {
    std::stringstream ss;
    ss.setf(std::ios::scientific,std::ios::floatfield);
    ss.precision(14);
    ss<<var; ss>>text;
    if(!(ss.good()||ss.eof()))
      throw Exception("stream error");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

template<> inline void XmlAttribute::setValue(const std::string &var) // Spezialization
{
  text = var;
}

/***********************************************/

template<typename T> void XmlNode::getValue(T &var) const
{
  try
  {
    std::stringstream ss(text_);
    ss>>var;
    if(!(ss.good()||ss.eof()))
      throw Exception("stream error");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

template<> inline void XmlNode::getValue(Double &var) const // Spezialization
{
  try
  {
    var = std::stod(text_);
  }
  catch(std::invalid_argument &)
  {
    throw(Exception("cannot read number: "+text_));
  }
}

template<> inline void XmlNode::getValue(Time &var) const // Spezialization
{
  Double mjd;
  getValue(mjd);
  var = mjd2time(mjd);
}

template<> inline void XmlNode::getValue(Angle &var) const // Spezialization
{
  Double x;
  getValue(x);
  var = Angle(x*DEG2RAD);
}

template<> inline void XmlNode::getValue(std::string &var) const
{
  var = getText();
}

/***********************************************/

template<typename T> void XmlNode::setValue(const T &var)
{
  try
  {
    std::stringstream ss;
    ss.setf(std::ios::scientific,std::ios::floatfield);
    ss.precision(14);
    ss<<var; ss>>text_;
    if(!(ss.good()||ss.eof()))
      throw Exception("stream error");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

template<> inline void XmlNode::setValue(const Time  &var) // Spezialization
{
  setValue(var.mjd());
}

template<> inline void XmlNode::setValue(const Angle &var) // Spezialization
{
  setValue(static_cast<Double>(var*RAD2DEG));
}

template<> inline void XmlNode::setValue(const std::string &var) // Spezialization
{
  text_ = var;
}

/***********************************************/

template<class T> inline XmlNodePtr readXml(const XmlNodePtr &parent, const std::string &name, T &var, Bool mustSet)
{
  XmlNodePtr child = getChild(parent, name, mustSet);
  if(child!=nullptr)
    child->getValue(var);
  return child;
}

/***********************************************/

template<class T> inline XmlNodePtr writeXml(const XmlNodePtr &parent, const std::string &name, const T &var)
{
  XmlNodePtr child = createXmlNode(parent, name);
  child->setValue(var);
  return child;
}

/***********************************************/

inline XmlNodePtr getChild(const XmlNodePtr &parent, const std::string &name, Bool mustSet)
{
  try
  {
    XmlNodePtr child = parent->getChild(name);
    if((child==nullptr) && mustSet)
      throw(Exception("'"+parent->getName()+"' must contain '"+name+"'"));
    return child;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline XmlNodePtr createXmlNode(const XmlNodePtr &parent, const std::string &name)
{
  XmlNodePtr child = XmlNode::create(name);
  parent->addChild(child);
  return child;
}

/***********************************************/

template<class T>
inline XmlAttrPtr readAttribute(const XmlNodePtr &node, const std::string &name, T &var, Bool mustSet)
{
  try
  {
    XmlAttrPtr attr = node->getAttribute(name);
    if((attr==nullptr) && mustSet)
      throw(Exception("'"+node->getName()+"' muss contain '"+name+"'"));
    if(attr!=nullptr)
      attr->getValue(var);
    return attr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<class T>
inline XmlAttrPtr writeAttribute(const XmlNodePtr &node, const std::string &name, const T &var)
{
  XmlAttrPtr attr = XmlAttrPtr(new XmlAttribute());
  attr->name = name;
  attr->setValue(var);
  node->addAttribute(attr);
  return attr;
}

/***********************************************/
/***********************************************/

#endif /* __GROOPS_XML__ */
