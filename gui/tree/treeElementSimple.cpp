/***********************************************/
/**
* @file treeElementSimple.cpp
*
* @brief Element without children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#include <QtDebug>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementSimple.h"

/***********************************************/

TreeElementSimple::TreeElementSimple(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                                     const QString &defaultOverride, XmlNodePtr xmlNode, Bool fromFile)
  : TreeElement(tree, parentElement, xsdElement, defaultOverride, xmlNode)
{
  try
  {
    if(!isLinked())
    {
      if(xmlNode && xmlNode->hasChildren())
        throw(Exception("xml node doesn't match with schema"));
      else if(xmlNode && !xmlNode->getText().isEmpty())
        insertNewValue(xmlNode->getText(), false);
      else if(fromFile || defaultValue().isEmpty())
        insertNewValue("", false);
    }
    if(!defaultValue().isEmpty())
      insertNewValue(defaultValue(), false);

    setSelectedIndex((isLinked() && selectedIndex() > 0) ? selectedIndex() : 0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlNodePtr TreeElementSimple::getXML(Bool withEmptyNodes) const
{
  try
  {
    if(selectedValue().isEmpty() && !withEmptyNodes)
       return XmlNodePtr(nullptr);
    XmlNodePtr xmlNode = TreeElement::getBaseXML();
    if((xmlNode==nullptr) || isLinked())
      return xmlNode;
    xmlNode->setText(selectedValue());
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
