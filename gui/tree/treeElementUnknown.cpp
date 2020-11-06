/***********************************************/
/**
* @file treeElementUnknown.cpp
*
* @brief Unknown element without XML schema entry.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementUnknown.h"

/***********************************************/

TreeElementUnknown::TreeElementUnknown(Tree *tree, TreeElementComplex *parentElement, XmlNodePtr xmlNode)
  : TreeElementComplex(tree, parentElement, XsdElementPtr(nullptr), "", xmlNode)
{
  try
  {
    if(!xmlNode)
      throw(Exception("no XML given"));

    _isEditable = !(xmlNode->hasChildren());
    addChoice(xmlNode->getText(), XsdElementPtr(nullptr), QJsonObject());
    createChildrenElements(0, xmlNode);
    setSelectedIndex(isLinked() ? selectedIndex() : 0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlNodePtr TreeElementUnknown::getXML(Bool withEmptyNodes) const
{
  try
  {
    XmlNodePtr xmlNode = TreeElement::getBaseXML();
    if(xmlNode==nullptr)
      return xmlNode;
    if(!isLinked())
      xmlNode->setText(selectedValue());
    getChildrenXML(xmlNode, withEmptyNodes);
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
