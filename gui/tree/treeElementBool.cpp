/***********************************************/
/**
* @file treeElementBool.cpp
*
* @brief Element with true or false.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementBool.h"

/***********************************************/

TreeElementBool::TreeElementBool(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                                 const QString &defaultOverride, XmlNodePtr xmlNode)
  : TreeElement(tree, parentElement, xsdElement, defaultOverride, xmlNode)
{
  try
  {
    insertNewValue("no",  false);
    insertNewValue("yes", false);

    int index=0;
    if(xmlNode)
    {
      if(xmlNode->getText()=="1")
        index=1;
      else if(!isLinked() && xmlNode->getText()!="0")
        throw(Exception("xml node doesn't match with schema"));
    }
    else if(defaultValue()=="1")
      index=1;

    setSelectedIndex(isLinked() ? selectedIndex() : index);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

XmlNodePtr TreeElementBool::getXML(Bool /*withEmptyNodes*/) const
{
  try
  {
    XmlNodePtr xmlNode = TreeElement::getBaseXML();
    if((xmlNode==nullptr) || isLinked())
      return xmlNode;
    if(selectedValue()=="yes")
      xmlNode->setText("1");
    else
      xmlNode->setText("0");
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
