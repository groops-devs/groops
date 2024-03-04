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
                                 const QString &defaultOverride, XmlNodePtr xmlNode, bool /*fillWithDefaults*/)
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
    else if(defaultValue=="1")
      index=1;

    setSelectedIndex(isLinked() ? selectedIndex() : index);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlNodePtr TreeElementBool::createXmlTree(bool /*createRootEvenIfEmpty*/) const
{
  try
  {
    XmlNodePtr xmlNode = TreeElement::createXmlBaseNode();
    if(!isLinked())
      xmlNode->setText(selectedValue()=="yes" ? "1" : "0");
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementBool::overwrite(const QString &type, XmlNodePtr xmlNode, bool contentOnly)
{
  try
  {
    if(!canOverwrite(type) || !xmlNode)
      return false;

    tree->undoStack->beginMacro("overwrite "+name());
    if(baseOverwrite(xmlNode, contentOnly))
      changeSelectedIndex(xmlNode->getText() == "1");
    tree->undoStack->endMacro();
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
