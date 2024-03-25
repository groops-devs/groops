/***********************************************/
/**
* @file treeElementAdd.cpp
*
* @brief The add element for unbounded elements.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#include <QPushButton>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementGlobal.h"
#include "tree/treeElementProgram.h"
#include "tree/treeElementAdd.h"
#include "addVariableDialog/addVariableDialog.h"

/***********************************************/

TreeElementAdd::TreeElementAdd(Tree *tree, TreeElementComplex *parentElement,
                               XsdElementPtr xsdElement, bool visible)
  : TreeElement(tree, parentElement, xsdElement, "", XmlNodePtr(nullptr)), unboundedCount(0), visible(visible)
{
  try
  {
    insertNewValue("...", false);
    setSelectedIndex(0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

TreeItem *TreeElementAdd::createItem(TreeItem *parent, TreeItem *after)
{
  try
  {
    if(!visible)
      return nullptr;
    return TreeElement::createItem(parent, after);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QWidget *TreeElementAdd::createEditor()
{
  try
  {
    // create Pushbutton
    QPushButton *pushButton = new QPushButton(QIcon(":/icons/scalable/edit-add.svg"), tr("&Add"), tree);
    pushButton->setDefault(true);

    // signals and slots connections
    connect(pushButton, SIGNAL(clicked()), this, SLOT(pushButtonClicked()));

    return pushButton;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementAdd::pushButtonClicked()
{
  try
  {
    if(dynamic_cast<TreeElementGlobal*>(parentElement))
      parentElement->addChild(this, "", "", XmlNodePtr(nullptr));
    else
      parentElement->addChild(this, type(), "", XmlNodePtr(nullptr));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/
