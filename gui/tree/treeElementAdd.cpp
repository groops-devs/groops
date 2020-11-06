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
#include "tree/treeElementAdd.h"
#include "addGlobalDialog/addGlobalDialog.h"

/***********************************************/

TreeElementAdd::TreeElementAdd(Tree *tree, TreeElementComplex *parentElement,
                               XsdElementPtr xsdElement, Bool visible)
  : TreeElement(tree, parentElement, xsdElement, "", XmlNodePtr(nullptr))
{
  try
  {
    this->unboundedCount = 0;
    this->visible = visible;

    insertNewValue("...", false);
    setSelectedIndex(0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElementAdd::~TreeElementAdd()
{
  try
  {
    if(unboundedCount!=0)
      throw(Exception("unboundedCount is not empty"));
  }
  catch(std::exception &e)
  {
    qDebug() << QString::fromStdString("Exception in destructor at "+_GROOPS_ERRORLINE+"\n"+e.what());
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
    if(parentElement && parentElement == tree->elementGlobal())
    {
      AddGlobalDialog dialog(tree->elementGlobal(), tree);
      if(dialog.exec())
      {
        QString label = dialog.elementName();
        tree->elementGlobal()->addChild(nullptr, dialog.elementType()->type, XmlNodePtr(nullptr), label, true);
      }
    }
    else if(parentElement)
      parentElement->addChild(this, type(), XmlNodePtr(nullptr));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/
