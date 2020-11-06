/***********************************************/
/**
* @file treeElementAdd.h
*
* @brief The add element for unbounded elements.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTADD__
#define __GROOPSGUI__TREEELEMENTADD__

#include "base/importGroops.h"
#include "tree/treeElement.h"

/***** CLASS ***********************************/

class TreeElementAdd : public TreeElement
{
  Q_OBJECT

public:
  UInt unboundedCount;
  Bool visible;

public:
  TreeElementAdd(Tree *tree, TreeElementComplex *parentElement,
                 XsdElementPtr xsdElement, Bool visible=true);
  virtual ~TreeElementAdd() override;

/** @brief Generate XML-tree. *
* the add button will not be saved. */
virtual XmlNodePtr getXML(Bool /*withEmptyNodes*/=false) const override {return XmlNodePtr(nullptr);}

/** @brief is this the add-element?.
* @return true. */
virtual Bool isElementAdd() const override {return true;}

/** @brief returns the add element.
* @return this */
virtual TreeElementAdd *elementAdd() const override {return const_cast<TreeElementAdd*>(this);}

/** @brief creates visual element (if visible). */
TreeItem *createItem(TreeItem *parent, TreeItem *after) override;

/** @brief creates a push button. */
QWidget *createEditor() override;

/** @brief Adds element. */
virtual void interact() override { pushButtonClicked(); }

public slots:
  void pushButtonClicked();
};

/***********************************************/

#endif
