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
  bool visible;

public:
  TreeElementAdd(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                 const QString &defaultOverride, bool visible=true);

  bool optional()          const override {return false;}
  bool unbounded()         const override {return true;}
  bool isRenamedInSchema() const override {return false;}

  /** @brief Generate XML-tree. *
  * the add button will not be saved. */
  XmlNodePtr createXmlTree(bool /*createRootEvenIfEmpty*/) const override {return XmlNodePtr(nullptr);}

  /** @brief returns the add element.
  * @return this */
  TreeElementAdd *elementAdd() const override {return const_cast<TreeElementAdd*>(this);}

  bool canSetLoop()      const override {return false;}
  bool canSetCondition() const override {return false;}
  bool canDisabled()     const override {return false;}
  bool canComment()      const override {return false;}

  /** @brief creates visual element (if visible). */
  TreeItem *createItem(TreeItem *parent, TreeItem *after) override;

  /** @brief creates a push button. */
  QWidget *createEditor() override;

  /** @brief Adds element. */
  void interact() override {pushButtonClicked();}

public slots:
  void pushButtonClicked();
};

/***********************************************/

#endif
