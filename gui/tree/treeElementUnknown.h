/***********************************************/
/**
* @file treeElementUnknown.h
*
* @brief Enknown element without XML schema entry.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTUNKNOWN__
#define __GROOPSGUI__TREEELEMENTUNKNOWN__

#include "base/importGroops.h"
#include "tree/treeElementComplex.h"

/***** CLASS ***********************************/

class TreeElementUnknown : public TreeElementComplex
{
  Q_OBJECT

  bool _isEditable;

public:
  TreeElementUnknown(Tree *tree, TreeElementComplex *parentElement, XmlNodePtr xmlNode);

  bool optional()          const override {return true;}
  bool unbounded()         const override {return true;}
  bool isRenamedInSchema() const override {return false;}
  bool canSetLoop()        const override {return false;}
  bool canSetCondition()   const override {return false;}

  /** @brief Values can be edited. */
  bool isEditable() const override {return _isEditable;}

  /** @brief Generate XML-tree.
  * recursivly called for all children. */
  XmlNodePtr createXmlTree(bool /*createRootEvenIfEmpty*/) const override;

  /** @brief creates an editable combo box. */
  QWidget *createEditor() override {return createComboBox(_isEditable);}
};

/***********************************************/

#endif /* __GROOPSGUI__TREEELEMENTSEQUENCE__ */
