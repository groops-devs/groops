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

  Bool _isEditable;

public:
  TreeElementUnknown(Tree *tree, TreeElementComplex *parentElement, XmlNodePtr xmlNode);
  virtual ~TreeElementUnknown() override {}

/** @brief Values can be edited. */
virtual Bool isEditable() const override {return _isEditable;}

/** @brief Generate XML-tree.
* recursivly called for all children. */
virtual XmlNodePtr getXML(Bool withEmptyNodes=false) const override;

/** @brief creates an editable combo box. */
virtual QWidget *createEditor() override {return createComboBox(_isEditable);}
};

/***********************************************/

#endif /* __GROOPSGUI__TREEELEMENTSEQUENCE__ */
