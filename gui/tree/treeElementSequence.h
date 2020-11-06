/***********************************************/
/**
* @file treeElementSequence.h
*
* @brief Element with children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTSEQUENCE__
#define __GROOPSGUI__TREEELEMENTSEQUENCE__

#include "base/importGroops.h"
#include "tree/treeElementComplex.h"

/***** CLASS ***********************************/

class TreeElementSequence : public TreeElementComplex
{
public:
  TreeElementSequence(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                      const QString &defaultOverride, XmlNodePtr xmlNode, Bool fromFile);
  virtual ~TreeElementSequence() override {}

/** @brief Generate XML-tree.
* recursively called for all children. */
virtual XmlNodePtr getXML(Bool withEmptyNodes=false) const override;

/** @brief creates an uneditable combo box. */
virtual QWidget *createEditor() override {return createComboBox(false);}
};

/***********************************************/

#endif
