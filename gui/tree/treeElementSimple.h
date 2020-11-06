/***********************************************/
/**
* @file treeElementSimple.h
*
* @brief Element without children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTSIMPLE__
#define __GROOPSGUI__TREEELEMENTSIMPLE__

#include "base/importGroops.h"
#include "tree/treeElement.h"

/***** CLASS ***********************************/

class TreeElementSimple : public TreeElement
{
  Q_OBJECT

public:
  TreeElementSimple(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                    const QString &defaultOverride, XmlNodePtr xmlNode, Bool fromFile);
  virtual ~TreeElementSimple() override {}

/** @brief Generate XML-tree. */
virtual XmlNodePtr getXML(Bool withEmptyNodes=false) const override;

/** @brief Values can be edited. */
virtual Bool isEditable() const override {return true;}

/** @brief creates an editable combo box. */
virtual QWidget *createEditor() override {return createComboBox(true);}
};

/***********************************************/

#endif
