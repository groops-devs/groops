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
                      const QString &defaultOverride, XmlNodePtr xmlNode, bool fillWithDefaults);

  /** @brief Generate XML-tree.
  * recursively called for all children. */
  XmlNodePtr createXmlTree(bool /*createRootEvenIfEmpty*/) const override;

  /** @brief Is it possible to overweite the element? */
  bool canOverwrite(const QString &type) override {return (this->type() == type) && parentElement;}

  /** @brief Copy the content of @a xmlNode into this.
  * Is undoable.
  * @return success */
  bool overwrite(const QString &type, XmlNodePtr xmlNode, bool contentOnly=false) override;

  /** @brief creates an uneditable combo box. */
  QWidget *createEditor() override {return createComboBox(false);}
};

/***********************************************/

#endif
