/***********************************************/
/**
* @file treeElementChoice.h
*
* @brief Choice element with children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTCHOICE__
#define __GROOPSGUI__TREEELEMENTCHOICE__

#include "base/importGroops.h"
#include "tree/treeElementComplex.h"

/***** CLASS ***********************************/

class TreeElementChoice : public TreeElementComplex
{
  Q_OBJECT

  QStringList schemaNameList;
  QStringList annotationList;

public:
  TreeElementChoice(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride,
                    XmlNodePtr xmlNode, bool fillWithDefaults, bool recieveAutoComments=false);

  /** @brief Is the selection without coresponding schema? */
  bool isSelectionUnknown(int index) const override;

  /** @brief Is the selection renamed (new name in the schema)? */
  bool isSelectionRenamedInSchema(int index) const override;

  /** @brief Generate XML-tree.
  * recursively called for all children. */
  XmlNodePtr createXmlTree(bool /*createRootEvenIfEmpty*/) const override;

  /** @brief changes the current index.
  * calls TreeElemenComplex::selectIndex
  * updates annotation */
  void setSelectedIndex(int index) override;

  /** @brief Is it possible to overweite the element? */
  bool canOverwrite(const QString &type) override {return (this->type() == type);}

  /** @brief Copy the content of @a xmlNode into this.
  * Is undoable.
  * @return success */
  bool overwrite(const QString &type, XmlNodePtr xmlNode, bool contentOnly=false) override;

  // ========================================================

  // Update name of elements in case of schema changes
  // -------------------------------------------------
  /** @brief Is it possible to update name of the element?
  * Only elements with name changes in the schema can be updated. */
  bool canUpdateName() const override;

  /** @brief Update name of the element.
  * Is undoable. */
  void updateName() override;

  // ========================================================

  /** @brief creates an uneditable combo box. */
  QWidget *createEditor() override;

  /** @brief event handler called by item when it gets selected. */
  void startSelected() override;

  /** @brief event handler called by item when it losts the selection. */
  void stopSelected() override;

protected slots:
  void comboBoxHighlighted(int index);
};

/***********************************************/

#endif
