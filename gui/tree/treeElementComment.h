/***********************************************/
/**
* @file treeElementComment.h
*
* @brief Comment element.
*
* @author Torsten Mayer-Guerr
* @date 2023-06-21
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTCOMMENT__
#define __GROOPSGUI__TREEELEMENTCOMMENT__

#include <QPointer>
#include "base/importGroops.h"
#include "tree/treeElement.h"

/***** TYPES ***********************************/

class QTextEdit;

/***** CLASS ***********************************/

class TreeElementComment : public TreeElement
{
  Q_OBJECT

  QPointer<QTextEdit> textEditor;

public:
  TreeElementComment(Tree *tree, TreeElementComplex *parentElement, const QString &text);
 ~TreeElementComment() {}

  QString type()              const override {return "COMMENT";}
  bool    optional()          const override {return true;}
  bool    unbounded()         const override {return true;}
  bool    isRenamedInSchema() const override {return false;}

  /** @brief Generate XML-tree. */
  XmlNodePtr createXmlTree(bool /*createRootEvenIfEmpty*/) const override;

  /** @brief Values can be edited. */
  bool isEditable() const override {return true;}

  /** @brief Updates editor. */
  void setSelectedIndex(int index) override;

   /** @brief Is it possible to overweite the element? */
  bool canOverwrite(const QString &type) override {return (this->type() == type);}

  /** @brief Copy the content of @a xmlNode into this.
  * Is undoable.
  * @return success */
  bool overwrite(const QString &type, XmlNodePtr xmlNode, bool contentOnly=false) override;

  bool canSetLoop()      const override {return false;}
  bool canSetCondition() const override {return false;}
  bool canDisabled()     const override {return false;}
  bool canComment()      const override {return false;}

  /** @brief creates an editable combo box. */
  QWidget *createEditor() override;

private slots:
  void textEditorTextChanged();
};

/***********************************************/

#endif
