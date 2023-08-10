/***********************************************/
/**
 * @file treeElementComment.cpp
 *
 * @brief Comment element.
 *
 * @author Torsten Mayer-Guerr
 * @date 2023-06-21
 */
/***********************************************/

#include <QtDebug>
#include <QTextEdit>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementGlobal.h"
#include "tree/treeItem.h"
#include "tree/treeElementComment.h"

/***********************************************/

TreeElementComment::TreeElementComment(Tree *tree, TreeElementComplex *parentElement, const QString &text)
  : TreeElement(tree, parentElement, nullptr, QString(), nullptr)
{
  try
  {
    insertNewValue(text, false);
    setSelectedIndex(0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlNodePtr TreeElementComment::createXmlTree(bool /*createRootEvenIfEmpty*/) const
{
  try
  {
    if(selectedValue().isEmpty())
       return XmlNodePtr(nullptr);
    XmlNodePtr xmlNode = XmlNode::create("COMMENT");
    xmlNode->setText(selectedValue());
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComment::setSelectedIndex(int index)
{
  TreeElement::setSelectedIndex(index);
  if(textEditor && (selectedValue() != textEditor->toPlainText()))
    textEditor->setPlainText(selectedValue());
}

/***********************************************/

bool TreeElementComment::overwrite(const QString &type, XmlNodePtr xmlNode, bool /*contentOnly*/)
{
  try
  {
    if(!canOverwrite(type) || !xmlNode)
      return false;

    tree->undoStack->beginMacro("overwrite "+name());
    changeSelectedValue(xmlNode->getText());
    tree->undoStack->endMacro();
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

QWidget *TreeElementComment::createEditor()
{
  try
  {
    textEditor = new QTextEdit(tree);
    textEditor->setContentsMargins(0,0,0,0);
    textEditor->setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed));
    textEditor->setWordWrapMode(QTextOption::NoWrap);
    textEditor->setTabChangesFocus(true);
    textEditor->setTextInteractionFlags(Qt::TextEditorInteraction);
    textEditor->setUndoRedoEnabled(false); // own implementation
    textEditor->setFocusPolicy(Qt::StrongFocus);
    textEditor->installEventFilter(tree);
    textEditor->setPlainText(selectedValue());

    QFontMetrics fontMetrics(textEditor->document()->defaultFont());
    textEditor->setFixedHeight(fontMetrics.size(0, textEditor->toPlainText()+"\n\nextraline").height() + 2*textEditor->frameWidth());

    connect(textEditor, SIGNAL(textChanged()), this, SLOT(textEditorTextChanged()));
    return textEditor;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComment::textEditorTextChanged()
{
  try
  {
    if(isLinked() || (textEditor->toPlainText() != selectedValue()))
      tree->undoStack->push(new UndoCommandEdit(this, textEditor->toPlainText()));

    // adjust size
    QFontMetrics fontMetrics(textEditor->document()->defaultFont());
    textEditor->setFixedHeight(fontMetrics.size(0, textEditor->toPlainText()+"\n\nextraline").height() + 2*textEditor->frameWidth());
    if(item())
      item()->setSizeHint(1, textEditor->size());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
