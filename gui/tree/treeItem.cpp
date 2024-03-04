/***********************************************/
/**
* @file treeItem.cpp
*
* @brief Visible representation of an element as an item in the tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#include <QtDebug>
#include <QLineEdit>
#include <QPainter>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementAdd.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementGlobal.h"
#include "tree/treeElementProgram.h"
#include "tree/treeElementLoopCondition.h"
#include "tree/treeElementComment.h"
#include "tree/treeElementUnknown.h"
#include "tree/treeItem.h"

/***********************************************/
/***********************************************/

TreeItem *TreeItem::newTreeItem(TreeElement *treeElement, TreeItem *parent, TreeItem *after)
{
  try
  {
    if(!treeElement)
      throw(Exception("TreeElement==nullptr"));
    TreeItem *item = nullptr;
    if(!parent)
      item = new TreeItem(treeElement, treeElement->tree->treeWidget);
    else if(after)
      item = new TreeItem(treeElement, parent, after);
    else
      item = new TreeItem(treeElement, parent);

    item->setFocus();

    return item;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeItem::init(TreeElement *treeElement)
{
  try
  {
    this->valueEditor   = nullptr;
    this->_treeElement  = treeElement;
    this->commentEditor = nullptr;

    if(dynamic_cast<TreeElementComment*>(treeElement))
      icon = iconDisabled = QIcon(":/icons/scalable/element-comment.svg");
    else if(dynamic_cast<TreeElementLoopCondition*>(treeElement) && (treeElement->type() == "loopType"))
      icon = iconDisabled = QIcon(":/icons/scalable/loop-set.svg");
    else if(dynamic_cast<TreeElementLoopCondition*>(treeElement) && (treeElement->type() == "conditionType"))
      icon = iconDisabled = QIcon(":/icons/scalable/condition-set.svg");
    else if(!treeElement->label().isEmpty())
    {
      icon         = QIcon(":/icons/scalable/link.svg");
      iconDisabled = QIcon(":/icons/scalable/link-disabled.svg");
    }
    else if(dynamic_cast<TreeElementProgram*>(treeElement))
    {
      icon         = QIcon(":/icons/scalable/program.svg");
      iconDisabled = QIcon(":/icons/scalable/program-disabled.svg");
    }
    else if(dynamic_cast<TreeElementAdd*>(treeElement))
      icon = iconDisabled = QIcon(":/icons/scalable/element-unbounded.svg");
    else if(dynamic_cast<TreeElementUnknown*>(treeElement))
    {
      icon         = QIcon(":/icons/scalable/element-unknown.svg");
      iconDisabled = QIcon(":/icons/scalable/element-unknown-disabled.svg");
    }
    else if(!treeElement->optional() && !treeElement->unbounded())
      icon = iconDisabled = QIcon(":/icons/scalable/element-mustset.svg");
    else if ((treeElement->optional()) && !treeElement->unbounded())
    {
      icon         = QIcon(":/icons/scalable/element.svg");
      iconDisabled = QIcon(":/icons/scalable/element-disabled.svg");
    }
    else if(!treeElement->optional() && treeElement->unbounded())
    {
      icon         = QIcon(":/icons/scalable/element-mustset-unbounded.svg");
      iconDisabled = QIcon(":/icons/scalable/element-mustset-unbounded-disabled.svg");
    }
    else if(treeElement->optional() && treeElement->unbounded())
    {
      icon         = QIcon(":/icons/scalable/element-unbounded.svg");
      iconDisabled = QIcon(":/icons/scalable/element-unbounded-disabled.svg");
    }

    updateIcon();
    updateName();
    updateValue();
    updateAnnotation(treeElement->annotation());
    updateComment();

    setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeItem::~TreeItem()
{
  lostCurrent();
}

/***********************************************/

QVariant TreeItem::data(int column, int role) const
{
  if(role == Qt::ToolTipRole)
  {
    QFontMetrics fm(font(column));
    if(treeWidget()->columnWidth(column) < fm.boundingRect(text(column)).width()+7)
      return QVariant(text(column).replace(" [=", "\n[="));
  }

  return QTreeWidgetItem::data(column, role);
}

/***********************************************/

void TreeItem::updateIcon()
{
  QIcon baseIcon = (treeElement()->disabled()) ? iconDisabled : icon;

  if(treeElement()->isRenamedInSchema()) // add renamed icon to base icon
  {
    int iconHeight = treeWidget()->iconSize().height();
    QIcon renamedIcon(":/icons/scalable/edit-rename.svg");
    QPixmap pixmap(2*iconHeight, iconHeight);
    pixmap.fill(QColor(0,0,0,0));
    QPainter painter(&pixmap);
    painter.drawPixmap(0, 0, baseIcon.pixmap(iconHeight));
    painter.drawPixmap(iconHeight+2, 0, renamedIcon.pixmap(iconHeight));
    setIcon(0, pixmap);
  }
  else
    setIcon(0, baseIcon);
}

/***********************************************/

void TreeItem::updateName()
{
  QString text = treeElement()->name();
  if(treeElement()->loop)      text += " [loop]";
  if(treeElement()->condition) text += " [condition]";
  setText(0, text);
}

/***********************************************/

void TreeItem::updateValue()
{
  if(valueEditor)
    return;

  QString text   = treeElement()->selectedValue();
  QString result = treeElement()->selectedResult();
  if(treeElement()->tree->showResults() && !result.isEmpty() && (result != text) && (result != "{"+text+"}"))
    text += "   [="+result+"]";
  setText(1, text);
  QIcon icon;
  if(treeElement()->isBrokenLinked())
    icon = QIcon(":/icons/scalable/link-broken.svg");
  else if(treeElement()->isLinked())
    icon = QIcon(":/icons/scalable/link.svg");
  else if(dynamic_cast<TreeElementChoice*>(treeElement()) && dynamic_cast<TreeElementChoice*>(treeElement())->isSelectionUnknown(treeElement()->selectedIndex()))
    icon = QIcon(":/icons/scalable/element-unknown.svg");
  else if(dynamic_cast<TreeElementChoice*>(treeElement()) && dynamic_cast<TreeElementChoice*>(treeElement())->isSelectionRenamedInSchema(treeElement()->selectedIndex()))
    icon = QIcon(":/icons/scalable/edit-rename.svg");
  setIcon(1, icon);
}

/***********************************************/

void TreeItem::updateAnnotation(const QString &text)
{
  setText(2, text);
}

/***********************************************/

void TreeItem::updateComment()
{
  if(commentEditor)
    return;
  setText(3, treeElement()->comment());
}

/***********************************************/

void TreeItem::becomeCurrent()
{
  try
  {
    treeWidget()->scrollToItem(this, QAbstractItemView::EnsureVisible);

    if(!valueEditor)
      valueEditor = treeElement()->createEditor();

    if(valueEditor)
    {
      // clear cell
      setText(1, "");
      setIcon(1, QIcon());
      // init widget
      treeWidget()->setItemWidget(this, 1, valueEditor);
      // treeElement()->tree->update();
    }

    treeElement()->startSelected();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeItem::lostCurrent()
{
  try
  {
    treeElement()->stopSelected();

    if(valueEditor)
    {
      treeWidget()->removeItemWidget(this, 1);
      valueEditor = nullptr;
      updateValue();
      setSizeHint(1, QSize());
      // treeElement()->tree->update();
    }
    if(commentEditor)
    {
      treeWidget()->removeItemWidget(this, 3);
      delete commentEditor;
      commentEditor = nullptr;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeItem::setFocus()
{
  try
  {
    if(valueEditor)
      valueEditor->setFocus();
    else
      treeWidget()->setFocus();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeItem::setSelection(int start, int length)
{
  try
  {
    if(!valueEditor || start < 0)
      return;

    QLineEdit *lineEdit = valueEditor->findChild<QLineEdit*>();
    if(lineEdit)
      lineEdit->setSelection(start, length);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeItem::selection(int &start, int &length) const
{
  try
  {
    start = -1;
    length = 0;
    if(!valueEditor)
      return;

    QLineEdit *lineEdit = valueEditor->findChild<QLineEdit*>();
    if(lineEdit)
    {
      start = lineEdit->selectionStart() >= 0 ? lineEdit->selectionStart() : lineEdit->cursorPosition();
      length = lineEdit->selectedText().size();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeItem::editComment()
{
  try
  {
    if(!treeElement()->canComment())
      return;
    setText(3, "");
    commentEditor = new QLineEdit(treeElement()->comment());
    connect(commentEditor, SIGNAL(editingFinished()), this, SLOT(editCommentFinished()));
    treeWidget()->setItemWidget(this, 3, commentEditor);
    commentEditor->setFocus();
    commentEditor->selectAll();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeItem::editCommentFinished()
{
  try
  {
    if(!commentEditor)
      return;
    treeElement()->setComment(commentEditor->text());
    treeWidget()->removeItemWidget(this, 3);
    commentEditor->deleteLater();
    commentEditor = nullptr;
    updateComment();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
