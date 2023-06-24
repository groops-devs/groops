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
#include "tree/treeElementComplex.h"
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
      item = new TreeItem(treeElement, treeElement->tree);
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

    if(!treeElement->xsdElement)
    {
      icon         = QIcon(":/icons/scalable/element-unknown.svg");
      iconDisabled = QIcon(":/icons/scalable/element-unknown-disabled.svg");
    }
    else if(treeElement->isProgram())
    {
      icon         = QIcon(":/icons/scalable/program.svg");
      iconDisabled = QIcon(":/icons/scalable/program-disabled.svg");
    }
    else if(treeElement->parentElement && treeElement->parentElement->isElementGlobal())
    {
      icon         = QIcon(":/icons/scalable/link.svg");
      iconDisabled = QIcon(":/icons/scalable/link-disabled.svg");
    }
    else if ((!treeElement->optional())  && (!treeElement->unbounded()))
      icon = iconDisabled = QIcon(":/icons/scalable/element-mustset.svg");
    else if ((treeElement->optional())  && (!treeElement->unbounded()))
    {
      icon         = QIcon(":/icons/scalable/element.svg");
      iconDisabled = QIcon(":/icons/scalable/element-disabled.svg");
    }
    else if ((!treeElement->optional())  && (treeElement->unbounded()))
    {
      icon         = QIcon(":/icons/scalable/element-mustset-unbounded.svg");
      iconDisabled = QIcon(":/icons/scalable/element-mustset-unbounded-disabled.svg");
    }
    else if ((treeElement->optional())  && (treeElement->unbounded()))
    {
      icon         = QIcon(":/icons/scalable/element-unbounded.svg");
      iconDisabled = QIcon(":/icons/scalable/element-unbounded-disabled.svg");
    }

    updateDisabled();
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
  try
  {
    lostCurrent();
  }
  catch(std::exception &e)
  {
    qDebug() << QString::fromStdString("Exception in destructor at "+_GROOPS_ERRORLINE+"\n"+e.what());
  }
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

void TreeItem::updateDisabled()
{
  QIcon baseIcon = (treeElement()->disabled()) ? iconDisabled : icon;

  if(treeElement()->isRenamed()) // add renamed icon to base icon
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
  QString loop      = treeElement()->loop();
  QString condition = treeElement()->condition();
  QString text      = treeElement()->name();

  if(treeElement()->isRenamed() && !treeElement()->parentElement->isElementGlobal())
    text = treeElement()->originalName();
  if(!loop.isEmpty())
    text += " [loop=" + loop + "]";
  if(!condition.isEmpty())
    text += " [condition=" + condition + "]";

  setText(0, text);
}

/***********************************************/

void TreeItem::updateValue()
{
  if(valueEditor)
    return;

  QString text = treeElement()->selectedValue();
  QString result = treeElement()->selectedResult();
  if(treeElement()->tree->showResults() && (!result.isEmpty()) && result != text && result != "{"+text+"}")
    text += "   [="+result+"]";
  setText(1, text);
  QIcon icon;
  if(treeElement()->isLinked())
    icon = QIcon(":/icons/scalable/link.svg");
  else if(treeElement()->isSelectionUnknown(treeElement()->selectedIndex()))
    icon = QIcon(":/icons/scalable/element-unknown.svg");
  else if(treeElement()->isSelectionRenamed(treeElement()->selectedIndex()))
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
    treeElement()->tree->scrollToItem(this, QAbstractItemView::EnsureVisible);

    if(!valueEditor)
      valueEditor = treeElement()->createEditor();

    if(valueEditor)
    {
      // clear cell
      setText(1, "");
      setIcon(1, QIcon());
      // init widget
      valueEditor->setContentsMargins(0,0,0,0);
      valueEditor->setFixedHeight(valueEditor->sizeHint().height());
      setSizeHint(1, valueEditor->sizeHint());
      treeElement()->tree->setItemWidget(this, 1, valueEditor);
      treeElement()->comboBoxSetToolTip();
      valueEditor->show();
      treeElement()->tree->update();
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
      valueEditor->hide();
      valueEditor->setFixedHeight(1);
      treeElement()->tree->removeItemWidget(this, 1);
      //delete valueEditor;
      valueEditor = nullptr;
      updateValue();
      setSizeHint(1, QSize());
      treeElement()->tree->update();
    }
    if(commentEditor)
    {
      treeElement()->tree->removeItemWidget(this, 3);
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
    if(commentEditor)
      delete commentEditor;
    setText(3, "");
    commentEditor = new QLineEdit(treeElement()->comment());
    connect(commentEditor, SIGNAL(editingFinished()), this, SLOT(editCommentFinished()));
    treeElement()->tree->setItemWidget(this, 3, commentEditor);
    commentEditor->setFocus();
    commentEditor->selectAll();

/*QKeyEventTransition *pressedEsc = new QKeyEventTransition(edit, QEvent::KeyPress, Qt::Key_Escape);
connect(pressedEsc, SIGNAL(triggered()), edit, SLOT(clear()));  }*/
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
    treeElement()->tree->removeItemWidget(this, 3);
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
