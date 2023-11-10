/***********************************************/
/**
* @file treeItem.h
*
* @brief Visible representation of an element as an item in the tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEITEM__
#define __GROOPSGUI__TREEITEM__

#include <QTreeWidgetItem>
#include <QPointer>
#include "base/importGroops.h"

/***** TYPES ***********************************/

class  TreeElement;
class  QLineEdit;

/***** CLASS ***********************************/

class TreeItem : public QObject, public QTreeWidgetItem
{
  Q_OBJECT

  TreeElement      *_treeElement;
  QIcon              icon, iconDisabled;
  QPointer<QWidget>  valueEditor;
  QLineEdit         *commentEditor;

  void init(TreeElement *treeElement);

  TreeItem(TreeElement *treeElement, QTreeWidget *parent)                : QTreeWidgetItem(parent)        {init(treeElement);}
  TreeItem(TreeElement *treeElement, TreeItem  *parent)                  : QTreeWidgetItem()              {parent->insertChild(0,this); init(treeElement);}
  TreeItem(TreeElement *treeElement, TreeItem  *parent, TreeItem *after) : QTreeWidgetItem(parent, after) {init(treeElement);}

  // adjust tooltip text
  QVariant data(int column, int role) const override;

public:
  virtual ~TreeItem();

  /** @brief Append a new item to the tree.
  * If parent==nullptr the element will be appended to root. */
  static TreeItem *newTreeItem(TreeElement *treeElement, TreeItem *parent=nullptr, TreeItem *after=nullptr);

  /** @brief The corresponding tree element. */
  TreeElement *treeElement() const {return _treeElement;}

  void updateIcon();
  void updateName();
  void updateValue();
  void updateAnnotation(const QString &text);
  void updateComment();

  void becomeCurrent();
  void lostCurrent();

  void setFocus();

  // for and search/replace
  void setSelection(int start, int length);
  void selection(int &start, int &length) const;

  void editComment();

public slots:
  void editCommentFinished();
};

/***********************************************/

#endif
