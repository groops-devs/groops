/***********************************************/
/**
* @file programDialog.h
*
* @brief Program List Widget and Dialog.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2015-05-15
*/
/***********************************************/

#ifndef __GROOPSGUI__PROGRAMDIALOG__
#define __GROOPSGUI__PROGRAMDIALOG__

#include <QWidget>
#include <QDialog>
#include <QList>
#include <QTreeWidget>
#include <QLineEdit>
#include "base/importGroops.h"
#include "tree/tree.h"

/***** TYPES ***********************************/

namespace Ui
{
  class ProgramDialog;
}

class QSettings;
class TreeElementProgram;

/***** CLASS ***********************************/

class ProgramListWidget : public QWidget
{
  Q_OBJECT

  class TreeWidget : public QTreeWidget
  {
  public:
    TreeWidget(ProgramListWidget *parent);

    int getItemIndex(QTreeWidgetItem *item);

  private:
    ProgramListWidget *programListWidget;
    QPoint               dragStartPosition;

    void mousePressEvent(QMouseEvent     *event);
    void mouseMoveEvent (QMouseEvent     *event);
  };

  friend class TreeWidget;

  QSettings                 *settings;
  TreeWidget                *treeWidget;
  QLineEdit                 *lineEdit;
  std::vector<XsdElementPtr> programList;
  Tree                      *tree;

public:
  ProgramListWidget(QWidget *parent=nullptr);
  ProgramListWidget(Tree *tree, int indexSelected=0, QWidget *parent=nullptr);
 ~ProgramListWidget();

  void init(Tree *tree, int indexSelected=0);
  int  columnWidth() const { return treeWidget ? treeWidget->columnWidth(0) : 350; }
  void setColumnWidth(int width);
  int  selectedIndex();
  void setUnknownProgram(const QString &name);
  int  programCount() { return programList.size(); }

protected:
  bool eventFilter(QObject *obj, QEvent *event);

signals:
  void programSelected(int index);
  void selectionChanged(bool isProgramSelected);

public slots:
  void itemActivated(QTreeWidgetItem *item, int column);
  void filterChanged(QString text);
  void itemSelectionChanged();
};

/***********************************************/

class ProgramDialog : public QDialog
{
  Q_OBJECT

  Ui::ProgramDialog  *ui;
  QSettings            *settings;
  TreeElementProgram *treeElement;
  ProgramListWidget  *programListWidget;

public:
  ProgramDialog(TreeElementProgram *treeElement, QWidget *parent=nullptr);
 ~ProgramDialog();

public slots:
  void programSelected(int index);
  void accept();

protected:
  void changeEvent(QEvent *e);
};

/***********************************************/

#endif
