/***********************************************/
/**
* @file sideBar.h
*
* @brief Side bar and all of its components.
*
* @author Sebastian Strasser
* @date 2017-06-17
*/
/***********************************************/

#ifndef __GROOPSGUI__SIDEBAR__
#define __GROOPSGUI__SIDEBAR__

#include <QStackedWidget>
#include <QPushButton>
#include <QButtonGroup>
#include <QPaintEvent>
#include <QVBoxLayout>
#include <QTreeWidget>
#include "base/importGroops.h"
#include "mainWindow/tabs.h"
#include "mainWindow/mainWindow.h"

/***** CLASS ***********************************/

/** @brief Side bar containing the buttons to switch side bar widgets, but not the side bar widgets themselves. */
class SideBar : public QWidget
{
  Q_OBJECT

  QStackedWidget       *_stackedWidget;
  QButtonGroup          buttonGroup;
  QVBoxLayout          *layout;
  QSettings             settings;
  int                   _lastWidth;

public:
  SideBar(QWidget *parent);
 ~SideBar();

  /** @brief Adds a @a button to the side bar and a @a widget to the side bar widgets stack. */
  void addSideBarWidget(const QString &buttonLabel, QWidget *widget);

  /** @brief Returns a pointer to the side bar widgets stack. */
  QStackedWidget *stackedWidget() const {return _stackedWidget;}

  /** @brief Returns the last set width of the widget. */
  int lastWidth() {return _lastWidth;}

public slots:
  /** @brief Updates the side bar selection and shows/hides the side bar widget. */
  void buttonClicked(int id);

  /** @brief Set last width of the widget, used for saving settings. */
  void widthChanged(int width, int index);
};

/***********************************************/

/** @brief PushButton that allows for vertical orientation. */
class PushButtonVertical : public QPushButton
{
  Q_OBJECT

public:
  PushButtonVertical(QString text, QWidget *parent=nullptr) : QPushButton(text, parent) {};

  QSize minimumSizeHint() const override {return QPushButton::minimumSizeHint().transposed();}
  QSize sizeHint() const override  {return QPushButton::sizeHint().transposed();}

protected:
  void  paintEvent(QPaintEvent *event) override; // Reimplemented to support vertical orientation.
};

/***********************************************/

/** @brief "Open Files" side bar widget with intelligent tree structure. */
class OpenFilesTreeWidget : public QTreeWidget
{
  Q_OBJECT

  ActionList      actionList;
  TabEnvironment *tabEnvironment;

public:
  /** @brief Initializes the widget and populates the tree with items from the @a tabEnvironment. */
  OpenFilesTreeWidget(TabEnvironment *tabEnvironment, ActionList &actionList);

  /** @brief Populates the tree with items from the tabEnvironment. */
  void populateTree();

  /** @brief Recursively check if an @a item has a "file" as child, grandchild, ...
   *  @return True if a "file" is found. */
  bool hasFileChild(const QTreeWidgetItem *item) const;

public slots:
  /** @brief Emits @a fileSelectionChanged(int) signal if another item is selected in the tree. */
  void currentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem *previous);

  /** @brief Updates the tree if anything changes (file opened/changed/closed/...). */
  void fileChanged(const QString &caption, const QString &fileName, bool isClean);

protected:
  /** @brief Catch clicks on "folders" to collapse/expand them instead of selecting them. */
  void mousePressEvent(QMouseEvent *event);

  /** @brief Create context menu if a file is right-clicked. */
  void contextMenuEvent(QContextMenuEvent *e);

  void dragEnterEvent (QDragEnterEvent *event);
  void dragMoveEvent  (QDragMoveEvent  *event);
  void dropEvent      (QDropEvent      *event);

signals:
  /** @brief This signal is emitted if another file is selected in the tree. */
  void fileSelectionChanged(int index);
};

/***********************************************/

/** @brief TreeWidgetItem that sorts "folders" above "files". */
class TreeWidgetItem : public QTreeWidgetItem
{
public:
  TreeWidgetItem(QTreeWidgetItem *parent=nullptr) : QTreeWidgetItem(parent) {}
  virtual bool operator<(const QTreeWidgetItem &other) const;
};

/***********************************************/

#endif // __GROOPSGUI__SIDEBAR__
