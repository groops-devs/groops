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
#include <QPaintEvent>
#include <QVBoxLayout>
#include <QTreeWidget>
#include "base/importGroops.h"
#include "mainWindow/tabs.h"
#include "mainWindow/mainWindow.h"

/***** CLASS ***********************************/

class PushButton;
class StackedWidget;

/** @brief Side bar containing the buttons to switch side bar widgets, but not the side bar widgets themselves. */
class SideBar : public QWidget
{
  Q_OBJECT

public:
  SideBar(QWidget *parent);
 ~SideBar();

  /** @brief Adds a @a button to the side bar and a @a widget to the side bar widgets stack. */
  void addSideBarWidget(QString buttonLabel, QWidget *widget);

  /** @brief Returns a pointer to the side bar widgets stack. */
  StackedWidget* stackedWidget() const { return _stackedWidget; }

  /** @brief Show or @a hide the side bar widgets. */
  void setSideBarWidgetsHidden(bool hide);

private:
  /** @brief Adds a new button with @a text and horizontal or vertical @a orientation to the side bar.
   *  @return Pointer to the button. */
  PushButton* addPushButton(QString text, Qt::Orientation orientation = Qt::Vertical);

  StackedWidget            *_stackedWidget;
  std::vector<QPushButton*> buttons;
  QVBoxLayout              *layout;
  QSettings                *settings;
  int                       currentIndex;

signals:
  /** @brief This signal is emitted when a button on the side bar is pressed.
   *  @param isCurrent True if the currently selected button was pressed again, false if another button was pressed. */
  void showHide(bool isCurrent = false);

  /** @brief This signal is emitted when the "page" of the side bar is changed.
   *  @param newIndex Index of the newly selected "page". */
  void currentChanged(int newIndex = 0);

public slots:
  /** @brief Updates the side bar selection and shows/hides the side bar widget. */
  void buttonClicked(bool checked = false);

  /** @brief Updates the side bar selection and shows the side bar widget. */
  void showWidget();
};

/***********************************************/

/** @brief Stacked widget containing the side bar widgets. */
class StackedWidget : public QStackedWidget
{
  Q_OBJECT

public:
  StackedWidget(int lastWidth = 200) : _lastWidth(lastWidth) {}

  /** @brief Returns the last set width of the widget. */
  int lastWidth() { return _lastWidth; }

public slots:
  /** @brief Shows side bar widget if it is hidden or hides it if @a isCurrent = true (active side bar button pressed again). */
  void showHide(bool isCurrent = false);

  /** @brief Set last width of the widget, used for saving settings. */
  void widthChanged(int width, int index);

private:
  int _lastWidth;
};

/***********************************************/

/** @brief PushButton that allows for vertical orientation. */
class PushButton : public QPushButton
{
  Q_OBJECT

public:
  PushButton(QString text, QWidget *parent=nullptr);

  /** @brief Set the @a orientation of the button to vertical or horizontal. */
  void setOrientation(Qt::Orientation orientation);

  /** @brief Returns the @a orientation of the button. */
  Qt::Orientation orientation() const;

  /** @brief Returns minimum @a size for the button. Reimplemented to support vertical orientation. */
  QSize minimumSizeHint() const;

  /** @brief Returns recommended @a size for the button. Reimplemented to support vertical orientation. */
  QSize sizeHint() const;

protected:
  /** @brief Paints the button. Reimplemented to support vertical orientation. */
  void  paintEvent(QPaintEvent *event);

private:
  Qt::Orientation _orientation;
};

/***********************************************/

/** @brief "Open Files" side bar widget with intelligent tree structure. */
class OpenFilesTreeWidget : public QTreeWidget
{
  Q_OBJECT

public:
  /** @brief Initializes the widget and populates the tree with items from the @a workspace. */
  void init(TabEnvironment *workspace, ActionList *actionList);

  /** @brief Populates the tree with items from the workspace. */
  void populateTree();

  /** @brief Recursively check if an @a item has a "file" as child, grandchild, ...
   *  @return True if a "file" is found. */
  bool hasFileChild(const QTreeWidgetItem *item) const;

public slots:
  /** @brief Emits @a fileSelectionChanged(int) signal if another item is selected in the tree. */
  void currentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem *previous);

  /** @brief Updates the tree if anything changes (file opened/changed/closed/...). */
  void fileChanged(const QString &fileName, bool changed);

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

private:
  ActionList      actionList;
  TabEnvironment *workspace;

  QTreeWidgetItem *updateItemStatus(const QString &fileName, bool changed);
};

/***********************************************/

/** @brief TreeWidgetItem that sorts "folders" above "files". */
class TreeWidgetItem : public QTreeWidgetItem
{
public:
  TreeWidgetItem(QTreeWidgetItem* parent=nullptr) : QTreeWidgetItem(parent) {}

private:
  virtual bool operator < (const QTreeWidgetItem &other) const;
};

/***********************************************/

#endif // __GROOPSGUI__SIDEBAR__
