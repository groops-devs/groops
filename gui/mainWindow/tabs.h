/***********************************************/
/**
* @file tabs.h
*
* @brief Tab environment.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2011-07-16
*/
/***********************************************/

#ifndef __GROOPSGUI__TABS__
#define __GROOPSGUI__TABS__

#include <QTabWidget>
#include <QTabBar>
#include <QFrame>
#include <QLabel>
#include <QSettings>
#include "base/importGroops.h"
#include "mainWindow/mainWindow.h"

/***** TYPES ***********************************/

class Tree;

/***** CLASS ***********************************/

class TabEnvironment : public QTabWidget
{
  Q_OBJECT

  ActionList  actionList;
  QSettings   settings;
  QUndoGroup *undoGroup;
  int         newFileCounter;
  bool        isFullyLoaded;

  friend class TabBar;

public:
  TabEnvironment(QWidget *parent, ActionList *actionList, QUndoGroup *undoGroup);
 ~TabEnvironment();

  /** @brief Get the current selected tree. */
  Tree *currentTree() const;

  /** @brief Makes tree the current tree. */
  void setCurrentTree(Tree *tree);

  /** @brief List of all trees. */
  QVector<Tree*> trees() const;

private:
  Tree *treeAt(int index) const;

  /** @brief Creates a new tab.
  * Opens file if @a fileName is given.
  * @return success? */
  bool newTab(const QString &fileName);

public slots:
  /** @brief Creates a new empty tab.
  * @return success? */
  bool fileNew() {return newTab(QString());}

  /** @brief Open the file with @a fileName.
  * If @a fileName is empty a file selector is opened.
  * The user is asked to save the old file before if necessary.
  * @return success? */
  bool fileOpen(const QString &fileName = QString());

  /** @brief Reload the file.
  * The user is asked to save the file before if necessary.
  * @return success? */
  bool fileReOpen();

  /** @brief Close the current tab.
  * The user is asked to save the file before if necessary.
  * @return success? */
  bool fileClose();

  /** @brief Close all but the current tab.
  * The user is asked to save the other files before if necessary. */
  void fileCloseOther();

  /** @brief Is it allowed to close the file?
  * The user is asked to save the file before if necessary.
  * @return allowed */
  bool okToAbandon();

signals:
  /** @brief This signal is emitted when tab is changed.
  * @param caption fileName without path or something like 'newX'.
  * @param fileName the current file (empty if file is new)
  * @param isClean no unsaved changes */
  void fileChanged(const QString &caption, const QString &fileName, bool isClean);

protected:
  void showEvent(QShowEvent *event) override;
  void resizeEvent(QResizeEvent *event) override;

private slots:
  void treeSectionResized(int logicalIndex, int /*oldSize*/, int newSize);
  void treeFileChanged(const QString &/*caption*/, const QString &/*fileName*/, bool /*isClean*/);
  void tabClose(int index);
  void tabChanged(int index);
};

/***********************************************/

// new tabBar implementation to implement context menu and double click
class TabBar : public QTabBar
{
  Q_OBJECT

public:
  TabBar(TabEnvironment *tabEnvironment);

public slots:
  void contextMenuEvent(QContextMenuEvent *e);
  void mouseDoubleClickEvent(QMouseEvent *e);

private:
  TabEnvironment *tabEnvironment;
  QPoint dragStartPosition;

  void dragEnterEvent(QDragEnterEvent *event);
  void dragMoveEvent(QDragMoveEvent *event);
  void dropEvent(QDropEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
};

/***********************************************/

#endif
