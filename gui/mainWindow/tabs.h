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
#include "base/importGroops.h"
#include "mainWindow/mainWindow.h"

/***** TYPES ***********************************/

class QSettings;
class Tree;

/***** CLASS ***********************************/

class TabEnvironment : public QTabWidget
{
  Q_OBJECT

friend class TabBar;

public:
  TabEnvironment(QWidget *parent, ActionList *actionList, QUndoGroup *undoGroup);
 ~TabEnvironment();

/**
* @brief Get the current selected tree.
*/
Tree *currentTree() const;

/**
 * @brief Get the file names of all opened and newly created files (with or without an associated XML file).
 */
QStringList allFileNames() const;

/**
 * @brief Get the file names of all opened files with an associated XML file.
 */
QStringList openedFileNames() const;

/**
 * @brief Get the names of all newly created files without an associated XML file.
 */
QStringList newFileNames() const;

/**
 * @brief Get the text of the currently active tab.
 */
QString currentTabText() const;

/**
 * @brief Resize columns of all trees.
 */
void resizeTreeColumns(const std::vector<int> &columnWidths);

public slots:
/**
* @brief show/hide column with the annotations depending on state.
*/
void setShowDescriptions(Bool state);

/**
* @brief show/hide results of links and expressions depending on state.
*/
void setShowResults(Bool state);

/**
* @brief Creates a new empty tab.
*
* @return success?
*/
Bool newFile();

/**
* @brief Close the current tab.
*
* The user is asked to save the file before if necessary.
* @return success?
*/
Bool closeFile();

/**
* @brief Close all but the current tab.
*
* The user is asked to save the other files before if necessary.
* @return success?
*/
Bool closeOtherFiles();

/**
* @brief Reload the file.
*
* The user is asked to save the file before if necessary.
* @return success?
*/
Bool reopenFile();

/**
* @brief Show file in external file manager.
*/
void showFileInManager() const;

/**
* @brief Open the file with @a fileName.
*
* If @a fileName is empty a file selector is opened.
* The user is asked to save the old file before if necessary.
* @return success?
*/
Bool openFile(const QString &fileName = QString(), bool emitFileOpened = true);

/**
* @brief Open the files from the previous session (from settings).
*/
void openInitialFiles();

/**
* @brief Save the file.
*
* A file selector is opened for trees without file name.
* @return success?
*/
Bool saveFile();

/**
* @brief Save the file with a new @a fileName.
*
* If @a fileName is empty a file selector is opened.
* @return success?
*/
Bool saveAsFile(const QString &fileName = QString());

/**
* @brief Exceute file.
*
* opens the execute dialog.
* @return success?
*/
Bool execFile();

public:
/**
* @brief Is it allowed to close the file?
*
* The user is asked to save the file before if necessary.
* @return allowed
*/
Bool okToAbandon();

/**
* @brief URL to documentation of currently selected element.
* @return URL to documentation file, or empty URL of no documentation exists.
*/
QUrl documentationUrl();

signals:
/**
* @brief This signal is emitted when the file is changed.
*
* @param fileName the current file (empty if file is new)
* @param changed current file has unsaved changes
*/
void fileChanged(const QString &fileName, bool changed);

/**
* @brief This signal is emitted when a new file is opened.
*
* @param fileName the current file (empty if file is new)
*/
void fileOpened(const QString &fileName);

/**
* @brief This signal is emitted when a selection in a tree is changed.
*
* @param selection Name of the selected item or empty if no selection.
*/
void treeSelectionChanged(const QString &selection);

/**
* @brief This signal is emitted when the documentation for a selection is requested.
*
* @param selection Name of the selected item or empty if no selection.
*/
void loadDocumentation(const QString &selection);

/**
* @brief This signal is emitted when the schema is changed.
*/
void schemaChanged();

private:
  ActionList  actionList;
  QSettings  *settings;
  QUndoGroup *undoGroup;
  int         newFileCounter;
  bool        isFullyLoaded;

  void saveOpenFiles();
  Tree *treeAt(int index) const;

protected:
  void showEvent(QShowEvent *event);

private slots:
  void tabClose(int index);
  void tabChanged(int index);
  void treeChanged(const QString &fileName, bool changed);
};

/***********************************************/

// new tabBar implementation to implement context menu and double click
class TabBar : public QTabBar
{
  Q_OBJECT

public:
  TabBar(TabEnvironment *tabEnvironment, QWidget *parent);
 ~TabBar();

public slots:
  void contextMenuEvent(QContextMenuEvent *e);
  void mouseDoubleClickEvent(QMouseEvent *e);

private:
  TabEnvironment *tabs;
  QPoint dragStartPosition;

  void dragEnterEvent(QDragEnterEvent *event);
  void dragMoveEvent(QDragMoveEvent *event);
  void dropEvent(QDropEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
};

/***********************************************/

class TabPage : public QWidget
{
  Q_OBJECT

  Tree *_tree;
  QFrame *barFileChanged;
  QFrame *barUnknownElements;
  QFrame *barRenamedElements;
  QLabel *labelUnknownElements;
  QLabel *labelRenamedElements;

public:
  TabPage(Tree *tree, TabEnvironment *parent);
 ~TabPage() {}

  Tree *tree() { return _tree; }

private slots:
  void treeFileChanged(const QString &fileName, bool changed);
  void unknownElementsChanged(int count);
  void renamedElementsChanged(int count);
  void fileChangeClickedReopen();
  void barClickedIgnore();
};

/***********************************************/

#endif
