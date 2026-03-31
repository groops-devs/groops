/***********************************************/
/**
* @file tree.h
*
* @brief Editable tree from XML schema.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#ifndef __GROOPSGUI__TREE__
#define __GROOPSGUI__TREE__

#include <QSettings>
#include <QTreeWidget>
#include <QAction>
#include <QDir>
#include <QUndoStack>
#include <QFileSystemWatcher>
#include "base/importGroops.h"
#include "base/schema.h"
#include "mainWindow/mainWindow.h"
#include "mainWindow/tabs.h"

/***** TYPES ***********************************/

class TreeItem;
class TreeElement;
class TreeElementGlobal;
class TreeElementComplex;
class TreeWidget;

/***** CLASS ***********************************/

class Tree: public QWidget
{
  Q_OBJECT

  QSettings      settings;
  ActionList     actionList;
  bool           _showResults;
  bool           _isCurrent;
  TreeItem      *_selectedItem;
  Schema         _schema;
  XsdComplexPtr  _xsdGlobal;
  QString        _fileNameSchema;
  bool           _isClean;         // file is unchanged
  QString        _caption;         // fileName without path or something like 'newX'.
  QString        _fileName;        // absolute path
  QDir           workingDirectory;

  void clearTree();
  bool readSchema();

  friend class TreeItem;

public:
  Tree(QWidget *parent, ActionList *actionList, TabEnvironment *tabEnvironment);
  virtual ~Tree();

  // public variables
  QUndoStack         *undoStack;
  TreeElementComplex *rootElement;
  TreeElementGlobal  *elementGlobal; // set by TreeElementGlobal
  XsdElementPtr       xsdElementLoop;
  XsdElementPtr       xsdElementCondition;

  /** @brief Is the tree selected?
  * Tree reacts in actions (like saveAs).  */
  bool isCurrent() {return _isCurrent;}

  /** @brief Sets the tree selected. */
  void setCurrent(bool isCurrent);

  /** @brief Sets the current TreeElement/TreeItem.
  * Creates the editor elements.
  * Updates the action list. */
  void setSelectedItem(TreeItem *item);

  /** @brief Gets the current TreeItem. */
  TreeItem *selectedItem() const {return _selectedItem;}

  /** @brief Gets the current TreeElement. */
  TreeElement *selectedElement() const;

  /** @brief Are results of links and expressions are shown?. */
  bool showResults() const {return _showResults;}

  /** @brief fileName without path or something like 'newX'. */
  QString caption() const {return _caption;}

  /** @brief the current fileName with absolute path. */
  QString fileName() const {return _fileName;}

  QString fileNameSchema() const {return _fileNameSchema;}

  /** @brief makes the fileName absolute. */
  QString addWorkingDirectory(const QString &fileName) const;

  /** @brief makes the fileName relative to the working directory if possible. */
  QString stripWorkingDirectory(const QString &fileName) const;

  /** @brief Does the current file have no unsaved changes? */
  bool isClean() const {return _isClean;}

  /** @brief Is it allowed to close the file?
  * The user is asked to save the file before if necessary.
  * @return allowed */
  bool okToAbandon();

  /** @brief Creates a new empty tree.
  * The user is asked to save the old file before if necessary.
  * @return success? */
  bool fileNew(const QString &caption);

  /** @brief Open the file with @a fileName.
  * If @a fileName is empty try to reload the current file.
  * The user is asked to save the old file before if necessary.
  * @return success? */
  bool fileOpen(QString fileName);

signals:
  /** @brief This signal is emitted when the file or status is changed.
  * @param caption fileName without path or something like 'newX'.
  * @param fileName the current file (empty if file is new)
  * @param isClean no unsaved changes */
  void fileChanged(const QString &caption, const QString &fileName, bool isClean);

public:
  const std::vector<XsdElementPtr> &xsdElements() const;
  XsdElementPtr xsdElement(const QString &type) const;
  QList<XsdElementPtr> programList() const;

public slots:
  void addProgram(const QString &name);

  // check state of file on disk
  // ---------------------------
private:
  QFileSystemWatcher *fileWatcher;

  void fileWatcherCreate();
  void fileWatcherClear();

private slots:
  void fileWatcherChangedExternally();

private:
  // enabling/disabling actions depending on the selected item
  void updateActions();

public slots:
  bool fileSave();
  bool fileSaveAs();
  void fileRun();
  void fileShowInManager();
  void editCut();
  void editCopy();
  void editPaste();
  void editPasteOverwrite();
  void editAdd();
  void editRemove();
  void editAddVariable();
  void editSetLoop();
  void editSetCondition();
  void editEnabled(bool checked);
  void editEnableAll();
  void editDisableAll();
  void editRename();
  void editUpdateName();
  void editAddComment();
  void editCollapseAll();
  void editOpenExternally();
  void helpShowDescriptions(bool);
  void helpShowResults(bool);
  void helpOpenDocumentation();

  // ========================================================

  // Layout
  // ------
  int  columnWidth(int column) const;
  void setColumnWidth(int column, int width);

signals:
  void sectionResized(int logicalIndex, int oldSize, int newSize);

private:
  QFrame     *barFileExternallyChanged;
  QFrame     *barBrokenLinks;
  QFrame     *barUnknownElements;
  QFrame     *barSchemaRenamedElements;
  QFrame     *barDeprecatedElements;
  QLabel     *labelBrokenLinks;
  QLabel     *labelUnknownElements;
  QLabel     *labelSchemaRenamedElements;
  QLabel     *labelDeprecatedElements;
  int         brokenLinkCount, unknownCount, renamedCount, deprecatedCount;
  TreeWidget *treeWidget;

public:
  void treeChanged();   // slot: called by TreeElements whenever unknowns or renamed could be changed

private slots:
  void undoStackCleanChanged(bool clean);
  void treeClipboardDataChanged();
  void treeContextMenuRequested(const QPoint &pos);
  void treeCurrentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem *previous);
  void treeItemSelectionChanged();
  void treeItemClicked       (QTreeWidgetItem *item, int column);
  void treeItemDoubleClicked (QTreeWidgetItem *item, int column);
  void barFileExternallyChangedReopen();
  void barBrokenLinksExpand();
  void barUnknownElementsExpand();
  void barUnknownElementsRemoveAll();
  void barSchemaRenamedElementsExpand();
  void barSchemaRenamedElementsUpdateAll();
  void barDeprecatedElementsExpand();
  void barClickedIgnore();
};

/***********************************************/

class TreeWidget: public QTreeWidget
{
  Q_OBJECT

  Tree           *tree;
  TabEnvironment *tabEnvironment;
  QPoint          dragStartPosition;
  TreeElement    *dragElement;

public:
  TreeWidget(Tree *tree, TabEnvironment *tabEnvironment);
  virtual ~TreeWidget() {}

protected:
  bool eventFilter(QObject *obj, QEvent *event) override;
  bool focusNextPrevChild(bool /*next*/) override {return false;}
  void keyPressEvent  (QKeyEvent       *event) override;
  void mousePressEvent(QMouseEvent     *event) override;
  void mouseMoveEvent (QMouseEvent     *event) override;
  void dragEnterEvent (QDragEnterEvent *event) override;
  void dragMoveEvent  (QDragMoveEvent  *event) override;
  void dropEvent      (QDropEvent      *event) override;
};

/***********************************************/

#endif /* __GROOPSGUI__TREE__ */
