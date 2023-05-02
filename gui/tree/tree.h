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
#include <QPointer>
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
class TreeElementProgram;
class QMimeData;
class QSettings;
class VariableList;

/***** CLASS ***********************************/

class Tree: public QTreeWidget
{
   Q_OBJECT

  friend class TreeElementGlobal;

  QSettings             *settings;
  ActionList             actionList;
  Schema                 schema;
  TreeElement          *_rootElement;
  TreeElementGlobal    *_elementGlobal; // set by TreeElementGlobal
  QDir                   xmlDir;
  QString                xmlFile;
  QString                xsdFile;
  Bool                   changed;
  TreeItem              *selectedItem;
  int                    heightSelectedItem;
  VariableList          _varList;  // values of elements in the global list
  Bool                  _showResults;
  QUndoStack           *_undoStack;
  TabEnvironment        *workspace;
  QSet<TreeElement*>     unknownElements;
  QSet<TreeElement*>     renamedElements;
  QString               _programType; // programType or programmeType (for backwards compatibility)

  void clearTree();
  Bool readSchema();

public:
  Tree(QWidget *parent, ActionList *actionList, TabEnvironment *workspace);
  virtual ~Tree();

  TreeElement        *rootElement()   const {return _rootElement;}
  TreeElementGlobal  *elementGlobal() const {return _elementGlobal;}
  const VariableList &varList()       const {return _varList;}
  VariableList       &setVarList()          {return _varList;}
  QString            schemaFile()     const {return xsdFile;}
  QString            programType()    const {return _programType;}

  QUndoStack *undoStack() const                      {return _undoStack;}
  void        pushUndoCommand(QUndoCommand *command) {_undoStack->push(command);}

/**
* @brief the current fileName with absolute path.
*/
QString fileName() const;

/**
* @brief Sets the current TreeElement/TreeItem.
*
* Creates the editor elements.
* Updates the action list.
*/
void setSelectedItem(TreeItem *item);

/**
* @brief Gets the current TreeItem.
*/
TreeItem *getSelectedItem() const {return selectedItem;}

/**
* @brief show/hide column with the annotations depending on state.
*/
void setShowDescriptions(Bool state);

/**
* @brief show/hide results of links and expressions depending on state.
*/
void setShowResults(Bool state);

/**
* @brief Are results of links and expressions are shown?.
*/
Bool showResults() const {return _showResults;}

/**
* @brief Update expressions of all elements in the tree based on the @p element.
*/
void updateExpressions(const TreeElement *element);


/**
* @brief Get list of available programs from the schema.
*/
std::vector<XsdElementPtr> programListFromSchema() const;

/**
* @brief Is it allowed to close the file?
*
* The user is asked to save the file before if necessary.
* @return allowed
*/
Bool okToAbandon();

/**
* @brief Creates a new empty tree.
*
* The user is asked to save the old file before if necessary.
* @return success?
*/
Bool newFile();

/**
* @brief Reload the file.
*
* The user is asked to save the file before if necessary.
* @return success?
*/
Bool reopenFile();

/**
* @brief Open the file with @a fileName.
*
* If @a fileName is empty a file selector is opened.
* The user is asked to save the old file before if necessary.
* @return success?
*/
Bool openFile(QString fileName = QString());

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
* @brief Execute file.
*
* opens the execute dialog.
* @return success?
*/
Bool execFile();

/**
* @brief Has the current file unsaved changes?
*/
Bool isChanged() const;

signals:
/**
* @brief This signal is emitted when the file or status is changed.
*
* @param fileName the current file (empty if file is new)
* @param changed current file has unsaved changes
*/
void fileChanged(const QString &fileName, bool changed);

/**
* @brief This signal is emitted when a selection in a tree is changed.
*
* @param selection Name of the selected item or empty if no selection.
*/
void treeSelectionChanged(const QString &selection);

/**
* @brief This signal is emitted when a file is changed.
*
* @param fileName of the changed file
* @param changed file has unsaved changes
*/
void treeFileChanged(const QString &fileName, bool changed);

// check state of file on disk
private:
  QFileSystemWatcher *fileWatcher;

  void createFileWatcher();
  void clearFileWatcher();

private slots:
  void fileChangedExternally();

// tracking of unknown elements
// ----------------------------
public:
void trackUnknownElement(TreeElement *element);
bool untrackUnknownElement(TreeElement *element);

signals:
void unknownElementsChanged(int count);

public slots:
void expandUnknownElements();
void removeAllUnknownElements();

// tracking of renamed elements
// ----------------------------
public:
void trackRenamedElement(TreeElement *element);
bool untrackRenamedElement(TreeElement *element);

signals:
void renamedElementsChanged(int count);

public slots:
void expandRenamedElements();
void updateAllRenamedElements();

// management of file names
// ------------------------
public:
/**
* @brief makes the fileName absolute.
*/
QString addXmlDirectory(const QString &fileName) const;

/**
* @brief makes the fileName relative to the working directory if possible.
*/
QString stripXmlDirectory(const QString &fileName) const;

private:
// create mime data (XML) of the element
QMimeData *createMimeData(const TreeElement *element);

// creates xmlNode and type from mimeData if possible
Bool fromMimeData(const QMimeData *mimeData, XmlNodePtr &xmlNode, QString &type);

// enabling/disabling actions depending on the selected item
void updateActions();

public slots:
  void editCut();
  void editCopy();
  void editPaste();
  void editPasteOverwrite();
  void editAdd();
  void editRemove();
  void editSetGlobal();
  void editSetLoop();
  void editRemoveLoop();
  void editSetCondition();
  void editRemoveCondition();
  void editEnabled(bool checked);
  void editEnableAll();
  void editDisableAll();
  void editRename();
  void editUpdateName();
  void editComment();
  void editCollapseAll();
  void editOpenExternally();

// Event-Handler
// -------------
private slots:
  void treeCleanChanged(bool clean);
  void treeClipboardDataChanged();
  void treeContextMenuRequested(const QPoint &pos);
  void treeCurrentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem *previous);
  void treeItemClicked       (QTreeWidgetItem *item, int column);
  void treeItemDoubleClicked (QTreeWidgetItem *item, int column);
  void addProgram(int index);
  void resizeColumn(int logicalIndex, int oldSize, int newSize);

protected:
  void resizeEvent(QResizeEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

// Key events
// -----------
private:
  bool focusNextPrevChild(bool /*next*/) {return false;}
  void keyPressEvent(QKeyEvent *event);

// Drag & drop
// -----------
private:
  QPoint dragStartPosition;

  void mousePressEvent(QMouseEvent     *event);
  void mouseMoveEvent (QMouseEvent     *event);
  void dragEnterEvent (QDragEnterEvent *event);
  void dragMoveEvent  (QDragMoveEvent  *event);
  void dropEvent      (QDropEvent      *event);
};

/***********************************************/

#endif /* __GROOPSGUI__TREE__ */
