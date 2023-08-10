/***********************************************/
/**
* @file tree.cpp
*
* @brief Editable tree from XML schema.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#include <QtDebug>
#include <QMimeData>
#include <QDrag>
#include <QSettings>
#include <QTreeWidget>
#include <QHeaderView>
#include <QClipboard>
#include <QProcess>
#include <QApplication>
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QMenu>
#include <QComboBox>
#include <QContextMenuEvent>
#include <QTextStream>
#include <QUrl>
#include <QDesktopServices>
#include <QVBoxLayout>
#include "base/importGroops.h"
#include "base/xml.h"
#include "base/schema.h"
#include "tree/treeItem.h"
#include "tree/treeElement.h"
#include "tree/treeElementAdd.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementFileName.h"
#include "tree/treeElementGlobal.h"
#include "tree/treeElementProgram.h"
#include "tree/treeElementComment.h"
#include "tree/treeElementUnknown.h"
#include "executeDialog/executeDialog.h"
#include "setLoopConditionDialog/setLoopConditionDialog.h"
#include "settingsDialog/settingsPathDialog.h"
#include "tree/tree.h"
#ifdef _WIN32
#undef  NOMINMAX
#define NOMINMAX 1
#include "windows.h"
#endif

/***********************************************/
/***********************************************/

Tree::Tree(QWidget *parent, ActionList *actionList, TabEnvironment *tabEnvironment) : QWidget(parent)
{
  try
  {
    this->_isClean      = true;
    this->_selectedItem = nullptr;
    this->rootElement   = nullptr;
    this->elementGlobal = nullptr;
    this->fileWatcher   = nullptr;
    this->undoStack     = new QUndoStack(this);
    this->_isCurrent    = true;

    // Layout
    // ======
    // Bar handling external file changes
    // ----------------------------------
    {
      QPushButton *buttonReopen = new QPushButton(QIcon(":/icons/scalable/view-refresh.svg"), "Reopen", this);
      QPushButton *buttonIgnore = new QPushButton(QIcon(":/icons/scalable/ignore.svg"), "Ignore", this);
      QLabel *iconLabel = new QLabel(this);
      iconLabel->setPixmap(QIcon(":/icons/scalable/warning.svg").pixmap(24,24));
      QHBoxLayout *layoutBar = new QHBoxLayout(this);
      layoutBar->addWidget(iconLabel);
      layoutBar->addWidget(new QLabel("File was modified externally. Reopen?", this), 1);
      layoutBar->addWidget(buttonReopen);
      layoutBar->addWidget(buttonIgnore);
      layoutBar->setContentsMargins(3, 3, 3, 3);
      barFileExternallyChanged = new QFrame(this);
      barFileExternallyChanged->setFrameStyle(QFrame::Box);
      barFileExternallyChanged->setLayout(layoutBar);
      barFileExternallyChanged->setVisible(false);
      const QString highlightColor = barFileExternallyChanged->palette().highlight().color().name().right(6);
      barFileExternallyChanged->setStyleSheet(".QFrame { color: #"+highlightColor+"; background-color: #4d"+highlightColor+" }");

      connect(buttonReopen, SIGNAL(clicked()), this, SLOT(barFileExternallyChangedReopen()));
      connect(buttonIgnore, SIGNAL(clicked()), this, SLOT(barClickedIgnore()));
    }

    // Bar handling unknown elements
    // -----------------------------
    {
      QPushButton *buttonShowAll   = new QPushButton(QIcon(":/icons/scalable/edit-find-replace.svg"), "Show all", this);
      QPushButton *buttonRemoveAll = new QPushButton(QIcon(":/icons/scalable/edit-delete.svg"), "Remove all", this);
      QPushButton *buttonIgnore    = new QPushButton(QIcon(":/icons/scalable/ignore.svg"), "Ignore", this);
      buttonRemoveAll->setMinimumWidth(95);
      QLabel *iconLabel = new QLabel(this);
      iconLabel->setPixmap(QIcon(":/icons/scalable/help-about.svg").pixmap(24,24));
      QHBoxLayout *layoutBar = new QHBoxLayout(this);
      labelUnknownElements = new QLabel(this);
      layoutBar->addWidget(iconLabel);
      layoutBar->addWidget(labelUnknownElements, 1);
      layoutBar->addWidget(buttonShowAll);
      layoutBar->addWidget(buttonRemoveAll);
      layoutBar->addWidget(buttonIgnore);
      layoutBar->setContentsMargins(3, 3, 3, 3);
      barUnknownElements = new QFrame(this);
      barUnknownElements->setFrameStyle(QFrame::Box);
      barUnknownElements->setLayout(layoutBar);
      barUnknownElements->setVisible(false);
      const QString highlightColor = barUnknownElements->palette().highlight().color().name().right(6);
      barUnknownElements->setStyleSheet(".QFrame { color: #"+highlightColor+"; background-color: #4d"+highlightColor+" }");

      connect(buttonShowAll,   SIGNAL(clicked()), this, SLOT(barUnknownElementsExpand()));
      connect(buttonRemoveAll, SIGNAL(clicked()), this, SLOT(barUnknownElementsRemoveAll()));
      connect(buttonIgnore,    SIGNAL(clicked()), this, SLOT(barClickedIgnore()));
    }

    // Bar handling schema renamed elements
    // ------------------------------------
    {
      QPushButton *buttonShowAll   = new QPushButton(QIcon(":/icons/scalable/edit-find-replace.svg"), "Show all", this);
      QPushButton *buttonUpdateAll = new QPushButton(QIcon(":/icons/scalable/edit-rename.svg"), "Update all", this);
      QPushButton *buttonIgnore    = new QPushButton(QIcon(":/icons/scalable/ignore.svg"), "Ignore", this);
      buttonUpdateAll->setMinimumWidth(95);
      QLabel *iconLabel = new QLabel(this);
      iconLabel->setPixmap(QIcon(":/icons/scalable/help-about.svg").pixmap(24,24));
      QHBoxLayout *layoutBar = new QHBoxLayout(this);
      labelSchemaRenamedElements = new QLabel(this);
      layoutBar->addWidget(iconLabel);
      layoutBar->addWidget(labelSchemaRenamedElements, 1);
      layoutBar->addWidget(buttonShowAll);
      layoutBar->addWidget(buttonUpdateAll);
      layoutBar->addWidget(buttonIgnore);
      layoutBar->setContentsMargins(3, 3, 3, 3);
      barSchemaRenamedElements = new QFrame(this);
      barSchemaRenamedElements->setFrameStyle(QFrame::Box);
      barSchemaRenamedElements->setLayout(layoutBar);
      barSchemaRenamedElements->setVisible(false);
      const QString highlightColor = barSchemaRenamedElements->palette().highlight().color().name().right(6);
      barSchemaRenamedElements->setStyleSheet(".QFrame { color: #"+highlightColor+"; background-color: #4d"+highlightColor+" }");

      connect(buttonShowAll,   SIGNAL(clicked()), this, SLOT(barSchemaRenamedElementsExpand()));
      connect(buttonUpdateAll, SIGNAL(clicked()), this, SLOT(barSchemaRenamedElementsUpdateAll()));
      connect(buttonIgnore,    SIGNAL(clicked()), this, SLOT(barClickedIgnore()));
    }

    this->treeWidget = new TreeWidget(this, tabEnvironment);

    QVBoxLayout *layout = new QVBoxLayout();
    layout->setSpacing(3);
    layout->setContentsMargins(0, 3, 0, 0);
    layout->addWidget(barFileExternallyChanged);
    layout->addWidget(barUnknownElements);
    layout->addWidget(barSchemaRenamedElements);
    layout->addWidget(treeWidget, 1);
    setLayout(layout);

    // QTreeWidget events
    // ------------------
    connect(QApplication::clipboard(),              SIGNAL(dataChanged()),                                          this, SLOT(treeClipboardDataChanged()));
    connect(qobject_cast<QTreeWidget*>(treeWidget), SIGNAL(customContextMenuRequested(const QPoint &)),             this, SLOT(treeContextMenuRequested (const QPoint&)));
    connect(qobject_cast<QTreeWidget*>(treeWidget), SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)), this, SLOT(treeCurrentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)));
    connect(qobject_cast<QTreeWidget*>(treeWidget), SIGNAL(itemSelectionChanged()),                                 this, SLOT(treeItemSelectionChanged()));
    connect(qobject_cast<QTreeWidget*>(treeWidget), SIGNAL(itemClicked       (QTreeWidgetItem*, int)),              this, SLOT(treeItemClicked       (QTreeWidgetItem*, int)));
    connect(qobject_cast<QTreeWidget*>(treeWidget), SIGNAL(itemDoubleClicked (QTreeWidgetItem*, int)),              this, SLOT(treeItemDoubleClicked (QTreeWidgetItem*, int)));
    connect(treeWidget->header(), SIGNAL(sectionResized(int, int, int)), this, SIGNAL(sectionResized(int, int, int)));

    // connection to actions
    // ---------------------
    this->actionList = *actionList;
    connect(this->actionList.fileSaveAction,              SIGNAL(triggered(bool)), this, SLOT(fileSave()));
    connect(this->actionList.fileSaveAsAction,            SIGNAL(triggered(bool)), this, SLOT(fileSaveAs()));
    connect(this->actionList.fileRunAction,               SIGNAL(triggered(bool)), this, SLOT(fileRun()));
    connect(this->actionList.fileShowInManagerAction,     SIGNAL(triggered(bool)), this, SLOT(fileShowInManager()));
    connect(this->actionList.editCutAction,               SIGNAL(triggered(bool)), this, SLOT(editCut()));
    connect(this->actionList.editCopyAction,              SIGNAL(triggered(bool)), this, SLOT(editCopy()));
    connect(this->actionList.editPasteAction,             SIGNAL(triggered(bool)), this, SLOT(editPaste()));
    connect(this->actionList.editPasteOverwriteAction,    SIGNAL(triggered(bool)), this, SLOT(editPasteOverwrite()));
    connect(this->actionList.editAddAction,               SIGNAL(triggered(bool)), this, SLOT(editAdd()));
    connect(this->actionList.editRemoveAction,            SIGNAL(triggered(bool)), this, SLOT(editRemove()));
    connect(this->actionList.editSetGlobalAction,         SIGNAL(triggered(bool)), this, SLOT(editSetGlobal()));
    connect(this->actionList.editSetLoopAction,           SIGNAL(triggered(bool)), this, SLOT(editSetLoop()));
    connect(this->actionList.editRemoveLoopAction,        SIGNAL(triggered(bool)), this, SLOT(editRemoveLoop()));
    connect(this->actionList.editSetConditionAction,      SIGNAL(triggered(bool)), this, SLOT(editSetCondition()));
    connect(this->actionList.editRemoveConditionAction,   SIGNAL(triggered(bool)), this, SLOT(editRemoveCondition()));
    connect(this->actionList.editEnabledAction,           SIGNAL(triggered(bool)), this, SLOT(editEnabled(bool)));
    connect(this->actionList.editEnableAllAction,         SIGNAL(triggered(bool)), this, SLOT(editEnableAll()));
    connect(this->actionList.editDisableAllAction,        SIGNAL(triggered(bool)), this, SLOT(editDisableAll()));
    connect(this->actionList.editRenameAction,            SIGNAL(triggered(bool)), this, SLOT(editRename()));
    connect(this->actionList.editUpdateNameAction,        SIGNAL(triggered(bool)), this, SLOT(editUpdateName()));
    connect(this->actionList.editAddCommentAction,        SIGNAL(triggered(bool)), this, SLOT(editAddComment()));
    connect(this->actionList.editCollapseAllAction,       SIGNAL(triggered(bool)), this, SLOT(editCollapseAll()));
    connect(this->actionList.editOpenExternallyAction,    SIGNAL(triggered(bool)), this, SLOT(editOpenExternally()));
    connect(this->actionList.helpShowDescriptionsAction,  SIGNAL(triggered(bool)), this, SLOT(helpShowDescriptions(bool)));
    connect(this->actionList.helpShowResultsAction,       SIGNAL(triggered(bool)), this, SLOT(helpShowResults(bool)));
    connect(this->actionList.helpOpenDocumentationAction, SIGNAL(triggered(bool)), this, SLOT(helpOpenDocumentation()));

    actionList->helpShowDescriptionsAction->setChecked(settings.value("misc/showDescriptions", true).toBool());
    actionList->helpShowResultsAction->setChecked(settings.value("misc/showResults", true).toBool());
    helpShowDescriptions(actionList->helpShowDescriptionsAction->isChecked());
    helpShowResults(actionList->helpShowResultsAction->isChecked());

    // track state of file via undoStack
    connect(undoStack,  SIGNAL(cleanChanged(bool)), this, SLOT(undoStackCleanChanged(bool)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Tree::~Tree()
{
  const QSignalBlocker blocker(this);
  clearTree();
}

/***********************************************/

void Tree::clearTree()
{
  try
  {
    fileWatcherClear();
    setSelectedItem(nullptr);
    _isClean = true;
    undoStack->clear();
    unknownCount = renamedCount = 0;
    barFileExternallyChanged->setHidden(true);
    barUnknownElements->setHidden(true);
    barSchemaRenamedElements->setHidden(true);

    if(rootElement)
      rootElement->removeItem();
    delete rootElement;
    rootElement   = nullptr;
    elementGlobal = nullptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool Tree::readSchema()
{
  try
  {
    QString fileNameSchema = settings.value("files/schemaFile").toString();
    Schema  schema;
    if(!schema.readFile(fileNameSchema))
    {
      QMessageBox::critical(this , tr("GROOPS"), tr("File '%1' seems not to be a valid XSD schema").arg(fileNameSchema));
      return false;
    }
    _fileNameSchema = fileNameSchema;
    _schema         = schema;
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void Tree::setCurrent(bool isCurrent)
{
  if(!isCurrent)
    setSelectedItem(nullptr);
  _isCurrent = isCurrent;
}

/***********************************************/

void Tree::setSelectedItem(TreeItem *item)
{
  try
  {
    if(item == _selectedItem)
      return;

    // old item lost selection
    if(_selectedItem)
      _selectedItem->lostCurrent();

    // new item gets selection
    _selectedItem = item;
    if(_selectedItem)
    {
      _selectedItem->becomeCurrent();
      updateActions();
    }

    if(_selectedItem != treeWidget->currentItem())
      treeWidget->setCurrentItem(item);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElement *Tree::selectedElement() const
{
  return _selectedItem ? _selectedItem->treeElement() : nullptr;
}

/***********************************************/
/***********************************************/

QString Tree::addWorkingDirectory(const QString &filename) const
{
  try
  {
    QFileInfo file(filename);
    QString s = file.filePath();
    if(file.isRelative())
      file.setFile(workingDirectory, s);
    return file.absoluteFilePath();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QString Tree::stripWorkingDirectory(const QString &filename) const
{
  try
  {
    QString path = workingDirectory.absolutePath();
    if(path.isEmpty())
      return filename;
    QFileInfo file(filename);
    if(file.isRelative())
      return filename;
    QString s = file.absoluteFilePath();
    path += QString("/");
    if(s.startsWith(path))
      s.remove(0, path.length());

    return QDir::toNativeSeparators(s);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool Tree::okToAbandon()
{
  try
  {
    if(isClean())
      return true;

    auto button = QMessageBox::question(this , tr("Close File - GROOPS"), tr("File '%1' was changed.\nDo you want to save the file?").arg(caption()),
                                        QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel, QMessageBox::Save);
    if(button == QMessageBox::Save)
      return fileSave();
    return (button == QMessageBox::Discard);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool Tree::fileNew(const QString &caption)
{
  try
  {
    if(!okToAbandon())
      return false;
    if(!readSchema())
      return false;

    // try to open template file, clear tree if not
    blockSignals(true); // do not emit fileChanged
    QString templateFile = settings.value("files/templateFile").toString();
    if(!QFileInfo(templateFile).isFile() || !fileOpen(templateFile))
    {
      clearTree();
      rootElement = dynamic_cast<TreeElementComplex*>(TreeElement::newTreeElement(this, nullptr, _schema.rootElement, "", XmlNodePtr(nullptr), true/*fillWithDefaults*/));
      TreeItem *item = rootElement->createItem(nullptr, nullptr);
      if(item)
        item->setExpanded(true);
    }
    blockSignals(false);

    _caption  = caption;
    _fileName = QString();
    workingDirectory.setPath(settings.value("files/workingDirectory").toString());
    emit fileChanged(caption, fileName(), isClean());

    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool Tree::fileOpen(QString fileName)
{
  try
  {
    if(!okToAbandon())
      return false;

    if(fileName.isEmpty())
      fileName = _fileName; // try to reopen
    if(fileName.isEmpty())
      fileNew(caption());

    if(!readSchema())
      return false;

    // read xml file
    // -------------
    QFile file(fileName);
    if(!file.open(QFile::ReadOnly | QFile::Text))
    {
      QMessageBox::critical(this , tr("Open file - GROOPS"), tr("Cannot open file '%1'.").arg(fileName));
      return false;
    }
    QString errorMessage;
    XmlNodePtr xmlNode = XmlNode::read(&file, errorMessage);
    if(!xmlNode)
    {
      QMessageBox::critical(this , tr("Open file - GROOPS"), errorMessage);
      return false;
    }

    // seems to be all ok
    // ------------------
    _caption          = QFileInfo(fileName).fileName();
    _fileName         = fileName;
    workingDirectory  = QFileInfo(fileName).dir();
    emit fileChanged(caption(), fileName, isClean());
    fileWatcherCreate();

    clearTree();
    TreeElement *element = TreeElement::newTreeElement(this, nullptr, _schema.rootElement, "", xmlNode, false/*fillWithDefaults*/);
    rootElement = dynamic_cast<TreeElementComplex*>(element); // set rootElement after complete initialization
    if(elementGlobal)
    {
      elementGlobal->informAboutGlobalElements(rootElement, true/*recursively*/);
      elementGlobal->updateVariableList();
    }
    rootElement->createItem(nullptr, nullptr)->setExpanded(true);
    treeChanged();
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QList<XsdElementPtr> Tree::programList() const
{
  return _schema.programList();
}

/***********************************************/

void Tree::addProgram(const QString &name)
{
  try
  {
    undoStack->beginMacro("add program "+name);
    TreeElement *elementProgram = rootElement->addChild(rootElement->children().back(), "programType", XmlNodePtr(nullptr));
    elementProgram->changeSelectedIndex(elementProgram->findValueIndex(name));
    undoStack->endMacro();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void Tree::fileWatcherCreate()
{
  try
  {
    fileWatcherClear();
    if(fileName().isEmpty())
      return;
    fileWatcher = new QFileSystemWatcher(this);
    fileWatcher->addPath(fileName());
    connect(fileWatcher, &QFileSystemWatcher::fileChanged, this, &Tree::fileWatcherChangedExternally);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::fileWatcherClear()
{
  if(fileWatcher)
  {
    disconnect(fileWatcher, &QFileSystemWatcher::fileChanged, this, &Tree::fileWatcherChangedExternally);
    delete fileWatcher;
    fileWatcher = nullptr;
  }
}

/***********************************************/

void Tree::fileWatcherChangedExternally()
{
  barFileExternallyChanged->setVisible(true);
}

/***********************************************/
/***********************************************/

// mime data <-> XML of the element
// --------------------------------
static const char *mimeFormatString = "application/x-groops"; //"text/plain";

static QMimeData *createMimeData(const TreeElement *element)
{
  try
  {
    if(element->type().isEmpty())
      return nullptr;
    XmlNodePtr xmlNode = element->createXmlTree(true/*createRootEvenIfEmpty*/);
    if(!xmlNode)
      return nullptr;
    writeAttribute(xmlNode, "xsdType", element->type());

    // create xml text
    QString xmlData;
    QTextStream stream(&xmlData, QIODevice::WriteOnly);
    XmlNode::write(stream, xmlNode, true/*writeCommentsAsElements*/);

    // create mime data
    QMimeData *mimeData = new QMimeData;
    mimeData->setData(mimeFormatString, xmlData.toUtf8());
    return mimeData;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

static bool fromMimeData(const QMimeData *mimeData, XmlNodePtr &xmlNode, QString &type)
{
  try
  {
    if(!mimeData->hasFormat(mimeFormatString))
      return false;
    // create xmlData from MIME data
    QString errorMessage;
    xmlNode = XmlNode::read(mimeData->data(mimeFormatString), errorMessage);
    if(!xmlNode)
      return false;
    readAttribute(xmlNode, "xsdType", type);
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void Tree::updateActions()
{
  try
  {
    if(!selectedElement())
      return;
    TreeElement        *element       = selectedElement();
    TreeElementComplex *parentElement = element->parentElement;

    // test content of clipboard
    QString    type;
    XmlNodePtr xmlNode;
    bool isClipboard = fromMimeData(QApplication::clipboard()->mimeData(), xmlNode, type);
    bool canCopy = !element->type().isEmpty() && element->createXmlTree(true);

    actionList.fileShowInManagerAction    ->setEnabled(!fileName().isEmpty() );
    actionList.editCutAction              ->setEnabled(canCopy && parentElement && parentElement->canRemoveChild(element));
    actionList.editCopyAction             ->setEnabled(canCopy);
    actionList.editPasteAction            ->setEnabled(isClipboard && ((parentElement && parentElement->canAddChild(element, type)) ||
                                                       (dynamic_cast<TreeElementGlobal*>(element) && dynamic_cast<TreeElementGlobal*>(element)->canAddChild(element, type))));
    actionList.editPasteOverwriteAction   ->setEnabled(isClipboard && element->canOverwrite(type));
    actionList.editAddAction              ->setEnabled((parentElement && parentElement->canAddChild(element, element->type())) || dynamic_cast<TreeElementGlobal*>(element));
    actionList.editRemoveAction           ->setEnabled(parentElement && parentElement->canRemoveChild(element));
    actionList.editSetGlobalAction        ->setEnabled(elementGlobal && elementGlobal->canSetGlobal(element));
    actionList.editSetLoopAction          ->setEnabled(element->canSetLoop());
    actionList.editRemoveLoopAction       ->setEnabled(!element->loop().isEmpty());
    actionList.editSetConditionAction     ->setEnabled(element->canSetCondition());
    actionList.editRemoveConditionAction  ->setEnabled(!element->condition().isEmpty());
    actionList.editEnabledAction          ->setEnabled(element->canDisabled());
    actionList.editEnabledAction          ->setChecked(!element->disabled());
    actionList.editEnableAllAction        ->setEnabled(true);
    actionList.editDisableAllAction       ->setEnabled(true);
    actionList.editRenameAction           ->setEnabled(dynamic_cast<TreeElementGlobal*>(parentElement) && dynamic_cast<TreeElementGlobal*>(parentElement)->canRenameChild(element));
    actionList.editUpdateNameAction       ->setEnabled(element->canUpdateName());
    actionList.editAddCommentAction       ->setEnabled(parentElement && parentElement->canAddChild(element, "COMMENT"));
    actionList.editCollapseAllAction      ->setEnabled(true);
    TreeElementFileName *fileNameElement = dynamic_cast<TreeElementFileName*>(element);
    actionList.editOpenExternallyAction   ->setEnabled(fileNameElement && QFileInfo(addWorkingDirectory(fileNameElement->selectedResult())).isFile());
    actionList.helpOpenDocumentationAction->setEnabled(true);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/**** slots ************************************/
/***********************************************/

bool Tree::fileSave()
{
  try
  {
    if(!isCurrent())
      return false;

    if(_fileName.isEmpty())
      return fileSaveAs();

    fileWatcherClear();
    XmlNodePtr xmlNode = rootElement->createXmlTree();
    XmlNode::writeFile(fileName(), xmlNode);
    undoStack->setClean();
    fileWatcherCreate();
    emit fileChanged(caption(), fileName(), isClean());
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool Tree::fileSaveAs()
{
  try
  {
    if(!isCurrent())
      return false;

    // Open file selector
    QString name = _fileName;
    if(name.isEmpty())
      name = settings.value("files/workingDirectory").toString(); // selecttor starts in default working directory
    name = QFileDialog::getSaveFileName(this, tr("Save file - GROOPS"), name, tr("XML files (*.xml)"));
    if(name.isEmpty())
      return false;
    if(!name.endsWith(".xml"))
      name += ".xml";

    _fileName         = name;
    _caption          = QFileInfo(name).fileName();
    workingDirectory  = QFileInfo(name).dir();

    return fileSave();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::fileRun()
{
  try
  {
    if(!isCurrent())
      return;

    // file must be saved first
    // ------------------------
    if(!isClean() || fileName().isEmpty())
    {
      QMessageBox::StandardButton button =
      QMessageBox::warning(this , tr("Run file - GROOPS"),
                           tr("File '%1' was changed.\nYou have to save it first.").arg(caption()),
                           QMessageBox::Save | QMessageBox::Cancel, QMessageBox::Save);
      if(button != QMessageBox::Save)
        return;
      if(!fileSave())
        return;
    }

    // execute dialog
    // --------------
    ExecuteDialog dialog(this, this);
    if(!dialog.exec())
      return;

    // create execute command
    // ----------------------
    QStringList commandList  = settings.value("execute/commands").toStringList();
    int         commandIndex = settings.value("execute/commandIndex", int(0)).toInt();
    if((commandIndex<0)||(commandIndex>=commandList.size()))
      return;
    QString command = commandList[commandIndex];
    if(command.isEmpty())
      return;

    QString optionString;
    // append log file option?
    if(settings.value("execute/useLogFile", bool(false)).toBool())
      optionString += " -l "+settings.value("execute/logFile", QString("groops.log")).toString()+" ";
    optionString += fileName();

    command.replace("%f", optionString);
    command.replace("%w", workingDirectory.absolutePath());

    qWarning()<<"run command:"<<command;

    // execute command
    // ---------------
    #ifdef _WIN32
    // https://stackoverflow.com/questions/42051405/qprocess-with-cmd-command-does-not-result-in-command-line-window
    class DetachableProcess : public QProcess
    {
    public:
      DetachableProcess(QObject *parent=0) : QProcess(parent) {}
      void detach()
      {
        waitForStarted();
        setProcessState(QProcess::NotRunning);
      }
    };

    DetachableProcess process;
    process.setCreateProcessArgumentsModifier([](QProcess::CreateProcessArguments *args)
    {args->flags |= CREATE_NEW_CONSOLE;
      args->startupInfo->dwFlags &=~ STARTF_USESTDHANDLES;});
    process.start(QString("cmd.exe"), QStringList({"/k", command}));
    process.detach();
    #else
    QProcess::startDetached(command);
    #endif

    return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::fileShowInManager()
{
  if(isCurrent() && !fileName().isEmpty())
    QDesktopServices::openUrl(QUrl::fromLocalFile(QFileInfo(fileName()).absoluteDir().path()));
}

/***********************************************/

void Tree::editCut()
{
  try
  {
    if(!selectedElement() || !selectedElement()->parentElement || !selectedElement()->parentElement->canRemoveChild(selectedElement()))
      return;
    QMimeData *mimeData = createMimeData(selectedElement());
    if(mimeData)
    {
      QApplication::clipboard()->setMimeData(mimeData);
      selectedElement()->parentElement->removeChild(selectedElement());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editCopy()
{
  try
  {
    if(!selectedElement())
      return;
    // copy to clipboard
    QMimeData *mimeData = createMimeData(selectedElement());
    if(mimeData)
      QApplication::clipboard()->setMimeData(mimeData);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editPaste()
{
  try
  {
    if(!selectedElement())
      return;
    // paste from clipboard
    QString    type;
    XmlNodePtr xmlNode;
    if(!fromMimeData(QApplication::clipboard()->mimeData(), xmlNode, type))
      return;
    if(dynamic_cast<TreeElementGlobal*>(selectedElement()))
      dynamic_cast<TreeElementGlobal*>(selectedElement())->addChild(nullptr, type, xmlNode);
    else if(selectedElement()->parentElement)
      selectedElement()->parentElement->addChild(selectedElement(), type, xmlNode);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editPasteOverwrite()
{
  try
  {
    if(!selectedElement())
      return;
    // paste from clipboard
    QString    type;
    XmlNodePtr xmlNode;
    if(!fromMimeData(QApplication::clipboard()->mimeData(), xmlNode, type))
      return;
    selectedElement()->overwrite(type, xmlNode);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editAdd()
{
  try
  {
    TreeElement *element = selectedElement();
    if(!element || !element->parentElement)
      return;
    if(dynamic_cast<TreeElementGlobal*>(element))
      dynamic_cast<TreeElementGlobal*>(element)->addNewChild();
    else
      element->parentElement->addChild(element, element->type(), element->createXmlTree(true/*createRootEvenIfEmpty*/));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editRemove()
{
  try
  {
    if(!selectedElement() || !selectedElement()->parentElement || !selectedElement()->parentElement->canRemoveChild(selectedElement()))
      return;

    if(QMessageBox::question(this , tr("Remove element - GROOPS"), tr("Do you really want to remove this element?"),
                             QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Ok) == QMessageBox::Ok)
      selectedElement()->parentElement->removeChild(selectedElement());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editSetGlobal()
{
  try
  {
    if(selectedElement() && elementGlobal && elementGlobal->canSetGlobal(selectedElement()))
      elementGlobal->setGlobal(selectedElement());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editSetLoop()
{
  try
  {
    TreeElement *element = selectedElement();
    if(!element || !elementGlobal || !element->canSetLoop())
      return;
    SetLoopConditionDialog dialog(elementGlobal, "loop", this);
    if(dialog.exec())
    {
      QString loopName = dialog.name();
      if(!elementGlobal->names().contains(loopName, Qt::CaseInsensitive))
        elementGlobal->addNewChild("loopType", loopName);
      element->setLoop(loopName);
      updateActions();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editRemoveLoop()
{
  try
  {
    if(!selectedElement() || selectedElement()->loop().isEmpty())
      return;
    selectedElement()->setLoop("");
    updateActions();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editSetCondition()
{
  try
  {
    TreeElement *element = selectedElement();
    if(!element || !elementGlobal || !element->canSetCondition())
      return;
    SetLoopConditionDialog dialog(elementGlobal, "condition", this);
    if(dialog.exec())
    {
      QString conditionName = dialog.name();
      if(!elementGlobal->names().contains(conditionName, Qt::CaseInsensitive))
        elementGlobal->addNewChild("conditionType", conditionName);
      element->setCondition(conditionName);
      updateActions();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editRemoveCondition()
{
  try
  {
    if(!selectedElement() || selectedElement()->condition().isEmpty())
      return;
    selectedElement()->setCondition("");
    updateActions();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editEnabled(bool checked)
{
  try
  {
    if(selectedElement())
      selectedElement()->setDisabled(!checked);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editEnableAll()
{
  try
  {
    if(!isCurrent())
      return;

    undoStack->beginMacro("enable all programs");
    for(auto &element : rootElement->children())
      if(element->disabled() && dynamic_cast<TreeElementProgram*>(element))
        element->setDisabled(false);
    undoStack->endMacro();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editDisableAll()
{
  try
  {
    if(!isCurrent())
      return;

    undoStack->beginMacro("disable all programs");
    for(auto &element : rootElement->children())
      if(!element->disabled() && dynamic_cast<TreeElementProgram*>(element))
        element->setDisabled(true);
    undoStack->endMacro();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editRename()
{
  try
  {
    if(!selectedElement() || !dynamic_cast<TreeElementGlobal*>(selectedElement()->parentElement) ||
       !dynamic_cast<TreeElementGlobal*>(selectedElement()->parentElement)->canRenameChild(selectedElement()))
      return;
    dynamic_cast<TreeElementGlobal*>(selectedElement()->parentElement)->renameChild(selectedElement(), QString());
    updateActions();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editUpdateName()
{
  try
  {
    if(!selectedElement() || !selectedElement()->canUpdateName())
      return;
    selectedElement()->updateName();
    updateActions();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editAddComment()
{
  try
  {
    TreeElement *element = selectedElement();
    if(!element || !element->parentElement)
      return;
    element->parentElement->addChild(element, "COMMENT", nullptr);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editCollapseAll()
{
  try
  {
    if(!isCurrent())
      return;
    treeWidget->collapseAll();
    treeWidget->expand(treeWidget->model()->index(0,0));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editOpenExternally()
{
  try
  {
    if(selectedElement())
      QDesktopServices::openUrl(QUrl::fromLocalFile(addWorkingDirectory(selectedElement()->selectedResult())));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::helpShowDescriptions(bool state)
{
  settings.setValue("misc/showDescriptions", state);
  treeWidget->setColumnHidden(2, !state);
  if(!treeWidget->isColumnHidden(2))
    treeWidget->setColumnWidth(2, std::max(treeWidget->columnWidth(2), 100));
}

/***********************************************/

void Tree::helpShowResults(bool state)
{
  settings.setValue("misc/showResults", state);
  _showResults = state;
  if(rootElement)
    rootElement->createItem(nullptr, nullptr);
}

/***********************************************/

void Tree::helpOpenDocumentation()
{
  try
  {
    if(!isCurrent())
      return;

    QString path     = settings.value("files/documentationDirectory").toString()+"/";
    QString fileName = "index.html";
    TreeElement *element = selectedElement();
    if(element)
    {
      fileName = (dynamic_cast<TreeElementProgram*>(element) ? element->selectedValue() : element->type())+".html";
      if(!QFileInfo::exists(path+fileName))
        fileName = "index.html";
    }
    if(QFileInfo::exists(path+fileName))
      QDesktopServices::openUrl(QUrl::fromLocalFile(path+fileName));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

int Tree::columnWidth(int column) const
{
  return treeWidget->columnWidth(column);
}

/***********************************************/

void Tree::setColumnWidth(int column, int width)
{
  treeWidget->setColumnWidth(column, width);
}

/***********************************************/
/**** Event-Handler ****************************/
/***********************************************/

static void countRenamesAndUnknowns(const TreeElement *element, int &unknownCount, int &renamedCount)
{
  try
  {
    if(!element)
      return;
    if(dynamic_cast<const TreeElementComplex*>(element))
      for(const auto &child : dynamic_cast<const TreeElementComplex*>(element)->children())
        countRenamesAndUnknowns(child, unknownCount, renamedCount);
    if(dynamic_cast<const TreeElementUnknown*>(element) || element->isSelectionUnknown(element->selectedIndex()))
      unknownCount++;
    if(element->isRenamedInSchema() || element->isSelectionRenamedInSchema(element->selectedIndex()))
      renamedCount++;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::treeChanged()
{
  try
  {
    if(!rootElement)
      return;
    int unknownCountOld = unknownCount;
    int renamedCountOld = renamedCount;
    unknownCount = renamedCount = 0;
    countRenamesAndUnknowns(rootElement, unknownCount, renamedCount);
    if(!unknownCount || (unknownCount > unknownCountOld))
      barUnknownElements->setVisible(unknownCount);
    if(!renamedCount || (renamedCount > renamedCountOld))
      barSchemaRenamedElements->setVisible(renamedCount);
    labelUnknownElements->setText(QString("File contains %1 unknown elements.").arg(unknownCount));
    labelSchemaRenamedElements->setText(QString("File contains %1 elements that were renamed in the schema.").arg(renamedCount));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::undoStackCleanChanged(bool clean)
{
  try
  {
    _isClean = clean;
    updateActions();
    emit fileChanged(caption(), fileName(), isClean());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void Tree::treeClipboardDataChanged()
{
  updateActions();
}

/***********************************************/

void Tree::treeContextMenuRequested(const QPoint &pos)
{
  try
  {
    TreeItem *item = dynamic_cast<TreeItem*>(treeWidget->itemAt(pos));
    if(!item)
      return;
    setSelectedItem(item);

    QMenu *contextMenu = new QMenu(this);
    contextMenu->addAction(actionList.helpOpenDocumentationAction);
    contextMenu->addSeparator();
    contextMenu->addAction(actionList.editCutAction);
    contextMenu->addAction(actionList.editCopyAction);
    contextMenu->addAction(actionList.editPasteAction);
    contextMenu->addAction(actionList.editPasteOverwriteAction);
    contextMenu->addSeparator();
    contextMenu->addAction(actionList.editAddAction);
    contextMenu->addAction(actionList.editRemoveAction);
    contextMenu->addSeparator();
    contextMenu->addAction(actionList.editSetGlobalAction);
    contextMenu->addSeparator();
    contextMenu->addAction(actionList.editSetLoopAction);
    contextMenu->addAction(actionList.editRemoveLoopAction);
    contextMenu->addAction(actionList.editSetConditionAction);
    contextMenu->addAction(actionList.editRemoveConditionAction);
    contextMenu->addSeparator();
    contextMenu->addAction(actionList.editEnabledAction);
    contextMenu->addAction(actionList.editRenameAction);
    contextMenu->addAction(actionList.editAddCommentAction);
    contextMenu->addAction(actionList.editOpenExternallyAction);
    contextMenu->exec(QWidget::mapToGlobal(pos));
    delete contextMenu;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::treeCurrentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem */*previous*/)
{
  try
  {
    setSelectedItem(dynamic_cast<TreeItem*>(current));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}


/***********************************************/

void Tree::treeItemSelectionChanged()
{
  // catch edge case where clicked item changes because editor is removed from item above and therefore vertical position changes
  const QSignalBlocker blocker(treeWidget);
  for(QTreeWidgetItem *item : treeWidget->selectedItems())
    item->setSelected(item == selectedItem());
  if(selectedItem())
    selectedItem()->setSelected(true);
}

/***********************************************/

void Tree::treeItemClicked(QTreeWidgetItem *item, int column)
{
  try
  {
    if(column != 1) // value column?
      return;
    if(item != _selectedItem)
      setSelectedItem(dynamic_cast<TreeItem*>(item));
    _selectedItem->setFocus();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::treeItemDoubleClicked(QTreeWidgetItem *item, int column)
{
  try
  {
    if(column != 3) // comment column?
      return;
    if(item != _selectedItem)
      setSelectedItem(dynamic_cast<TreeItem*>(item));
    if(_selectedItem->treeElement()->canComment())
      _selectedItem->editComment();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::barFileExternallyChangedReopen()
{
  if(fileOpen(fileName()))
    barFileExternallyChanged->setHidden(true);
}

/***********************************************/

void Tree::barUnknownElementsExpand()
{
  // recursive call
  std::function<bool(TreeElement*)> expand = [&expand](TreeElement *element)
  {
    if(dynamic_cast<TreeElementComplex*>(element))
      for(auto &child : dynamic_cast<TreeElementComplex*>(element)->children())
        if(expand(child))
          return true;
    if(!element->item() || (!dynamic_cast<TreeElementUnknown*>(element) && !element->isSelectionUnknown(element->selectedIndex())))
      return false;
    TreeElementComplex *parentElement = element->parentElement;
    while(parentElement && parentElement->item())
    {
      parentElement->item()->setExpanded(true);
      parentElement = parentElement->parentElement;
    }
    return true;
  };

  expand(rootElement);
}

/***********************************************/

void Tree::barUnknownElementsRemoveAll()
{
  try
  {
    if(QMessageBox::question(this , tr("Remove all unknown elements - GROOPS"),
                             tr("Do you really want to remove all unknown elements?\n\nNotice: This will not affect unknown choices or programs."),
                             QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Ok) != QMessageBox::Ok)
      return;

    // recursive call
    std::function<void(TreeElement*)> removeUnknown = [&removeUnknown](TreeElement *element)
    {
      if(dynamic_cast<TreeElementUnknown*>(element) && element->parentElement)
        element->parentElement->removeChild(element);
      else if(dynamic_cast<TreeElementComplex*>(element))
        for(UInt i=dynamic_cast<TreeElementComplex*>(element)->children().size(); i-->0;)
          removeUnknown(dynamic_cast<TreeElementComplex*>(element)->children().at(i));
    };

    removeUnknown(rootElement);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::barSchemaRenamedElementsExpand()
{
  // recursive call
  std::function<bool(TreeElement*)> expand = [&expand](TreeElement *element)
  {
    if(dynamic_cast<TreeElementComplex*>(element))
      for(auto &child : dynamic_cast<TreeElementComplex*>(element)->children())
        if(expand(child))
          return true;
    if(!element->item() || (!element->isRenamedInSchema() && !element->isSelectionRenamedInSchema(element->selectedIndex())))
      return false;
    TreeElementComplex *parentElement = element->parentElement;
    while(parentElement && parentElement->item())
    {
      parentElement->item()->setExpanded(true);
      parentElement = parentElement->parentElement;
    }
    return true;
  };

  expand(rootElement);
}

/***********************************************/

void Tree::barSchemaRenamedElementsUpdateAll()
{
  // recursive call
  std::function<void(TreeElement*)> updateName = [&updateName](TreeElement *element)
  {
    if(dynamic_cast<TreeElementComplex*>(element))
      for(auto &child : dynamic_cast<TreeElementComplex*>(element)->children())
        updateName(child);
    if(element->canUpdateName())
      element->updateName();
  };

  updateName(rootElement);
}

/***********************************************/

void Tree::barClickedIgnore()
{
  QFrame *bar = dynamic_cast<QFrame*>(sender()->parent());
  if(bar)
    bar->setHidden(true);
}

/***********************************************/
/**** TreeWidget *******************************/
/***********************************************/

TreeWidget::TreeWidget(Tree *tree, TabEnvironment *tabEnvironment)
  : QTreeWidget(tree), tree(tree), tabEnvironment(tabEnvironment), dragElement(nullptr)
{
  try
  {
    setlocale(LC_NUMERIC, "en_US.UTF-8"); // force . as decimal separator (QLocale behavior is strange)

    // init QTreeWidget
    // ----------------
    QStringList headerLabel;
    setHeaderLabels(headerLabel<<tr("Type")<<tr("Value")<<tr("Description")<<tr("Comment"));
    setAlternatingRowColors(true);
    setRootIsDecorated(false);
    setItemsExpandable(true);
    setSortingEnabled(false);
    setContextMenuPolicy(Qt::CustomContextMenu);
    setSelectionMode(QAbstractItemView::SingleSelection);
    header()->setStretchLastSection(true);
    header()->setSectionsMovable(false);
    viewport()->setAcceptDrops(true); // internal & external drag

    // set rectangular icon size
    // -------------------------
    QLabel *label = new QLabel(this);
    int height = label->sizeHint().height();
    setIconSize(QSize(2*height, height));
    delete label;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeWidget::eventFilter(QObject *obj, QEvent *event)
{
  Qt::KeyboardModifiers modifiers = QApplication::keyboardModifiers();

  if(event->type() == QEvent::KeyPress)
  {
    QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
    if(keyEvent->key() == Qt::Key_Tab || keyEvent->key() == Qt::Key_Backtab  ||
      (modifiers == Qt::ControlModifier && keyEvent->key() == Qt::Key_Space))
    {
      keyPressEvent(keyEvent);
      return true;
    }
  }

  if(event->type() == QEvent::Wheel && modifiers != Qt::ControlModifier)
  {
    QComboBox *combo = qobject_cast<QComboBox*>(obj);
    if(combo && !combo->hasFocus())
      return true;
  }

  // standard event processing
  return QObject::eventFilter(obj, event);
}

/***********************************************/

// Key events
// ----------
void TreeWidget::keyPressEvent(QKeyEvent *event)
{
  try
  {
    Qt::KeyboardModifiers modifiers = QApplication::keyboardModifiers();
    QModelIndex index = currentIndex();
    QModelIndex newIndex;
    const int columnEditor = 1;

    auto setNewIndex = [&] ()
    {
      if(newIndex.isValid())
      {
        setCurrentIndex(newIndex);
        TreeItem *item = dynamic_cast<TreeItem*>(currentItem());
        if(item)
          item->setFocus();
      }
    };

    // Enter: Switch focus from tree to input field of selected row
    if(event->key() == Qt::Key_Enter || event->key() == Qt::Key_Return)
    {
      newIndex = index.sibling(index.row(), columnEditor);
      setNewIndex();
    }

    // Escape: Switch focus from input field back to tree
    else if(event->key() == Qt::Key_Escape)
      setFocus();

    // Ctrl+Space: Interact with the element (fileName/program: open dialog, time: switch focus)
    else if(modifiers.testFlag(Qt::ControlModifier) && event->key() == Qt::Key_Space)
    {
      TreeItem *item = dynamic_cast<TreeItem*>(currentItem());
      if(item)
        item->treeElement()->interact();
    }

    // Ctrl+Tab: Next tab
    else if(event->matches(QKeySequence::NextChild))
      tabEnvironment->setCurrentIndex((tabEnvironment->currentIndex()+1) % tabEnvironment->count());
    // Ctrl+Shift+Tab: Previous tab
    else if(event->matches(QKeySequence::PreviousChild))
      tabEnvironment->setCurrentIndex((tabEnvironment->currentIndex()-1+tabEnvironment->count()) % tabEnvironment->count());

    // Tab: Next sibling element (or next sibling of parent if there is no next sibling, or next child otherwise)
    else if(event->key() == Qt::Key_Tab && !event->matches(QKeySequence::NextChild))
    {
      if(index.column() != columnEditor)
        newIndex = index.sibling(index.row(), columnEditor);
      else
      {
        newIndex = index.sibling(index.row()+1, columnEditor);
        if(!newIndex.isValid())
        {
          newIndex = index.parent().sibling(index.parent().row()+1, columnEditor);
          if(!newIndex.isValid())
            newIndex = model()->index(0, columnEditor, index);
        }
      }
      setNewIndex();
    }

    // Shift+Tab: Previous sibling element (or parent if there is no previous sibling)
    else if(event->key() == Qt::Key_Backtab && !event->matches(QKeySequence::PreviousChild))
    {
      if(index.row() == 0)
        newIndex = index.parent().sibling(index.parent().row(), 1);
      else
        newIndex = index.sibling(index.row()-1, 1);
      setNewIndex();
    }

    // Ctrl+Shift+Up/Down: Move element
    else if(modifiers.testFlag(Qt::ControlModifier) && modifiers.testFlag(Qt::ShiftModifier) && (event->key() == Qt::Key_Up || event->key() == Qt::Key_Down))
    {
      TreeItem *item = dynamic_cast<TreeItem*>(currentItem());
      if(item && item->treeElement()->parentElement)
      {
        TreeElement *element = item->treeElement();
        auto children = element->parentElement->children();
        auto iter = std::find(children.begin(), children.end(), element);
        if((event->key() == Qt::Key_Up) && (iter != children.begin()))
          element->parentElement->moveChild(*(--iter), element);
        else if((event->key() == Qt::Key_Down) && (++iter != children.end()) && (++iter != children.end()))
          element->parentElement->moveChild(*iter, element);
      }
    }

    // other key
    else
      QTreeWidget::keyPressEvent(event);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

// Drag & Drop events
// ------------------
void TreeWidget::mousePressEvent(QMouseEvent *event)
{
  try
  {
    if(event->button() == Qt::LeftButton)
      dragStartPosition = event->pos();
    QTreeWidget::mousePressEvent(event); // the orginal event handler
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeWidget::mouseMoveEvent(QMouseEvent *event)
{
  try
  {
    // is there a movement with pressed left button?
    if((!(event->buttons() & Qt::LeftButton)) ||
       ((event->pos()-dragStartPosition).manhattanLength() < QApplication::startDragDistance()))
    {
      QTreeWidget::mouseMoveEvent(event);  // the orginal event handler
      return;
    }

    // is position not over an item?
    TreeItem *item = dynamic_cast<TreeItem*>(itemAt(dragStartPosition));
    if(!item || !item->treeElement()->parentElement)
    {
      QTreeWidget::mouseMoveEvent(event);  // the orginal event handler
      return;
    }

    // start drag
    // ----------
    QMimeData *mimeData = createMimeData(item->treeElement());
    if(!mimeData)
    {
      QTreeWidget::mouseMoveEvent(event);  // the orginal event handler
      return;
    }

    // perform drag
    dragElement = item->treeElement();
    QDrag *drag = new QDrag(this);
    drag->setMimeData(mimeData);
    Qt::DropAction result;
    if(dragElement->parentElement->canRemoveChild(dragElement))
      result = drag->exec(Qt::MoveAction|Qt::CopyAction, Qt::MoveAction);
    else
      result = drag->exec(Qt::CopyAction);

    // delete element from source if its moved
    if(dragElement && (result == Qt::MoveAction))
      dragElement->parentElement->removeChild(dragElement);
    dragElement = nullptr;

    event->accept();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeWidget::dragEnterEvent(QDragEnterEvent *event)
{
  try
  {
    // file names?
    if(event->mimeData()->hasUrls())
    {
      event->acceptProposedAction();
      return;
    }

    if(event->mimeData()->hasFormat(mimeFormatString))
      event->acceptProposedAction();
    else
      event->ignore();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeWidget::dragMoveEvent(QDragMoveEvent *event)
{
  try
  {
    // file names?
    if(event->mimeData()->hasUrls())
    {
      event->acceptProposedAction();
      return;
    }

    // is position over an item?
    TreeItem *item = dynamic_cast<TreeItem*>(itemAt(event->pos()));
    XmlNodePtr xmlNode;
    QString    type;
    if(item && fromMimeData(event->mimeData(), xmlNode, type) &&
       item->treeElement()->parentElement && item->treeElement()->parentElement->canAddChild(item->treeElement(), type) &&
       !((event->proposedAction() == Qt::CopyAction) && !(event->keyboardModifiers() & Qt::ControlModifier)))
    {
      // check if element tried to move into itself or its children
      if(dragElement && (event->source() == this) && (event->possibleActions() & Qt::MoveAction))
      {
        TreeElement *parent = item->treeElement();
        while(parent && (parent != dragElement))
          parent = parent->parentElement;
        if(parent == dragElement)
        {
          event->ignore();
          return;
        }
      }

      event->acceptProposedAction();
      return;
    }

    event->ignore();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeWidget::dropEvent(QDropEvent *event)
{
  try
  {
    // file names?
    if(event->mimeData()->hasUrls())
    {
      QList<QUrl> urls = event->mimeData()->urls();
      if((urls.size() == 1) && (event->keyboardModifiers() == Qt::ShiftModifier))
        tree->fileOpen(urls.at(0).toLocalFile()); // replace current tab
      else
        for(int i=0; i<urls.size(); i++)
          tabEnvironment->fileOpen(urls.at(i).toLocalFile());
      event->acceptProposedAction();
      return;
    }

    // is position over an item?
    TreeItem *item = dynamic_cast<TreeItem*>(itemAt(event->pos()));
    XmlNodePtr xmlNode;
    QString    type;
    if(item && item->treeElement()->parentElement && fromMimeData(event->mimeData(), xmlNode, type))
    {
      TreeElement *targetElement = item->treeElement();
      // can directly moved?
      if(dragElement && (event->source() == this) && (event->dropAction() & Qt::MoveAction))
      {
        // move within same unbounded list
        if(targetElement->parentElement->canMoveChild(targetElement, dragElement))
        {
          targetElement->parentElement->moveChild(targetElement, dragElement);
          dragElement = nullptr;
          event->acceptProposedAction();
          return;
        }

        tree->undoStack->beginMacro("move "+dragElement->name());
        dragElement->parentElement->removeChild(dragElement);
        targetElement->parentElement->addChild(targetElement, type, xmlNode);
        tree->undoStack->endMacro();
        dragElement = nullptr;
        event->acceptProposedAction();
        return;
      }

      // external source
      if(targetElement->parentElement->addChild(targetElement, type, xmlNode))
      {
        event->acceptProposedAction();
        return;
      }
    }
    event->ignore();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
