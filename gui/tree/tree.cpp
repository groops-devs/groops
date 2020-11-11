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
#include <QLabel>
#include <QComboBox>
#include <QContextMenuEvent>
#include <QTextStream>
#include <QDomDocument>
#include <QUrl>
#include <QDesktopServices>
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
#include "addGlobalDialog/addGlobalDialog.h"
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

Tree::Tree(QWidget *parent, ActionList *actionList, TabEnvironment *workspace) : QTreeWidget(parent)
{
  try
  {
    this->settings          = new QSettings(this);
    this->changed           = false;
    this->_showResults      = true;
    this->selectedItem      = nullptr;
    this->_rootElement      = nullptr;
    this->_elementGlobal    = nullptr;
    this->fileWatcher       = nullptr;
    this->_undoStack        = new QUndoStack(this);
    this->workspace         = workspace;

    setlocale(LC_NUMERIC, "en_US.UTF-8"); // force . as decimal separator (QLocale behavior is strange)

    // init QTreeWidget
    // ----------------
    QStringList headerLabel;
    setHeaderLabels(headerLabel << tr("Type") << tr("Value") << tr("Description") << tr("Comment"));
    setAlternatingRowColors(true);
    setRootIsDecorated(false);
    setItemsExpandable(true);
    setSortingEnabled(false);
    setContextMenuPolicy(Qt::CustomContextMenu);
    setSelectionMode(QAbstractItemView::SingleSelection);
    header()->setSectionsMovable(false);
    viewport()->setAcceptDrops(true); // internal & external drag

    // set rectangular icon size
    // -------------------------
    QLabel *label = new QLabel(this);
    int height = label->sizeHint().height();
    setIconSize(QSize(2*height, height));
    delete label;

    // connection to actions
    // ---------------------
    this->actionList = *actionList;
    connect(this->actionList.editCutAction,             SIGNAL(triggered(bool)), this, SLOT(editCut()));
    connect(this->actionList.editCopyAction,            SIGNAL(triggered(bool)), this, SLOT(editCopy()));
    connect(this->actionList.editPasteAction,           SIGNAL(triggered(bool)), this, SLOT(editPaste()));
    connect(this->actionList.editPasteOverwriteAction,  SIGNAL(triggered(bool)), this, SLOT(editPasteOverwrite()));
    connect(this->actionList.editAddAction,             SIGNAL(triggered(bool)), this, SLOT(editAdd()));
    connect(this->actionList.editRemoveAction,          SIGNAL(triggered(bool)), this, SLOT(editRemove()));
    connect(this->actionList.editSetGlobalAction,       SIGNAL(triggered(bool)), this, SLOT(editSetGlobal()));
    connect(this->actionList.editSetExternalLinkAction, SIGNAL(triggered(bool)), this, SLOT(editSetExternalLink()));
    connect(this->actionList.editSetLoopAction,         SIGNAL(triggered(bool)), this, SLOT(editSetLoop()));
    connect(this->actionList.editRemoveLoopAction,      SIGNAL(triggered(bool)), this, SLOT(editRemoveLoop()));
    connect(this->actionList.editSetConditionAction,    SIGNAL(triggered(bool)), this, SLOT(editSetCondition()));
    connect(this->actionList.editRemoveConditionAction, SIGNAL(triggered(bool)), this, SLOT(editRemoveCondition()));
    connect(this->actionList.editEnabledAction,         SIGNAL(triggered(bool)), this, SLOT(editEnabled(bool)));
    connect(this->actionList.editEnableAllAction,       SIGNAL(triggered(bool)), this, SLOT(editEnableAll()));
    connect(this->actionList.editDisableAllAction,      SIGNAL(triggered(bool)), this, SLOT(editDisableAll()));
    connect(this->actionList.editRenameAction,          SIGNAL(triggered(bool)), this, SLOT(editRename()));
    connect(this->actionList.editUpdateNameAction,      SIGNAL(triggered(bool)), this, SLOT(editUpdateName()));
    connect(this->actionList.editCommentAction,         SIGNAL(triggered(bool)), this, SLOT(editComment()));
    connect(this->actionList.editCollapseAllAction,     SIGNAL(triggered(bool)), this, SLOT(editCollapseAll()));
    connect(this->actionList.openExternallyAction,      SIGNAL(triggered(bool)), this, SLOT(openExternally()));

    // QTreeWidget events
    // ------------------
    connect(QApplication::clipboard(), SIGNAL(dataChanged()),                                   this, SLOT(treeClipboardDataChanged()));
    connect(qobject_cast<QTreeWidget*>(this), SIGNAL(customContextMenuRequested(const QPoint &)),             this, SLOT(treeContextMenuRequested (const QPoint&)));
    connect(qobject_cast<QTreeWidget*>(this), SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)), this, SLOT(treeCurrentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)));
    connect(qobject_cast<QTreeWidget*>(this), SIGNAL(itemClicked       (QTreeWidgetItem*, int)),              this, SLOT(treeItemClicked       (QTreeWidgetItem*, int)));
    connect(qobject_cast<QTreeWidget*>(this), SIGNAL(itemDoubleClicked (QTreeWidgetItem*, int)),              this, SLOT(treeItemDoubleClicked (QTreeWidgetItem*, int)));

    // track state of file
    connect(_undoStack,  SIGNAL(cleanChanged(bool)), this, SLOT(treeCleanChanged(bool)));

    connect(this, SIGNAL(treeSelectionChanged(const QString &)), this->workspace, SIGNAL(treeSelectionChanged(const QString &)));
    connect(header(), &QHeaderView::sectionResized, this, &Tree::resizeColumn);

    newFile();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Tree::~Tree()
{
  try
  {
    {
      const QSignalBlocker blocker(this);
      clearTree();
    }
  }
  catch(std::exception &e)
  {
    qDebug() << QString::fromStdString("Exception in destructor at "+_GROOPS_ERRORLINE+"\n"+e.what());
  }
}

/***********************************************/

void Tree::clearTree()
{
  try
  {
    setSelectedItem(nullptr);
    changed = false;
    _undoStack->clear();

    if(rootElement())
    {
      rootElement()->removeItem();
      delete _rootElement;
    }
    _rootElement   = nullptr;
    _elementGlobal = nullptr;
    _varList.clear();
    unknownElements.clear();
    renamedElements.clear();
    emit unknownElementsChanged(unknownElements.size());
    emit renamedElementsChanged(renamedElements.size());
    clearFileWatcher();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::setSelectedItem(TreeItem *item)
{
  try
  {
    if(item == selectedItem)
      return;

    // old item lost selection
    if(selectedItem)
      selectedItem->lostCurrent();

    // new item gets selection
    selectedItem = item;
    if(selectedItem)
    {
      selectedItem->becomeCurrent();
      updateActions();
      emit treeSelectionChanged(selectedItem->treeElement()->selectedValue());
    }

    if(selectedItem != QTreeWidget::currentItem())
      QTreeWidget::setCurrentItem(item);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::setShowDescriptions(Bool state)
{
  setColumnHidden(2, !state);
  if(!isColumnHidden(2))
    setColumnWidth(2, std::max(columnWidth(2), 100));
}

/***********************************************/

void Tree::setShowResults(Bool state)
{
  _showResults = state;
  if(rootElement())
    rootElement()->createItem(nullptr, nullptr);
}

/***********************************************/

void Tree::updateExpressions(const TreeElement *element)
{
  try
  {
    if(!element)
      return;

    if(element->parentElement && element->parentElement == elementGlobal())
    {
      for(auto &&var : _varList)
        var.second->setValue(var.second->getText()); // reset as TEXT
      dynamic_cast<TreeElementComplex*>(rootElement())->updateExpression();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

std::vector<XsdElementPtr> Tree::programListFromSchema() const
{
  try
  {
    return schema.programList();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::addProgram(int index)
{
  if(programType().isEmpty())
    return;

  TreeElementComplex *root = dynamic_cast<TreeElementComplex*>(rootElement());
  root->addChild(root->childAt(root->childrenCount()-1), programType(), XmlNodePtr(nullptr));
  root->childAt(root->childrenCount()-2)->changeSelectedIndex(index);
}

void Tree::resizeColumn(int /*logicalIndex*/, int /*oldSize*/, int /*newSize*/)
{
  if(this != workspace->currentTree())
    return;

  std::vector<int> columnWidths;
  for(int i = 0; i < 3; i++)
    columnWidths.push_back(columnWidth(i));
  workspace->resizeTreeColumns(columnWidths);
}

/***********************************************/

// Wird die Datei zum Abbruch freigegeben?
Bool Tree::okToAbandon()
{
  try
  {
    if(!isChanged())
      return true;

    QString name = xmlFile;
    if(name.isEmpty())
      name = workspace->currentTabText().toUtf8();

    QMessageBox::StandardButton button =
      QMessageBox::question(this , tr("Close File - GROOPS"),
                            tr("File '%1' was changed.\nDo you want to save the file?").arg(name),
                            QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel, QMessageBox::Save);
    if(button == QMessageBox::Save)
      return saveFile();
    else if(button == QMessageBox::Discard)
      return true;
    return false;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool Tree::readSchema()
{
  try
  {
    QString schemaFile = settings->value("files/schemaFile").toString();
    while(Schema::validateSchema(schemaFile) == false)
    {
        if(schemaFile.isEmpty())
          QMessageBox::information(this , tr("GROOPS"), tr("GROOPS seems not to be configured yet. You should set at least the XSD schema file."));
        else
          QMessageBox::critical(this , tr("GROOPS"), tr("File '%1' seems not to be a valid XSD schema").arg(schemaFile));
        SettingsPathDialog dialog(this);
        if(dialog.exec())
          emit workspace->schemaChanged();
        else
          return false;
        schemaFile = settings->value("files/schemaFile").toString();
    }
    schema = Schema(schemaFile);
    _programType = schema.programType();
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool Tree::newFile()
{
  try
  {
    if(!okToAbandon())
      return false;

    if(!readSchema())
      return false;

    if(xmlFile.isEmpty())
      xmlDir.setPath(settings->value("files/workingDirectory").toString());
    xmlFile = QString();
    xsdFile = settings->value("files/schemaFile").toString();
    clearTree();
    emit fileChanged(fileName(), changed);

    QString templateFile = settings->value("files/templateFile").toString();
    if(QFileInfo(templateFile).isFile() && openFile(templateFile))
    {
      xmlFile = QString();
      xmlDir.setPath(settings->value("files/workingDirectory").toString());
      return true;
    }

    _rootElement = TreeElement::newTreeElement(this, nullptr, schema.rootElement, "", XmlNodePtr(nullptr), false);
    TreeItem *item = rootElement()->createItem(nullptr, nullptr);
    if(item)
      item->setExpanded(true);
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool Tree::reopenFile()
{
  try
  {
    if(!xmlFile.isEmpty())
      return openFile(QFileInfo(xmlDir, xmlFile).absoluteFilePath());
    return openFile();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool Tree::openFile(QString fileName)
{
  try
  {
    if(!okToAbandon())
      return false;

    if(fileName.isEmpty() && !xmlFile.isEmpty())
      fileName = QFileInfo(xmlDir, xmlFile).absoluteFilePath();

    // read xml file
    // -------------
    XmlNodePtr xmlNode;
    if(!fileName.isEmpty())
    {
      QFile file(fileName);
      if(!file.open(QFile::ReadOnly | QFile::Text))
      {
        QMessageBox::critical(this , tr("Open file - GROOPS"), tr("Cannot open file '%1'.").arg(fileName));
        return false;
      }
      int          errRow, errCol;
      QString      errStr;
      QDomDocument doc;
      if(!doc.setContent(&file, false, &errStr, &errRow, &errCol))
      {
        QMessageBox::critical(this , tr("Open file - GROOPS"),
                              tr("Error Reading File '%1' (row=%2, col=%3):\n%4").arg(fileName).arg(errRow).arg(errCol).arg(errStr));
        return false;
      }
      xmlNode = XmlNode::create(doc.documentElement());
    }

    if(!readSchema())
      return false;

    // seems to be all ok
    // ------------------
    xmlFile = QFileInfo(fileName).fileName();
    xmlDir  = QFileInfo(fileName).dir();
    xsdFile = settings->value("files/schemaFile").toString();;
    clearTree();
    emit fileChanged(this->fileName(), changed);
    _rootElement = TreeElement::newTreeElement(this, nullptr, schema.rootElement, "", xmlNode, true);
    rootElement()->createItem(nullptr, nullptr)->setExpanded(true);
    createFileWatcher();
    for(auto &&var : _varList)
      var.second->setValue(var.second->getText()); // reset as TEXT
    dynamic_cast<TreeElementComplex*>(rootElement())->updateExpression();

    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool Tree::saveFile()
{
  try
  {
    if(xmlFile.isEmpty())
      return saveAsFile();

    clearFileWatcher();
    XmlNodePtr xmlNode = rootElement()->getXML();
    XmlNode::writeFile(fileName(), xmlNode);
    _undoStack->setClean();
    treeCleanChanged(true);
    emit treeFileChanged(fileName(), false);
    createFileWatcher();
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool Tree::saveAsFile(const QString &fileName)
{
  try
  {
    // get file name: open file selector?
    // ----------------------------------
    QString name = fileName;
    if(name.isEmpty())
    {
      // Open file selector
      if(xmlFile.isEmpty())
        xmlDir.setPath(settings->value("files/workingDirectory").toString());
      name = QFileInfo(xmlDir, xmlFile).absoluteFilePath();
      name = QFileDialog::getSaveFileName(this, tr("Save file - GROOPS"), name, tr("XML files (*.xml)"));
      if(name.isEmpty())
        return false;
      if(!name.endsWith(".xml"))
        name += ".xml";
    }

    xmlFile = QFileInfo(name).fileName();
    xmlDir  = QFileInfo(name).dir();

    return saveFile();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool Tree::execFile()
{
  try
  {
    // file must be saved first
    // ------------------------
    if(isChanged() || xmlFile.isEmpty())
    {
      QString name = xmlFile;
      if(name.isEmpty())
        name = workspace->currentTabText().toUtf8();

      QMessageBox::StandardButton button =
        QMessageBox::warning(this , tr("Run file - GROOPS"),
                             tr("File '%1' was changed.\nYou have to save it first.").arg(name),
                             QMessageBox::Save | QMessageBox::Cancel, QMessageBox::Save);
      if(button != QMessageBox::Save)
        return false;
      if(!saveFile())
        return false;
    }

    // execute dialog
    // --------------
    if(!rootElement())
      return false;
    ExecuteDialog dialog(this, this);
    if(!dialog.exec())
      return false;

    // create execute command
    // ----------------------
    QStringList commandList  = settings->value("execute/commands").toStringList();
    int         commandIndex = settings->value("execute/commandIndex", int(0)).toInt();
    if((commandIndex<0)||(commandIndex>=commandList.size()))
      return false;
    QString command = commandList[commandIndex];
    if(command.isEmpty())
      return false;

    QString optionString;
    // append log file option?
    if(settings->value("execute/useLogFile", bool(false)).toBool())
      optionString += " -l "+settings->value("execute/logFile", QString("groops.log")).toString()+" ";
    optionString += xmlFile;

    command.replace("%f", optionString);
    command.replace("%w", xmlDir.absolutePath());

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

    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool Tree::isChanged() const
{
  return changed;
}

/***********************************************/
/***********************************************/

void Tree::createFileWatcher()
{
  try
  {
    clearFileWatcher();
    if(fileName().isEmpty())
      return;
    fileWatcher = new QFileSystemWatcher(this);
    fileWatcher->addPath(QFileInfo(xmlDir, xmlFile).absoluteFilePath());
    connect(fileWatcher, &QFileSystemWatcher::fileChanged, this, &Tree::fileChangedExternally);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::clearFileWatcher()
{
  if(fileWatcher)
  {
    disconnect(fileWatcher, &QFileSystemWatcher::fileChanged, this, &Tree::fileChangedExternally);
    delete fileWatcher;
    fileWatcher = nullptr;
  }
}

/***********************************************/

void Tree::fileChangedExternally()
{
  try
  {
    emit treeFileChanged(fileName(), true);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void Tree::trackUnknownElement(TreeElement *element)
{
  try
  {
    if(!element || element->isElementAdd() || unknownElements.contains(element))
      return;

    unknownElements.insert(element);
    emit unknownElementsChanged(unknownElements.size());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool Tree::untrackUnknownElement(TreeElement *element)
{
  try
  {
    if(!element || !unknownElements.remove(element))
      return false;

    emit unknownElementsChanged(unknownElements.size());
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::expandUnknownElements()
{
  try
  {
    for(auto &&element : unknownElements)
      if(element && element->item())
      {
        TreeElementComplex* parentElement = element->parentElement;
        while(parentElement && parentElement->item())
        {
          if(!parentElement->isUnknown() && !parentElement->isSelectionUnknown(parentElement->selectedIndex()))
            parentElement->item()->setExpanded(true);
          parentElement = parentElement->parentElement;
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::removeAllUnknownElements()
{
  try
  {
    QMessageBox::StandardButton button =
        QMessageBox::question(this , tr("Remove all unknown elements - GROOPS"),
                            tr("Do you really want to remove all unknown elements?\n\nNotice: This will not affect unknown choices or programs."),
                            QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Ok);
    if(button != QMessageBox::Ok)
      return;

    QTreeWidgetItemIterator it(this);
    while(*it)
    {
      TreeItem *item = dynamic_cast<TreeItem*>(*it);
      if(item && item->treeElement())
      {
        if(item->treeElement()->isUnknown())
        {
          bool hasParentWithUnknownSelection = false;
          TreeElementComplex* parentElement = item->treeElement()->parentElement;
          while(parentElement && parentElement->item())
          {
            if(!parentElement->isUnknown() && parentElement->isSelectionUnknown(parentElement->selectedIndex()))
            {
              hasParentWithUnknownSelection = true;
              break;
            }
            parentElement = parentElement->parentElement;
          }

          if(!hasParentWithUnknownSelection && item->treeElement()->parentElement->removeChild(item->treeElement()))
            it = QTreeWidgetItemIterator(this); // restart at beginning to revalidate iterators
        }
      }

      ++it;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void Tree::trackRenamedElement(TreeElement *element)
{
  try
  {
    if(!element || renamedElements.contains(element) || (!element->isRenamed() && !element->isSelectionRenamed(element->selectedIndex())))
      return;

    renamedElements.insert(element);
    emit renamedElementsChanged(renamedElements.size());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool Tree::untrackRenamedElement(TreeElement *element)
{
  try
  {
    if(!element || !renamedElements.remove(element))
      return false;

    emit renamedElementsChanged(renamedElements.size());
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::expandRenamedElements()
{
  try
  {
    for(auto &&element : renamedElements)
      if(element && element->item())
      {
        TreeElementComplex* parentElement = element->parentElement;
        while(parentElement && parentElement->item())
        {
          parentElement->item()->setExpanded(true);
          parentElement = parentElement->parentElement;
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::updateAllRenamedElements()
{
  try
  {
    for(auto &&element : renamedElements.values())
      if(element && element->parentElement)
        element->updateName();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

QString Tree::fileName() const
{
  if(xmlFile.isEmpty())
    return QString();
  return QFileInfo(xmlDir, xmlFile).absoluteFilePath();
}

/***********************************************/


QString Tree::addXmlDirectory(const QString &filename) const
{
  try
  {
    QFileInfo file(filename);
    QString s = file.filePath();
    if(file.isRelative())
      file.setFile(xmlDir, s);
    return file.absoluteFilePath();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QString Tree::stripXmlDirectory(const QString &filename) const
{
  try
  {
    QString path = xmlDir.absolutePath();
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
/***********************************************/

static const char *mimeFormatString = "application/x-groops"; //"text/plain";

/***********************************************/

QMimeData *Tree::createMimeData(const TreeElement *element)
{
  try
  {
    if(element->type().isEmpty())
      return nullptr;
    XmlNodePtr xmlNode = element->getXML(true);
    if(xmlNode==nullptr)
      return nullptr;
    writeAttribute(xmlNode, "xsdType", element->type());

    // create xml text
    QString xmlData;
    QTextStream stream(&xmlData, QIODevice::WriteOnly);
    XmlNode::write(stream, xmlNode);

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

Bool Tree::fromMimeData(const QMimeData *mimeData, XmlNodePtr &xmlNode, QString &type)
{
  try
  {
    if(!mimeData->hasFormat(mimeFormatString))
      return false;
    // create xmlData from MIME data
    QByteArray xmlData = mimeData->data(mimeFormatString);
    QDomDocument doc;
    if(!doc.setContent(xmlData))
      return false;
    xmlNode = XmlNode::create(doc.documentElement());
    if(!xmlNode)
      return false;
    readAttribute(xmlNode, "xsdType", type, false);
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::updateActions()
{
  try
  {
    if(!selectedItem)
      return;
    TreeElement        *element       = selectedItem->treeElement();
    TreeElementComplex *parentElement = element->parentElement;


    // test content of clipboard
    QString    type;
    XmlNodePtr xmlNode;
    Bool isClipboard = fromMimeData(QApplication::clipboard()->mimeData(), xmlNode, type);
    Bool canAdd      = element == elementGlobal() || (parentElement && parentElement->canAddChild(element, element->type()));
    Bool canRemove   = parentElement && parentElement->canRemoveChild(element);

    actionList.editCutAction->setEnabled( canRemove && (!element->type().isEmpty()) && (element->getXML(true)!=nullptr) );
    actionList.editCopyAction->setEnabled( (!element->type().isEmpty()) && (element->getXML(true)!=nullptr) );
    actionList.editPasteAction->setEnabled( isClipboard && canAdd );
    actionList.editPasteOverwriteAction->setEnabled( isClipboard && (type==element->type()) && (!element->isElementAdd()) );
    actionList.editAddAction->setEnabled( canAdd );
    actionList.editRemoveAction->setEnabled( canRemove );
    actionList.editSetGlobalAction->setEnabled( elementGlobal() && element != elementGlobal() && elementGlobal()->canSetGlobal(element) );
    actionList.editSetExternalLinkAction->setEnabled( elementGlobal() && element != elementGlobal() && elementGlobal()->canSetGlobal(element)  );
    actionList.editSetLoopAction->setEnabled( element->canSetLoop() );
    actionList.editRemoveLoopAction->setEnabled( !element->loop().isEmpty() );
    actionList.editSetConditionAction->setEnabled( element->canSetCondition() );
    actionList.editRemoveConditionAction->setEnabled( !element->condition().isEmpty() );
    actionList.editEnabledAction->setEnabled( element->canDisabled() );
    actionList.editEnabledAction->setChecked( !element->disabled() );
    actionList.editEnableAllAction->setEnabled( true );
    actionList.editDisableAllAction->setEnabled( true );
    actionList.editRenameAction->setEnabled( element->canRename());
    actionList.editUpdateNameAction->setEnabled( element->canUpdateName());
    actionList.editCommentAction->setEnabled(true);
    actionList.editCollapseAllAction->setEnabled( true );
    TreeElementFileName *fileNameElement = dynamic_cast<TreeElementFileName*>(element);
    actionList.openExternallyAction->setEnabled( fileNameElement && QFileInfo(addXmlDirectory(fileNameElement->selectedResult())).isFile() );
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

// slots
// -----

void Tree::editCut()
{
  try
  {
    if(!selectedItem)
      return;
    TreeElement        *element       = selectedItem->treeElement();
    TreeElementComplex *parentElement = element->parentElement;
    if(parentElement && parentElement->canRemoveChild(element))
    {
      // copy to clipboard
      QMimeData *mimeData = createMimeData(element);
      if(mimeData)
      {
        QApplication::clipboard()->setMimeData(mimeData);
        parentElement->removeChild(element);
      }
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
    if(!selectedItem)
      return;
    // copy to clipboard
    QMimeData *mimeData = createMimeData(selectedItem->treeElement());
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    // paste from clipboard
    QString    type;
    XmlNodePtr xmlNode;
    if(!fromMimeData(QApplication::clipboard()->mimeData(), xmlNode, type))
      return;
    if(element->parentElement)
      element->parentElement->addChild(element, type, xmlNode);
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    // paste from clipboard
    QString    type;
    XmlNodePtr xmlNode;
    if(!fromMimeData(QApplication::clipboard()->mimeData(), xmlNode, type))
      return;
    if(element->parentElement)
      element->parentElement->overwrite(element, type, xmlNode);
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();

    if(element == elementGlobal() || (element && element->parentElement == elementGlobal() && element->isElementAdd()))
    {
      AddGlobalDialog dialog(elementGlobal(), this);
      if(dialog.exec())
      {
        QString label = dialog.elementName();
        elementGlobal()->addChild(elementGlobal()->elementAdd(), dialog.elementType()->type, XmlNodePtr(nullptr), label);
      }
    }
    else if(element->parentElement)
      element->parentElement->addChild(element, element->type(), element->getXML());
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
    if(!selectedItem)
      return;
    TreeElement        *element       = selectedItem->treeElement();
    TreeElementComplex *parentElement = element->parentElement;
    if(parentElement && parentElement->canRemoveChild(element))
    {
      QMessageBox::StandardButton button =
          QMessageBox::question(this , tr("Remove element - GROOPS"),
                              tr("Do you really want to remove this element?"),
                              QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Ok);
      if(button != QMessageBox::Ok)
        return;

       // remove element
      parentElement->removeChild(element);
    }
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    if(!(elementGlobal() && elementGlobal()->canSetGlobal(element)))
      return;

    elementGlobal()->setGlobal(element);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editSetExternalLink()
{
  try
  {
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    if(!(elementGlobal() && elementGlobal()->canSetGlobal(element)))
      return;

    QStringList existingNames = elementGlobal()->getChildrenNames();
    existingNames.removeAt(existingNames.indexOf(element->label()));

    bool ok;
    QString label = QInputDialog::getText(this, tr("Set external link - GROOPS"), tr("Name of global variable in external file:"), QLineEdit::Normal, "", &ok);
    QRegExp regex("[a-zA-Z]([a-zA-Z0-9])*");
    while(ok && (label.isEmpty() || existingNames.contains(label) || !regex.exactMatch(label)))
      label = QInputDialog::getText(this, tr("Set external link - GROOPS"), tr("Name is already used by global variable in this file or is invalid (only letters and digits allowed)!\nChoose another name:"), QLineEdit::Normal, label, &ok);

    if(ok)
      element->setNewExternalLink(label);
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    if(!element->canSetLoop())
      return;

    SetLoopConditionDialog dialog(elementGlobal(), "loop", this);
    if(dialog.exec())
    {
      QStringList globalElements = elementGlobal()->getChildrenNames();
      QString loopName = dialog.name();
      if(!globalElements.contains(loopName, Qt::CaseInsensitive))
        elementGlobal()->addChild(nullptr, "loopType", XmlNodePtr(nullptr), loopName, true);
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    if(element->loop().isEmpty())
      return;

    element->setLoop("");
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    if(!element->canSetCondition())
      return;

    SetLoopConditionDialog dialog(elementGlobal(), "condition", this);
    if(dialog.exec())
    {
      QStringList globalElements = elementGlobal()->getChildrenNames();
      QString conditionName = dialog.name();
      if(!globalElements.contains(conditionName, Qt::CaseInsensitive))
        elementGlobal()->addChild(nullptr, "conditionType", XmlNodePtr(nullptr), conditionName, true);
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    if(element->condition().isEmpty())
      return;

    element->setCondition("");
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
    if(!selectedItem)
      return;
    selectedItem->treeElement()->setDisabled(!checked);
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
    TreeElementComplex* rootElement = dynamic_cast<TreeElementComplex*>(this->rootElement());
    if(workspace->currentTree() != this || !rootElement)
      return;

    undoStack()->beginMacro("enable all programs");
    for(int i = 0; i < rootElement->childrenCount(); i++)
      if(rootElement->childAt(i)->disabled() && rootElement->childAt(i)->isProgram())
        rootElement->childAt(i)->setDisabled(false);
    undoStack()->endMacro();
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
    TreeElementComplex* rootElement = dynamic_cast<TreeElementComplex*>(this->rootElement());
    if(workspace->currentTree() != this || !rootElement)
      return;

    undoStack()->beginMacro("disable all programs");
    for(int i = 0; i < rootElement->childrenCount(); i++)
      if(!rootElement->childAt(i)->disabled() && rootElement->childAt(i)->isProgram())
        rootElement->childAt(i)->setDisabled(true);
    undoStack()->endMacro();
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    if(!element || !element->canRename())
      return;

    QStringList existingNames = elementGlobal()->getChildrenNames();
    existingNames.removeAt(existingNames.indexOf(element->label()));

    bool ok;
    QString label = QInputDialog::getText(this, tr("Rename global element - GROOPS"), tr("New name of global element:"), QLineEdit::Normal, element->label(), &ok);
    QRegExp regex("[a-zA-Z]([a-zA-Z0-9])*");
    while(ok && (label.isEmpty() || existingNames.contains(label) || !regex.exactMatch(label)))
      label = QInputDialog::getText(this, tr("Rename global element - GROOPS"), tr("Name already exists or is invalid (only letters and digits allowed)!\nChoose another name:"), QLineEdit::Normal, label, &ok);

    if(ok)
      element->rename(label);
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
    if(!selectedItem)
      return;
    TreeElement *element = selectedItem->treeElement();
    if(!element || !element->canUpdateName())
      return;

    element->updateName();
    updateActions();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::editComment()
{
  try
  {
    if(selectedItem)
      selectedItem->editComment();
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
    if(workspace->currentTree() != this)
      return;

    collapseAll();
    expand(model()->index(0,0));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::openExternally()
{
  try
  {
    if(workspace->currentTree() != this && !selectedItem)
      return;

    QDesktopServices::openUrl(QUrl::fromLocalFile(addXmlDirectory(selectedItem->treeElement()->selectedResult())));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}


/***********************************************/
/**** Event-Handler ****************************/
/***********************************************/

void Tree::treeCleanChanged(bool clean)
{
  try
  {
    changed = !clean;
    updateActions();
    emit fileChanged(fileName(), changed);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

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
    QTreeWidgetItem *item = QTreeWidget::itemAt(pos);
    if(item==nullptr)
      return;
    setSelectedItem(dynamic_cast<TreeItem*>(item));

    QMenu *contextMenu = new QMenu(this);
    contextMenu->addAction(actionList.editCutAction);
    contextMenu->addAction(actionList.editCopyAction);
    contextMenu->addAction(actionList.editPasteAction);
    contextMenu->addAction(actionList.editPasteOverwriteAction);
    contextMenu->addSeparator();
    contextMenu->addAction(actionList.editAddAction);
    contextMenu->addAction(actionList.editRemoveAction);
    contextMenu->addSeparator();
    contextMenu->addAction(actionList.editSetGlobalAction);
    contextMenu->addAction(actionList.editSetExternalLinkAction);
    contextMenu->addSeparator();
    contextMenu->addAction(actionList.editSetLoopAction);
    contextMenu->addAction(actionList.editRemoveLoopAction);
    contextMenu->addAction(actionList.editSetConditionAction);
    contextMenu->addAction(actionList.editRemoveConditionAction);
    contextMenu->addSeparator();
    contextMenu->addAction(actionList.editEnabledAction);
    contextMenu->addAction(actionList.editRenameAction);
    contextMenu->addAction(actionList.editUpdateNameAction);
    contextMenu->addAction(actionList.editCommentAction);
    contextMenu->addAction(actionList.openExternallyAction);
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

void Tree::treeItemClicked(QTreeWidgetItem *item, int column)
{
  try
  {
    // value column?
    if(column == 1)
    {
      if(item != selectedItem)
        setSelectedItem(dynamic_cast<TreeItem*>(item));
      selectedItem->setFocus();
    }
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
    // comment column?
    if(column == 3)
    {
      if(item != selectedItem)
        setSelectedItem(dynamic_cast<TreeItem*>(item));
      selectedItem->editComment();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

// Key events
// ------------

void Tree::keyPressEvent(QKeyEvent *event)
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
        {
          item->setFocus();
        }
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
      workspace->setCurrentIndex((workspace->currentIndex()+1) % workspace->count());

    // Ctrl+Shift+Tab: Previous tab
    else if(event->matches(QKeySequence::PreviousChild))
      workspace->setCurrentIndex((workspace->currentIndex()-1+workspace->count()) % workspace->count());

    // Tab: Next sibling element (or next sibling of parent if there is no next sibling, or next child otherwise)
    else if(event->key() == Qt::Key_Tab)
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
    else if(event->key() == Qt::Key_Backtab)
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
      int idxShift = (event->key() == Qt::Key_Up) ? -1 : 1;

      TreeItem *item = dynamic_cast<TreeItem*>(currentItem());
      if(item && item->treeElement() && item->treeElement()->parentElement)
      {
        TreeElement *element = item->treeElement();
        const int childrenCount = element->parentElement->childrenCount();
        int idx = 0;
        for(; idx < childrenCount; idx++)
          if(element->parentElement->childAt(idx) == element)
            break;

        if(idx < childrenCount && idx+idxShift >= 0 && idx+idxShift < childrenCount-1 && !element->parentElement->childAt(idx+idxShift)->isElementAdd())
          element->parentElement->moveChild(element, element->parentElement->childAt(idx+idxShift));
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

void Tree::mousePressEvent(QMouseEvent *event)
{
  try
  {
    if(event->button() == Qt::LeftButton)
      dragStartPosition = event->pos();
    QTreeWidget::mousePressEvent(event); // the orginal event handler

    // catch edge case where clicked item changes because editor is removed from item above and therefore vertical position changes
    if(selectedItem && QTreeWidget::selectedItems().size() > 1)
    {
      for(const auto &item : QTreeWidget::selectedItems())
          item->setSelected(false);
      selectedItem->setSelected(true);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::mouseMoveEvent(QMouseEvent *event)
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
    if(!item)
    {
      QTreeWidget::mouseMoveEvent(event);  // the orginal event handler
      return;
    }
    TreeElement        *element       = item->treeElement();
    TreeElementComplex *parentElement = element->parentElement;

    // start drag
    // ----------
    QMimeData *mimeData = createMimeData(element);
    if(mimeData==nullptr)
    {
      QTreeWidget::mouseMoveEvent(event);  // the orginal event handler
      return;
    }
    // perform drag
    QDrag *drag = new QDrag(viewport());
    drag->setMimeData(mimeData);
    Qt::DropAction result;
    if(parentElement->canRemoveChild(element))
      result = drag->exec(Qt::MoveAction|Qt::CopyAction, Qt::MoveAction);
    else
      result = drag->exec(Qt::CopyAction);

    // delete element from source if its moved
    if(result == Qt::MoveAction && parentElement)
    {
      // check if element was moved into itself or its children
      TreeItem *targetItem = dynamic_cast<TreeItem*>(itemAt(this->mapFromGlobal(QCursor::pos())));
      TreeElement *parent = targetItem ? targetItem->treeElement()->parentElement : nullptr;
      while(parent && parent != element)
        parent = parent->parentElement;

      if(parent == element)
      {
        // abort and reverse drag and drop
        undoStack()->undo();
        undoStack()->push(new QUndoCommand);
        undoStack()->undo();
        event->ignore();
        return;
      }
      else
        parentElement->removeChild(element);
    }

    event->accept();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void Tree::dragEnterEvent(QDragEnterEvent *event)
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

void Tree::dragMoveEvent(QDragMoveEvent *event)
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
    if(item)
    {
      XmlNodePtr xmlNode;
      QString    type;
      if(fromMimeData(event->mimeData(), xmlNode, type))
      {
        TreeElement        *element       = item->treeElement();
        TreeElementComplex *parentElement = element->parentElement;
        if(parentElement && parentElement->canAddChild(element, type) &&
           !(event->proposedAction() == Qt::CopyAction && !(event->keyboardModifiers() & Qt::ControlModifier)))
        {
          event->acceptProposedAction();
          return;
        }
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

void Tree::dropEvent(QDropEvent *event)
{
  try
  {
    // file names?
    if(event->mimeData()->hasUrls())
    {
      QList<QUrl> urls = event->mimeData()->urls();
      if(urls.size() == 1 && event->keyboardModifiers() == Qt::ShiftModifier)
        openFile(urls.at(0).toLocalFile()); // replace current tab
      else
      {
        int i = 0;
        // if only a single, unchanged new tab exists, replace it with the first file to open and open the other files as additional tabs
        if(workspace->count() == 1 && workspace->currentTree() && !workspace->currentTree()->fileName().startsWith("/") && !workspace->currentTree()->isChanged())
          openFile(urls.at(i++).toLocalFile());

        for(; i < urls.size(); i++)
          workspace->openFile(urls.at(i).toLocalFile());
      }
      event->acceptProposedAction();
      return;
    }

    // is position over an item?
    TreeItem *item = dynamic_cast<TreeItem*>(itemAt(event->pos()));
    if(item)
    {
      XmlNodePtr xmlNode;
      QString    type;
      if(fromMimeData(event->mimeData(), xmlNode, type))
      {
        TreeElement        *element        = item->treeElement();
        TreeElementComplex *parentElement  = element->parentElement;

        bool moved = event->dropAction() == Qt::MoveAction;
        if(parentElement && parentElement->addChild(element, type, xmlNode, moved))
        {
          event->acceptProposedAction();
          return;
        }
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

void Tree::resizeEvent(QResizeEvent *event)
{
  Double oldWidth = event->oldSize().width();
  if(oldWidth > 0)
  {
    int newWidth = event->size().width();
    Double totalWidth = columnWidth(0) + columnWidth(1) + columnWidth(2) + columnWidth(3);
    if(totalWidth < oldWidth)
      totalWidth = oldWidth;

    std::vector<int> columnWidths;
    for(int i = 0; i < 3; i++)
      columnWidths.push_back(static_cast<int>(std::round(columnWidth(i)*newWidth/totalWidth)));

    header()->blockSignals(true);
    workspace->resizeTreeColumns(columnWidths);
    header()->blockSignals(false);
  }

  QTreeWidget::resizeEvent(event);
}

/***********************************************/

bool Tree::eventFilter(QObject *obj, QEvent *event)
{
  Qt::KeyboardModifiers modifiers = QApplication::keyboardModifiers();

  if (event->type() == QEvent::KeyPress)
  {
    QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);

    if(keyEvent->key() == Qt::Key_Tab || keyEvent->key() == Qt::Key_Backtab  ||
       (modifiers == Qt::ControlModifier && keyEvent->key() == Qt::Key_Space))
    {
      keyPressEvent(keyEvent);
      return true;
    }
  }

  if(event->type() == QEvent::Wheel && modifiers != Qt::ControlModifier)
  {
      QComboBox* combo = qobject_cast<QComboBox*>(obj);
      if(combo && !combo->hasFocus())
          return true;
  }

  // standard event processing
  return QObject::eventFilter(obj, event);
}

/***********************************************/
/***********************************************/
