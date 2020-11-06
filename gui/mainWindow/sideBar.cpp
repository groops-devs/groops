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

#include "sideBar.h"
#include "tree/tree.h"
#include <QDebug>
#include <QDesktopServices>
#include <QFileInfo>
#include <QPainter>
#include <QStyleOption>
#include <QMenu>
#include <QHeaderView>
#include <QUrl>
#include <QMimeData>

/***********************************************/

SideBar::SideBar(QWidget *parent)
{
  try
  {
    if(parent)
      this->setParent(parent);

    settings = new QSettings(this);

    layout = new QVBoxLayout(this);
    layout->setAlignment(this, Qt::AlignTop);
    layout->setSpacing(0);
    layout->setContentsMargins(0,0,0,0);

    currentIndex = -1;

    _stackedWidget = new StackedWidget(settings->value("sideBar/width", parentWidget()->width()/5).toInt());

    connect(this, SIGNAL(currentChanged(int)), _stackedWidget, SLOT(setCurrentIndex(int)));
    connect(this, SIGNAL(showHide(bool)),      _stackedWidget, SLOT(showHide(bool)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

SideBar::~SideBar()
{
  if(layout)
    delete layout;
  if(_stackedWidget)
    delete _stackedWidget;
}

/***********************************************/

void SideBar::addSideBarWidget(QString buttonLabel, QWidget *widget)
{
  addPushButton(buttonLabel);
  _stackedWidget->addWidget(widget);
}

/***********************************************/

void SideBar::setSideBarWidgetsHidden(bool hide)
{
  stackedWidget()->setHidden(hide);

  if(hide)
    for(size_t i = 0; i < buttons.size(); i++)
      buttons.at(i)->setChecked(false);
}

/***********************************************/

PushButton* SideBar::addPushButton(QString text, Qt::Orientation orientation)
{
  PushButton *button = new PushButton(text, this);

  if(orientation == Qt::Vertical)
    button->setMaximumWidth(button->minimumSizeHint().height());
  button->setOrientation(orientation);
  button->setFlat(true);
  button->setCheckable(true);

  connect(button, SIGNAL(clicked(bool)), this, SLOT(buttonClicked(bool)));

  buttons.push_back(button);
  layout->addWidget(button);

  if(currentIndex == -1 && !settings->value("sideBar/isHidden", false).toBool())
  {
    currentIndex = 0;
    button->setChecked(true);
  }

  return button;
}

/***********************************************/

void SideBar::buttonClicked(bool /*checked*/)
{
  QObject* obj = sender();
  for(size_t i = 0; i < buttons.size(); i++)
  {
    if(obj == buttons.at(i) || obj == _stackedWidget->widget(static_cast<int>(i)))
    {
      if(static_cast<int>(i) == currentIndex)
      {
        emit showHide(true);
        currentIndex = -1;
        buttons.at(i)->setChecked(false);
      }
      else
      {
        emit currentChanged(static_cast<int>(i));
        emit showHide(false);
        currentIndex = static_cast<int>(i);
        buttons.at(i)->setChecked(true);
      }
    }
    else
      buttons.at(i)->setChecked(false);
  }
}

/***********************************************/

void SideBar::showWidget()
{
  QObject* obj = sender();
  for(size_t i = 0; i < buttons.size(); i++)
  {
    if(obj == buttons.at(i) || obj == _stackedWidget->widget(static_cast<int>(i)))
    {
      emit currentChanged(static_cast<int>(i));
      emit showHide(false);
      currentIndex = static_cast<int>(i);
      buttons.at(i)->setChecked(true);
    }
    else
      buttons.at(i)->setChecked(false);
  }
}

/***********************************************/

void StackedWidget::widthChanged(int width, int /*index*/)
{
  if(!isHidden())
    _lastWidth = width;
}

/***********************************************/

void StackedWidget::showHide(bool isCurrent)
{
  if(!isHidden() && isCurrent)
    hide();
  else
    show();
}

/***********************************************/

PushButton::PushButton(QString text, QWidget *parent)
{
  setText(text);
  if(parent)
    setParent(parent);
}

/***********************************************/

void PushButton::setOrientation(Qt::Orientation orientation)
{
  _orientation = orientation;
}

/***********************************************/

Qt::Orientation PushButton::orientation() const
{
  return _orientation;
}

/***********************************************/

void PushButton::paintEvent(QPaintEvent *event)
{
  if(orientation() == Qt::Vertical)
  {
    QPainter painter(this);

    // support focus/hover styles
    QStyleOptionButton options;
    options.initFrom(this);
    options.icon = icon();
    options.iconSize = iconSize();
    if(isFlat())
        options.features |= QStyleOptionButton::Flat;
    if(menu())
        options.features |= QStyleOptionButton::HasMenu;
    if(autoDefault() || isDefault())
        options.features |= QStyleOptionButton::AutoDefaultButton;
    if(isDefault())
        options.features |= QStyleOptionButton::DefaultButton;
    if(isDown() || (menu() && menu()->isVisible()))
        options.state |= QStyle::State_Sunken;
    if(isChecked())
        options.state |= QStyle::State_On;
    if(!isFlat() && !isDown())
        options.state |= QStyle::State_Raised;

    style()->drawControl(QStyle::CE_PushButton, &options, &painter, this);

    // text render hint
    painter.setRenderHints(QPainter::TextAntialiasing);
    painter.rotate(-90);

    QRect r(0, 0, -height(), width());
    painter.drawText(r, Qt::AlignHCenter | Qt::AlignVCenter, QPushButton::text());
  }
  else  // default drawing for horizontal buttons
    QPushButton::paintEvent(event);
}

/***********************************************/

QSize PushButton::minimumSizeHint() const
{
  QSize size = QPushButton::minimumSizeHint();
  if(orientation() == Qt::Vertical)
    size.transpose();
  return size;
}

/***********************************************/

QSize PushButton::sizeHint() const
{
  QSize size = QPushButton::sizeHint();
  if(orientation() == Qt::Vertical)
    size.transpose();
  return size;
}

/***********************************************/

void OpenFilesTreeWidget::init(TabEnvironment *workspace, ActionList *actionList)
{
  if(workspace)
    this->workspace = workspace;
  else
    throw Exception("workspace uninitialized, cannot initialize Open Files side bar widget");

  this->actionList = *actionList;

  setColumnCount(1);
  header()->close();
  setSortingEnabled(true);
  sortByColumn(0, Qt::AscendingOrder);
  viewport()->setAcceptDrops(true); // internal & external drag

  populateTree();

  // select initial item in tree
  Tree *tree = workspace->currentTree();
  if(tree)
  {
    QTreeWidgetItemIterator it(this);
    while(*it)
    {
      if((*it)->data(0, Qt::UserRole+1).toString() == "file" && (*it)->toolTip(0) == tree->fileName())
      {
        setCurrentItem((*it));
        break;
      }
      ++it;
    }
  }

  connect(workspace, SIGNAL(fileChanged(const QString &, bool)), this, SLOT(fileChanged(const QString &, bool)));
  connect(this, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)), this, SLOT(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)));
  connect(this, SIGNAL(fileSelectionChanged(int)), workspace, SLOT(tabChanged(int)));
}

/***********************************************/

void OpenFilesTreeWidget::populateTree()
{
  // remember collapsed status of folders and changed status of files
  QStringList changedFiles;
  QStringList collapsedFolders;
  QTreeWidgetItemIterator it(this);
  while(*it)
  {
    if((*it)->data(0, Qt::UserRole+1).toString() == "folder" && !(*it)->isExpanded())
      collapsedFolders.push_back((*it)->toolTip(0));
    else if((*it)->data(0, Qt::UserRole+1).toString() == "file" && (*it)->data(0, Qt::UserRole+2).toBool() == true)
      changedFiles.push_back((*it)->toolTip(0));
    ++it;
  }

  // clear tree
  clear();

  // populate tree with opened files
  TreeWidgetItem *topLevelItem = nullptr;
  foreach(const QString &fileName, workspace->openedFileNames())
  {
    QStringList splitFileName = fileName.split("/");

    // add root folder as top level item if treeWidget doesn't already have it
    if(findItems(splitFileName[0], Qt::MatchFixedString).isEmpty())
    {
      topLevelItem = new TreeWidgetItem;
      topLevelItem->setText(0, splitFileName[0]);
      topLevelItem->setIcon(0, QIcon(":/icons/scalable/folder.svg"));
      addTopLevelItem(topLevelItem);
    }

    QTreeWidgetItem *parentItem = topLevelItem;

    // iterate through non-root directories (file name comes after)
    for(int i = 1; i < splitFileName.size() - 1; ++i)
    {
      // iterate through children of parentItem to see if this directory exists
      bool thisDirectoryExists = false;
      for(int j = 0; j < parentItem->childCount(); ++j)
      {
        if(splitFileName[i] == parentItem->child(j)->text(0))
        {
          thisDirectoryExists = true;
          parentItem = parentItem->child(j);
          break;
        }
      }

      if(!thisDirectoryExists)
      {
        parentItem = new TreeWidgetItem(parentItem);
        parentItem->setText(0, splitFileName[i].toUtf8());
        QStringList fileName;
        for(int j = 0; j <= i; j++)
          fileName << splitFileName[j];
        parentItem->setToolTip(0, fileName.join("/").toUtf8());
        parentItem->setData(0, Qt::UserRole+1, QVariant(QString("folder")));
        parentItem->setIcon(0, QIcon(":/icons/scalable/folder.svg"));
        parentItem->setFlags(Qt::ItemIsEnabled);
      }
    }

    QTreeWidgetItem *childItem = new TreeWidgetItem(parentItem);
    childItem->setText(0, splitFileName.last().toUtf8());
    childItem->setToolTip(0, fileName.toUtf8());
    childItem->setData(0, Qt::UserRole+1, QVariant(QString("file")));
    if(changedFiles.contains(fileName))
    {
      childItem->setIcon(0, QIcon(":/icons/scalable/document-save.svg"));
      childItem->setData(0, Qt::UserRole+2, QVariant(true));
    }
    else
    {
      childItem->setIcon(0, QIcon(":/icons/scalable/text-xml.svg"));
      childItem->setData(0, Qt::UserRole+2, QVariant(false));
    }
  }

  // reduce branches to a maximum of one folder above the first file
  it = QTreeWidgetItemIterator(this);
  while(*it)
  {
    if((*it)->data(0, Qt::UserRole+1).toString() == "folder" && (*it)->childCount() == 1 &&
       (*it)->child(0)->data(0, Qt::UserRole+1).toString() == "folder")
    {
      QTreeWidgetItem *parentParent = (*it)->parent();
      QTreeWidgetItem *child = (*it)->takeChild(0);
      if(parentParent)
        parentParent->addChild(child);
      else
        addTopLevelItem(child);
      it = QTreeWidgetItemIterator(this);
    }
    else
      ++it;
  }

  // delete branches that don't contain any files
  it = QTreeWidgetItemIterator(this);
  while(*it)
  {
    if((*it)->data(0, Qt::UserRole+1).toString() == "folder" && !hasFileChild((*it)))
    {
      delete (*it);
      it = QTreeWidgetItemIterator(this);
    }
    else
      ++it;
  }

  // if only one top level item is left and it has no files as direct children, remove it and make children top level items
  while(topLevelItemCount() == 1)
  {
    QTreeWidgetItem *topLevelItem = this->topLevelItem(0);
    int directFileChildren = 0;
    for(int i = 0; i < topLevelItem->childCount(); i++)
      if(topLevelItem->child(i)->data(0, Qt::UserRole+1).toString() == "file")
        directFileChildren++;

    if(directFileChildren == 0)
    {
      this->insertTopLevelItems(topLevelItemCount(), topLevelItem->takeChildren());
      delete topLevelItem;
    }
    else
      break;
  }

  // restore expanded/collapsed status of folders
  expandAll();
  it = QTreeWidgetItemIterator(this);
  while(*it)
  {
    if((*it)->data(0, Qt::UserRole+1).toString() == "folder" && collapsedFolders.contains((*it)->toolTip(0)))
      (*it)->setExpanded(false);
    ++it;
  }

  // add newly created files (without associated XML file) to top level
  foreach(const QString &fileName, workspace->newFileNames())
  {
    QTreeWidgetItem *item = new TreeWidgetItem();
    item->setText(0, fileName.toUtf8());
    item->setToolTip(0, fileName.toUtf8());
    item->setData(0, Qt::UserRole+1, QVariant(QString("file")));
    if(changedFiles.contains(fileName))
    {
      item->setIcon(0, QIcon(":/icons/scalable/document-save.svg"));
      item->setData(0, Qt::UserRole+2, QVariant(true));
    }
    else
    {
      item->setIcon(0, QIcon(":/icons/scalable/text-xml.svg"));
      item->setData(0, Qt::UserRole+2, QVariant(false));
    }
    this->insertTopLevelItem(topLevelItemCount(), item);
  }
}

/***********************************************/

bool OpenFilesTreeWidget::hasFileChild(const QTreeWidgetItem *item) const
{
  for(int i = 0; i < item->childCount(); i++)
    if(item->child(i)->data(0, Qt::UserRole+1).toString() == "file" || hasFileChild(item->child(i)))
      return true;

  return false;
}

/***********************************************/

void OpenFilesTreeWidget::currentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem */*previous*/)
{
  Tree *tree = workspace->currentTree();
  if(workspace && tree && current->data(0, Qt::UserRole+1) == "file" && current->toolTip(0) != tree->fileName())
  {
    disconnect(workspace, SIGNAL(fileChanged(const QString &, bool)), this, SLOT(fileChanged(const QString &, bool)));
    emit fileSelectionChanged(workspace->allFileNames().indexOf(current->toolTip(0)));
    connect(workspace, SIGNAL(fileChanged(const QString &, bool)), this, SLOT(fileChanged(const QString &, bool)));
  }
}

/***********************************************/

QTreeWidgetItem *OpenFilesTreeWidget::updateItemStatus(const QString &fileName, bool changed)
{
  QTreeWidgetItem *item = nullptr;
  // opened file with associated XML file
  if(!fileName.isEmpty())
  {
    QTreeWidgetItemIterator it(this);
    while(*it)
    {
      if((*it)->toolTip(0) == fileName && (*it)->data(0, Qt::UserRole+1) == "file")
      {
        item = (*it);
        break;
      }
      ++it;
    }
  }
  // newly created file without associated XML file
  else
  {
    for(int i = 0; i < topLevelItemCount(); i++)
      if(topLevelItem(i)->data(0, Qt::UserRole+1) == "file" && topLevelItem(i)->text(0) == workspace->currentTabText())
      {
        item = topLevelItem(i);
        break;
      }
  }

  if(item)
  {
    if(changed)
    {
      item->setIcon(0, QIcon(":/icons/scalable/document-save.svg"));
      item->setData(0, Qt::UserRole+2, QVariant(true));
    }
    else
    {
      item->setIcon(0, QIcon(":/icons/scalable/text-xml.svg"));
      item->setData(0, Qt::UserRole+2, QVariant(false));
    }
  }

  return item;
}

/***********************************************/

void OpenFilesTreeWidget::fileChanged(const QString &fileName, bool changed)
{
  this->blockSignals(true);
  populateTree();
  this->blockSignals(false);

  QTreeWidgetItem *item = updateItemStatus(fileName, changed);
  if(item)
    this->setCurrentItem(item);
}

/***********************************************/

void OpenFilesTreeWidget::mousePressEvent(QMouseEvent *event)
{
  QTreeWidgetItem *item = itemAt(event->pos());
  if(item && item->data(0, Qt::UserRole+1) == "folder")
    item->setExpanded(!item->isExpanded());
  else
    QTreeWidget::mousePressEvent(event);
}

/***********************************************/

void OpenFilesTreeWidget::contextMenuEvent(QContextMenuEvent *e)
{
  QTreeWidgetItem *item = itemAt(e->pos());
  if(!item || item->data(0, Qt::UserRole+1) == "folder")
    return;
  setCurrentItem(item);

  QMenu *contextMenu = new QMenu(this);
  contextMenu->addAction(actionList.fileNewAction);
  contextMenu->addAction(actionList.fileOpenAction);
  contextMenu->addAction(actionList.fileReOpenAction);
  contextMenu->addAction(actionList.fileShowInManagerAction);
  contextMenu->addSeparator();
  contextMenu->addAction(actionList.fileSaveAction);
  contextMenu->addAction(actionList.fileSaveAsAction);
  contextMenu->addSeparator();
  contextMenu->addAction(actionList.fileRunAction);
  contextMenu->addSeparator();
  contextMenu->addAction(actionList.fileCloseAction);
  contextMenu->addAction(actionList.fileCloseOtherAction);
  contextMenu->exec(e->globalPos());
  delete contextMenu;
}

/***********************************************/

void OpenFilesTreeWidget::dragEnterEvent(QDragEnterEvent *event)
{
  try
  {
    // file names?
    if(event->mimeData()->hasUrls())
    {
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

void OpenFilesTreeWidget::dragMoveEvent(QDragMoveEvent *event)
{
  try
  {
    // file names?
    if(event->mimeData()->hasUrls())
    {
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


void OpenFilesTreeWidget::dropEvent(QDropEvent *event)
{
  // file names?
  if(event->mimeData()->hasUrls())
  {
    QList<QUrl> urls = event->mimeData()->urls();
    for(int i=0; i<urls.size(); i++)
      workspace->openFile(urls.at(i).toLocalFile());
    event->acceptProposedAction();
  }
}

/***********************************************/

bool TreeWidgetItem::operator < (const QTreeWidgetItem &other) const
{
  int column = this->treeWidget()->sortColumn();

  if(this->data(column, Qt::UserRole+1) == "file" && other.data(column, Qt::UserRole+1) == "folder")
    return false;
  else if(this->data(column, Qt::UserRole+1) == "folder" && other.data(column, Qt::UserRole+1) == "file")
    return true;
  else
    return QTreeWidgetItem::operator<(other);
}

/***********************************************/
