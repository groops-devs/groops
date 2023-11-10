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

SideBar::SideBar(QWidget *parent) : QWidget(parent)
{
  try
  {
    layout = new QVBoxLayout(this);
    layout->setAlignment(this, Qt::AlignTop);
    layout->setSpacing(0);
    layout->setContentsMargins(0,0,0,0);

    _stackedWidget = new QStackedWidget(this);
    _stackedWidget->setHidden(settings.value("sideBar/isHidden", false).toBool());
    _lastWidth = settings.value("sideBar/width", parentWidget()->width()/5).toInt();

    buttonGroup.setExclusive(false);
    connect(&buttonGroup, &QButtonGroup::idClicked, this, &SideBar::buttonClicked);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

SideBar::~SideBar()
{
  settings.setValue("sideBar/width",    lastWidth());
  settings.setValue("sideBar/isHidden", stackedWidget()->isHidden());
}

/***********************************************/

void SideBar::addSideBarWidget(const QString &buttonLabel, QWidget *widget)
{
  // add push button
  PushButtonVertical *button = new PushButtonVertical(buttonLabel, this);
  button->setMaximumWidth(button->minimumSizeHint().height()); // orientation = Qt::Vertical;
  button->setFlat(true);
  button->setCheckable(true);
  button->setChecked(false);
  buttonGroup.addButton(button, _stackedWidget->count());
  layout->addWidget(button);

  _stackedWidget->addWidget(widget);

  buttonGroup.button(0)->setChecked(!_stackedWidget->isHidden());
  _stackedWidget->setCurrentIndex(0);
}

/***********************************************/

void SideBar::buttonClicked(int id)
{
  // deselect all others
  for(int i=0; i<_stackedWidget->count(); i++)
    if(i != id)
      buttonGroup.button(i)->setChecked(false);

  _stackedWidget->setVisible(buttonGroup.button(id)->isChecked());
  if(buttonGroup.button(id)->isChecked())
    _stackedWidget->setCurrentIndex(id);
}

/***********************************************/

void SideBar::widthChanged(int width, int /*index*/)
{
  if(!_stackedWidget->isHidden())
    _lastWidth = width;
}

/***********************************************/

void PushButtonVertical::paintEvent(QPaintEvent */*event*/)
{
  // support focus/hover styles
  QStyleOptionButton options;
  options.initFrom(this);
  options.icon     = icon();
  options.iconSize = iconSize();
  if(isFlat())                                    options.features |= QStyleOptionButton::Flat;
  if(menu())                                      options.features |= QStyleOptionButton::HasMenu;
  if(autoDefault() || isDefault())                options.features |= QStyleOptionButton::AutoDefaultButton;
  if(isDefault())                                 options.features |= QStyleOptionButton::DefaultButton;
  if(isDown() || (menu() && menu()->isVisible())) options.state    |= QStyle::State_Sunken;
  if(isChecked())                                 options.state    |= QStyle::State_On;
  if(!isFlat() && !isDown())                      options.state    |= QStyle::State_Raised;

  QPainter painter(this);
  style()->drawControl(QStyle::CE_PushButton, &options, &painter, this);

  // text render hint
  painter.setRenderHints(QPainter::TextAntialiasing);
  painter.rotate(-90);
  painter.drawText(QRect(0, 0, -height(), width()), Qt::AlignHCenter | Qt::AlignVCenter, QPushButton::text());
}

/***********************************************/

OpenFilesTreeWidget::OpenFilesTreeWidget(TabEnvironment *tabEnvironment, ActionList &actionList) : actionList(actionList), tabEnvironment(tabEnvironment)
{
  setColumnCount(1);
  header()->close();
  setSortingEnabled(true);
  sortByColumn(0, Qt::AscendingOrder);
  viewport()->setAcceptDrops(true); // internal & external drag
  populateTree();

  connect(tabEnvironment, SIGNAL(fileChanged(const QString&, const QString&, bool)), this, SLOT(fileChanged(const QString&, const QString&, bool)));
  connect(this, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)), this, SLOT(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)));
  connect(this, SIGNAL(fileSelectionChanged(int)), tabEnvironment, SLOT(tabChanged(int)));
}

/***********************************************/

void OpenFilesTreeWidget::populateTree()
{
  // remember collapsed status of folders
  QStringList collapsedFolders;
  for(auto it = QTreeWidgetItemIterator(this); *it; it++)
    if(((*it)->data(0, Qt::UserRole).toString() == "folder") && !(*it)->isExpanded())
      collapsedFolders.push_back((*it)->toolTip(0));

  clear(); // clear tree

  // populate tree with opened files
  TreeWidgetItem *topLevelItem = nullptr;
  for(auto &tree : tabEnvironment->trees())
    if(!tree->fileName().isEmpty())
    {
      QStringList splitFileName = tree->fileName().split("/");

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
      for(int i=1; i<splitFileName.size()-1; i++)
      {
        // iterate through children of parentItem to see if this directory exists
        bool thisDirectoryExists = false;
        for(int k=0; k<parentItem->childCount(); k++)
          if(splitFileName[i] == parentItem->child(k)->text(0))
          {
            thisDirectoryExists = true;
            parentItem = parentItem->child(k);
            break;
          }

        if(!thisDirectoryExists)
        {
          parentItem = new TreeWidgetItem(parentItem);
          parentItem->setText(0, splitFileName[i].toUtf8());
          QStringList fileName;
          for(int k=0; k<=i; k++)
            fileName<<splitFileName[k];
          parentItem->setToolTip(0, fileName.join("/").toUtf8());
          parentItem->setData(0, Qt::UserRole, QVariant(QString("folder")));
          parentItem->setIcon(0, QIcon(":/icons/scalable/folder.svg"));
          parentItem->setFlags(Qt::ItemIsEnabled);
        }
      }

      QTreeWidgetItem *childItem = new TreeWidgetItem(parentItem);
      childItem->setText(0,    tree->caption());
      childItem->setToolTip(0, tree->fileName());
      childItem->setData(0, Qt::UserRole, QVariant(QString("file")));
      childItem->setIcon(0, tree->isClean() ? QIcon(":/icons/scalable/text-xml.svg") : QIcon(":/icons/scalable/document-save.svg"));
    }

  // reduce branches to a maximum of one folder above the first file
  auto it = QTreeWidgetItemIterator(this);
  while(*it)
  {
    if(((*it)->data(0, Qt::UserRole).toString() == "folder") && ((*it)->childCount() == 1) &&
       ((*it)->child(0)->data(0, Qt::UserRole).toString() == "folder"))
    {
      QTreeWidgetItem *child = (*it)->takeChild(0);
      if((*it)->parent())
        (*it)->parent()->addChild(child);
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
    if((*it)->data(0, Qt::UserRole).toString() == "folder" && !hasFileChild((*it)))
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
    for(int i=0; i<topLevelItem->childCount(); i++)
      if(topLevelItem->child(i)->data(0, Qt::UserRole).toString() == "file")
        directFileChildren++;
    if(directFileChildren > 0)
      break;
    this->insertTopLevelItems(topLevelItemCount(), topLevelItem->takeChildren());
    delete topLevelItem;
  }

  // restore expanded/collapsed status of folders
  expandAll();
  for(auto it=QTreeWidgetItemIterator(this); *it; it++)
    if(((*it)->data(0, Qt::UserRole).toString() == "folder") && collapsedFolders.contains((*it)->toolTip(0)))
      (*it)->setExpanded(false);

  // add newly created files (without associated XML file) to top level
  for(auto &tree : tabEnvironment->trees())
    if(tree->fileName().isEmpty())
    {
      QTreeWidgetItem *item = new TreeWidgetItem();
      item->setText(0, tree->caption());
      // item->setToolTip(0, tree->caption());
      item->setData(0, Qt::UserRole, QVariant(QString("file")));
      item->setIcon(0, tree->isClean() ? QIcon(":/icons/scalable/text-xml.svg") : QIcon(":/icons/scalable/document-save.svg"));
      insertTopLevelItem(topLevelItemCount(), item);
    }
}

/***********************************************/

bool OpenFilesTreeWidget::hasFileChild(const QTreeWidgetItem *item) const
{
  for(int i=0; i<item->childCount(); i++)
    if((item->child(i)->data(0, Qt::UserRole).toString() == "file") || hasFileChild(item->child(i)))
      return true;
  return false;
}

/***********************************************/

void OpenFilesTreeWidget::currentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem */*previous*/)
{
  for(auto &tree : tabEnvironment->trees())
    if((tree != tabEnvironment->currentTree()) && (current->text(0) == tree->caption()) && (current->toolTip(0) == tree->fileName()))
    {
      disconnect(tabEnvironment, SIGNAL(fileChanged(const QString&, const QString&, bool)), this, SLOT(fileChanged(const QString&, const QString&, bool)));
      tabEnvironment->setCurrentTree(tree);
      connect(tabEnvironment, SIGNAL(fileChanged(const QString&, const QString&, bool)), this, SLOT(fileChanged(const QString&, const QString&, bool)));
      break;
    }
}

/***********************************************/

void OpenFilesTreeWidget::fileChanged(const QString &caption, const QString &fileName, bool /*isClean*/)
{
  this->blockSignals(true);
  populateTree();
  for(auto it=QTreeWidgetItemIterator(this); *it; it++)
    if(((*it)->text(0) == caption) && ((*it)->toolTip(0) == fileName))
      this->setCurrentItem(*it);
  this->blockSignals(false);
}

/***********************************************/

void OpenFilesTreeWidget::mousePressEvent(QMouseEvent *event)
{
  QTreeWidgetItem *item = itemAt(event->pos());
  if(item && item->data(0, Qt::UserRole).toString() == "folder")
    item->setExpanded(!item->isExpanded());
  else
    QTreeWidget::mousePressEvent(event);
}

/***********************************************/

void OpenFilesTreeWidget::contextMenuEvent(QContextMenuEvent *e)
{
  QTreeWidgetItem *item = itemAt(e->pos());
  if(!item || item->data(0, Qt::UserRole).toString() == "folder")
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
  if(event->mimeData()->hasUrls()) // file names?
    event->acceptProposedAction();
  else
    event->ignore();
}

/***********************************************/

void OpenFilesTreeWidget::dragMoveEvent(QDragMoveEvent *event)
{
  if(event->mimeData()->hasUrls()) // file names?
    event->acceptProposedAction();
  else
    event->ignore();
}

/***********************************************/

void OpenFilesTreeWidget::dropEvent(QDropEvent *event)
{
  if(event->mimeData()->hasUrls())  // file names?
  {
    for(auto &url : event->mimeData()->urls())
      tabEnvironment->fileOpen(url.toLocalFile());
    event->acceptProposedAction();
  }
}

/***********************************************/

bool TreeWidgetItem::operator<(const QTreeWidgetItem &other) const
{
  if((this->data(0, Qt::UserRole).toString() == "file") && (other.data(0, Qt::UserRole).toString() == "folder"))
    return false;
  else if((this->data(0, Qt::UserRole).toString() == "folder") && (other.data(0, Qt::UserRole).toString() == "file"))
    return true;
  else
    return QTreeWidgetItem::operator<(other);
}

/***********************************************/
