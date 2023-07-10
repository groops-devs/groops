/***********************************************/
/**
* @file findReplaceDock.cpp
*
* @brief Find/replace dock.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2013-10-30
*/
/***********************************************/

#include "ui_findReplaceDock.h"
#include <QtWidgets>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeItem.h"
#include "mainWindow/mainWindow.h"
#include "findReplaceDock.h"

/***********************************************/

FindReplaceDock::FindReplaceDock(MainWindow *parent) : QDockWidget(parent), ui(new Ui::FindReplaceDock)
{
  try
  {
    mainWindow = parent;
    ui->setupUi(this);
    ui->labelMessage->setText("");

    this->settings = new QSettings(this);
    resize(settings->value("findReplaceDock/size", size()).toSize());

    connect(ui->editFind, SIGNAL(textEdited(const QString &)), this, SLOT(editFindTextChanged(const QString &)));
    connect(ui->checkBoxSearchType, &QCheckBox::stateChanged, this, &FindReplaceDock::searchTypeChanged);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

FindReplaceDock::~FindReplaceDock()
{
  settings->setValue("findReplaceDock/size", size());

  delete ui;
}

/***********************************************/

void FindReplaceDock::setFocusOnFind()
{
  try
  {
    ui->editFind->selectAll();
    ui->editFind->setFocus();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void FindReplaceDock::setFindText(QString text)
{
  try
  {
    ui->editFind->setText(text);
    editFindTextChanged(text);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QString FindReplaceDock::getFindText() const
{
  try
  {
    return ui->editFind->text();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}


/***********************************************/
/***********************************************/

void FindReplaceDock::editFindTextChanged(const QString &text)
{
  try
  {
    bool enable = !text.isEmpty();
    bool searchType = ui->checkBoxSearchType->isChecked();
    ui->buttonNext->setEnabled(enable);
    ui->buttonPrevious->setEnabled(enable);
    ui->buttonReplace->setEnabled(enable && !searchType);
    ui->buttonReplaceAll->setEnabled(enable && !searchType);
    ui->labelMessage->setText("");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void FindReplaceDock::searchTypeChanged()
{
  try
  {
    bool enable = !ui->editFind->text().isEmpty();
    bool searchType = ui->checkBoxSearchType->isChecked();
    ui->editReplace->setEnabled(!searchType);
    ui->buttonReplace->setEnabled(enable && !searchType);
    ui->buttonReplaceAll->setEnabled(enable && !searchType);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}


/***********************************************/
/***********************************************/

void FindReplaceDock::clickedNext()
{
  findNext(TRUE/*forwards*/);
}

/***********************************************/

void FindReplaceDock::clickedPrevious()
{
  findNext(FALSE/*forwards*/);
}

/***********************************************/

void FindReplaceDock::clickedReplace()
{
  try
  {
    Tree *tree = mainWindow->getCurrentTree();
    if(!tree)
      return;
    TreeItem *treeItem = tree->getSelectedItem();
    int start, length;
    treeItem->selection(start, length);
    QString text = treeItem->treeElement()->selectedValue();
    if((start >= 0) && (length > 0) && (getExpression().match(text.mid(start, length)).capturedLength() == length))
    {
      text.replace(getExpression(), ui->editReplace->text());
      //text.replace(start, length, ui->editReplace->text());
      if(text == treeItem->treeElement()->selectedValue()) // no change?
        return;
      treeItem->treeElement()->changeSelectedValue(text);
      treeItem->setSelection(start, 0);
    }
    clickedNext();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void FindReplaceDock::clickedReplaceAll()
{
  try
  {
    QRegularExpression regExp = getExpression();
    QString replaceStr = ui->editReplace->text();

    Tree *tree = mainWindow->getCurrentTree();
    if(!tree)
      return;

    int count = 0;
    for(QTreeWidgetItemIterator itemIter(tree->rootElement()->item(), QTreeWidgetItemIterator::All); *itemIter; itemIter++)
    {
      TreeElement *treeElement = dynamic_cast<TreeItem*>(*itemIter)->treeElement();
      if(treeElement->isEditable())
      {
        QString text = treeElement->selectedValue();
        text.replace(regExp, replaceStr);
        if(text != treeElement->selectedValue()) // change?
        {
          treeElement->changeSelectedValue(text);
          count++;
        }
      }
    }

    if(count)
      ui->labelMessage->setText(QString("%1 matches replaced.").arg(count));
    else
      messageNotFound();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void FindReplaceDock::findNext(Bool forwards)
{
  try
  {
    QRegularExpression regExp = getExpression();
    bool searchType = ui->checkBoxSearchType->isChecked();

    Tree *tree = mainWindow->getCurrentTree();
    if(!tree)
      return;
    TreeItem *treeItemStart = tree->getSelectedItem();
    if(!treeItemStart)
      treeItemStart = tree->rootElement()->item();

    QTreeWidgetItemIterator itemIter(treeItemStart, QTreeWidgetItemIterator::All);
    if(searchType)
      (forwards) ? itemIter++ : itemIter--; // start with next item to prevent lock

    for(;;)
    {
      if(!*itemIter) // last element? -> warp search
      {
        if(forwards)
          itemIter = QTreeWidgetItemIterator(tree->rootElement()->item(), QTreeWidgetItemIterator::All);
        else // find last item
          for(QTreeWidgetItemIterator i(treeItemStart, QTreeWidgetItemIterator::All); *i; i++)
            itemIter = i;
      }

      TreeItem *treeItem = dynamic_cast<TreeItem*>(*itemIter);
      if(treeItem && treeItem->treeElement())
      {
        auto treeElement = treeItem->treeElement();
        if(searchType)
        {
          if(regExp.match(treeElement->name(), 0).hasMatch())
          {
            tree->setSelectedItem(treeItem);
            return;
          }
        } // search value
        else if(treeElement->isEditable() && !treeElement->isLinked() && treeElement->selectedValue().size())
        {
          int start, length, startMatch = -1;
          treeItem->selection(start, length);
          QRegularExpressionMatch match;
          if(forwards)
            startMatch = treeElement->selectedValue().indexOf(regExp, std::max(start, 0)+length, &match);
          else if(start != 0)
            startMatch = treeElement->selectedValue().lastIndexOf(regExp, ((start>1) ? start : treeElement->selectedValue().size())-1, &match);
          if(startMatch >= 0)
          {
            tree->setSelectedItem(treeItem);
            treeItem->setSelection(startMatch, match.capturedLength());
            ui->labelMessage->setText("");
            return;
          }
        } // if(searchValue)
      } // if(treeItem)

      // next item
      (forwards) ? itemIter++ : itemIter--;
      if(*itemIter == treeItemStart)
        break;
    } // for(;;)

    messageNotFound();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void FindReplaceDock::messageNotFound()
{
  ui->labelMessage->setText(QString("No matches found."));
  QApplication::beep();
}

/***********************************************/

QRegularExpression FindReplaceDock::getExpression()
{
  try
  {
    QString pattern = ui->editFind->text();
    if(!ui->checkBoxRegExp->checkState())
      pattern.replace(QRegularExpression(R"(([\.\*\?\+\/\^\$\|\(\)\[\]\{\}\\]))"), R"(\\1)");
    QRegularExpression regExp(pattern);
    if(ui->checkBoxMatchCase->checkState())
      regExp.setPatternOptions(QRegularExpression::CaseInsensitiveOption);
    return regExp;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void FindReplaceDock::closeEvent(QCloseEvent *e)
{
  try
  {
    hide();
    e->ignore();
    //e->accept();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void FindReplaceDock::changeEvent(QEvent *e)
{
  QDockWidget::changeEvent(e);
  switch(e->type())
  {
    case QEvent::LanguageChange:
      ui->retranslateUi(this);
      break;
    default:
      break;
  }
}

/***********************************************/
