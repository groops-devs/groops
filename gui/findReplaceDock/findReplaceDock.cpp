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


void FindReplaceDock::clickedNext()
{
  try
  {
    QRegExp regExp = getExpression();
    bool searchType = ui->checkBoxSearchType->isChecked();

    Tree *tree = mainWindow->getCurrentTree();
    if(!tree)
      return;
    TreeItem *treeItemStart = tree->getSelectedItem();
    if(!treeItemStart)
      treeItemStart = tree->rootElement()->item();

    QTreeWidgetItemIterator itemIter(treeItemStart, QTreeWidgetItemIterator::All);
    while(*itemIter)
    {
      if(searchType && *itemIter == treeItemStart)
        itemIter++; // to prevent lock

      TreeItem *treeItem = dynamic_cast<TreeItem*>(*itemIter);

      if(searchType)
      {
        if(treeItem && treeItem->treeElement())
          regExp.indexIn(treeItem->treeElement()->name(), 0);
        if(regExp.matchedLength() >= 0)
        {
          tree->setSelectedItem(treeItem);
          return;
        }
      }
      else // search value
      {
        int start, length;
        treeItem->selection(start, length);
        start += length;
        if(start < 0)
          start = 0;
        length = find(treeItem->treeElement(), regExp, false, start);
        if(length > 0)
        {
          tree->setSelectedItem(treeItem);
          treeItem->setSelection(start, length);
          ui->labelMessage->setText("");
          return;
        }
      }

      // Wrap search
      itemIter++;
      if(!*itemIter)
        itemIter = QTreeWidgetItemIterator(tree->rootElement()->item(), QTreeWidgetItemIterator::All);
      if(*itemIter == treeItemStart)
        break;
    } // while(*itemIter)

    messageNotFound();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void FindReplaceDock::clickedPrevious()
{
  try
  {
    QRegExp regExp = getExpression();
    bool searchType = ui->checkBoxSearchType->isChecked();

    Tree *tree = mainWindow->getCurrentTree();
    if(!tree)
      return;
    TreeItem *treeItemStart = tree->getSelectedItem();
    if(!treeItemStart)
      treeItemStart = tree->rootElement()->item();

    QTreeWidgetItemIterator itemIter(treeItemStart, QTreeWidgetItemIterator::All);
    while(*itemIter)
    {
      if(searchType && *itemIter == treeItemStart)
        itemIter--; // to prevent lock

      TreeItem *treeItem = dynamic_cast<TreeItem*>(*itemIter);

      if(searchType)
      {
        if(treeItem && treeItem->treeElement())
          regExp.indexIn(treeItem->treeElement()->name(), 0);
        if(regExp.matchedLength() >= 0)
        {
          tree->setSelectedItem(treeItem);
          return;
        }
      }
      else // search value
      {
        int start, length;
        treeItem->selection(start, length);
        if(start < 0)
          start = treeItem->treeElement()->selectedValue().size();
        if(--start >= 0)
        {
          length = find(treeItem->treeElement(), regExp, true, start);
          if(length > 0)
          {
            tree->setSelectedItem(treeItem);
            treeItem->setSelection(start, length);
            ui->labelMessage->setText("");
            return;
          }
        }
      }


      // Wrap search
      itemIter--;
      if(!*itemIter)
      {
        // find last item
       QTreeWidgetItemIterator itemIterTmp(treeItemStart, QTreeWidgetItemIterator::All);
       while(*itemIterTmp)
         itemIter = itemIterTmp++;
      }
      if(*itemIter == treeItemStart)
        break;
    } // while(*itemIter)

    messageNotFound();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
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
    if(start >= 0 && length > 0 && getExpression().exactMatch(text.mid(start, length)))
    {
      text.replace(start, length, ui->editReplace->text());
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
    QRegExp regExp     = getExpression();
    QString replaceStr = ui->editReplace->text();

    Tree *tree = mainWindow->getCurrentTree();
    if(!tree)
      return;

    int count = 0;
    TreeItem *treeItem = tree->rootElement()->item();
    QTreeWidgetItemIterator itemIter(treeItem, QTreeWidgetItemIterator::All);
    while(*itemIter)
    {
      TreeElement *treeElement = dynamic_cast<TreeItem*>(*itemIter)->treeElement();
      if(treeElement->isEditable() && replace(treeElement, regExp, replaceStr))
        count++;
      itemIter++;
    } // while(*itemIter)

    if(count)
    {
      ui->labelMessage->setText(QString("%1 matches replaced.").arg(count));
      return;
    }

    messageNotFound();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

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

void FindReplaceDock::messageNotFound()
{
  ui->labelMessage->setText(QString("No matches found."));
  QApplication::beep();
}

/***********************************************/

QRegExp FindReplaceDock::getExpression()
{
  try
  {
    QRegExp regExp(ui->editFind->text());
    regExp.setCaseSensitivity(ui->checkBoxMatchCase->checkState() ? Qt::CaseSensitive : Qt::CaseInsensitive);
    if(!ui->checkBoxRegExp->checkState())
      regExp.setPatternSyntax(QRegExp::FixedString);
    //regExp.setPatternSyntax(QRegExp::RegExp2);
// QRegExp::RegExp       0 A rich Perl-like pattern matching syntax. This is the default.
// QRegExp::RegExp2      3 Like RegExp, but with greedy quantifiers. This will be the default in Qt 5. (Introduced in Qt 4.2.)
// QRegExp::Wildcard     1 This provides a simple pattern matching syntax similar to that used by shells (command interpreters) for "file globbing". See Wildcard Matching.
// QRegExp::WildcardUnix 4 This is similar to Wildcard but with the behavior of a Unix shell. The wildcard characters can be escaped with the character "\".
// QRegExp::FixedString  2 The pattern is a fixed string. This is equivalent to using the RegExp pattern on a string in which all metacharacters are escaped using escape().
// QRegExp::W3CXmlSchema11
    return regExp;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

int FindReplaceDock::find(TreeElement *treeElement, const QRegExp &regExp, bool backwards, int &start)
{
  try
  {
    if(!treeElement || !treeElement->isEditable() || treeElement->isLinked())
      return -1;
    if(backwards)
      start = regExp.lastIndexIn(treeElement->selectedValue(), start);
    else
      start = regExp.indexIn(treeElement->selectedValue(), start);
    return regExp.matchedLength();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool FindReplaceDock::replace(TreeElement *treeElement, const QRegExp &regExp, const QString &replaceStr)
{
  try
  {
    QString text = treeElement->selectedValue();
    text.replace(regExp, replaceStr);
    if(text == treeElement->selectedValue()) // no change?
      return false;
    treeElement->changeSelectedValue(text);
    return true;
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
