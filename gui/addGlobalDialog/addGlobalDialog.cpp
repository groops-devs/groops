/***********************************************/
/**
* @file addGlobalDialog.cpp
*
* @brief Dialog to add global variable.
*
* @author Sebastian Strasser
* @date 2016-10-17
*/
/***********************************************/

#include <QPushButton>
#include <QScreen>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementGlobal.h"
#include "addGlobalDialog.h"
#include "ui_addGlobalDialog.h"

/***********************************************/

AddGlobalDialog::AddGlobalDialog(TreeElementGlobal *globalRoot, QWidget *parent)
  : QDialog(parent), ui(new Ui::AddGlobalDialog)
{
  try
  {
    ui->setupUi(this);

    QSize screenSize = QGuiApplication::primaryScreen()->size();
    setMinimumSize(screenSize.width()/5, screenSize.height()/10);
    resize(minimumSizeHint());

    // only allow variable names with letters and numbers
    QRegExpValidator *validator = new QRegExpValidator(QRegExp("[a-zA-Z]([a-zA-Z0-9])*"), this);
    ui->lineEditName->setValidator(validator);

    // fill combo box with types
    xsdElements = globalRoot->xsdComplexGlobalTypes()->element;
    for(auto element : xsdElements)
      ui->comboBoxType->addItem(element->type);
    if(ui->comboBoxType->count() < 1)
      throw(Exception("no global element types found in schema"));

    // disable Ok button if line edit is empty
    ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);

    // list of existing global variable names
    globalElementNames = globalRoot->getChildrenNames();

    // signal-slot connections
    connect(ui->lineEditName, SIGNAL(textEdited(const QString &)), this, SLOT(nameEdited(const QString &)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

AddGlobalDialog::~AddGlobalDialog()
{
  delete ui;
}

/***********************************************/

void AddGlobalDialog::nameEdited(const QString &name)
{
  try
  {
    ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(!name.isEmpty() && !globalElementNames.contains(name, Qt::CaseInsensitive));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QString AddGlobalDialog::elementName()
{
  try
  {
    return ui->lineEditName->text();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XsdElementPtr AddGlobalDialog::elementType()
{
  try
  {
    return xsdElements.at(static_cast<size_t>(ui->comboBoxType->currentIndex()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

