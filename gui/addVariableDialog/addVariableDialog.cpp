/***********************************************/
/**
* @file addVariableDialog.cpp
*
* @brief Dialog to add variables.
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
#include "addVariableDialog.h"
#include "ui_addVariableDialog.h"

/***********************************************/

AddVariableDialog::AddVariableDialog(Tree *tree, const QString &type, const QString &label, bool disablePlace, bool disableCreateLink, QWidget *parent)
  : QDialog(parent), sourceType(type), disableCreateLink(disableCreateLink), ui(new Ui::AddVariableDialog)
{
  try
  {
    ui->setupUi(this);

    QSize screenSize = QGuiApplication::primaryScreen()->size();
    setMinimumSize(screenSize.width()/5, screenSize.height()/10);
    resize(minimumSizeHint());

    // fill combo box with types
    for(auto xsdElement : tree->xsdElements())
      ui->comboBoxType->addItem(xsdElement->type);
    if(ui->comboBoxType->count() < 1)
      throw(Exception("no global element types found in schema"));
    ui->comboBoxType->setCurrentText(type);

    bool fromTemplate = !label.isEmpty() && (ui->comboBoxType->currentText() == type);

    if(fromTemplate)
      ui->lineEditName->setText(label);
    // only allow variable names with letters and numbers
    QRegularExpressionValidator *validator = new QRegularExpressionValidator(QRegularExpression("[a-zA-Z]([a-zA-Z0-9])*"), this);
    ui->lineEditName->setValidator(validator);


    ui->radioButtonGlobal->setEnabled(!disablePlace);
    ui->radioButtonLocal->setEnabled(!disablePlace);

    ui->checkBoxCreateLink->setEnabled(!disableCreateLink && fromTemplate);
    ui->checkBoxCreateLink->setChecked(!disableCreateLink && fromTemplate);

    // disable Ok button if line edit is empty
    ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(!this->label().isEmpty());

    ui->lineEditName->setFocus();

    // signal-slot connections
    connect(ui->lineEditName, SIGNAL(textEdited(const QString &)), this, SLOT(nameEdited(const QString &)));
    connect(ui->comboBoxType, SIGNAL(currentIndexChanged(int)),    this, SLOT(typeChanged(int)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

AddVariableDialog::~AddVariableDialog()
{
  delete ui;
}

/***********************************************/

void AddVariableDialog::nameEdited(const QString &label)
{
  try
  {
    ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(!label.isEmpty());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}


/***********************************************/

void AddVariableDialog::typeChanged(int /*index*/)
{
  try
  {
    if(ui->comboBoxType->currentText() == sourceType)
      ui->checkBoxCreateLink->setEnabled(!disableCreateLink);
    else
    {
      ui->checkBoxCreateLink->setEnabled(false);
      ui->checkBoxCreateLink->setChecked(false);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}


/***********************************************/

QString AddVariableDialog::label() const
{
  return ui->lineEditName->text();
}

/***********************************************/

QString AddVariableDialog::type() const
{
  return ui->comboBoxType->currentText();
}

/***********************************************/

Bool AddVariableDialog::inGlobal() const
{
  return ui->radioButtonGlobal->isEnabled() && ui->radioButtonGlobal->isChecked();
}

/***********************************************/

Bool AddVariableDialog::setLink() const
{
  return ui->checkBoxCreateLink->isChecked();
}

/***********************************************/

