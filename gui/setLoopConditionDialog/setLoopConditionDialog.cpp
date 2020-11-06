/***********************************************/
/**
* @file setLoopConditionDialog.cpp
*
* @brief Set loop/condition for an unbounded element.
*
* @author Sebastian Strasser
* @date 2017-02-09
*/
/***********************************************/

#include "setLoopConditionDialog.h"
#include "ui_setLoopConditionDialog.h"
#include <QPushButton>
#include <QScreen>

/***********************************************/

SetLoopConditionDialog::SetLoopConditionDialog(TreeElementGlobal *globalRoot, const QString &type, QWidget *parent)
  : QDialog(parent), ui(new Ui::SetLoopConditionDialog)
{
  try
  {
    ui->setupUi(this);

    QSize screenSize = QGuiApplication::primaryScreen()->size();
    setMinimumSize(screenSize.width()/5, screenSize.height()/10);
    resize(minimumSizeHint());

    if(type != "loop" && type != "condition")
      throw(Exception("invalid type: "+type.toStdString()));
    setWindowTitle("Set " + type + " for element - GROOPS");

    // only allow variable names with letters and numbers
    QRegExpValidator *validator = new QRegExpValidator(QRegExp("[a-zA-Z]([a-zA-Z0-9])*"), this);
    ui->lineEditName->setValidator(validator);

    // fill combo box with existing loops/conditions
    std::vector<TreeElement> a;
    for(int i = 0; i < globalRoot->childrenCount(); i++)
      if(globalRoot->childAt(i)->type() == type+"Type" && !globalRoot->childAt(i)->label().isEmpty())
        ui->comboBox->addItem(globalRoot->childAt(i)->name());
    if(ui->comboBox->count() > 0)
    {
      ui->radioButtonExisting->setChecked(true);
      ui->lineEditName->setEnabled(false);
    }
    else
    {
      ui->radioButtonExisting->setEnabled(false);
      ui->radioButtonNew->setChecked(true);
      ui->comboBox->setEnabled(false);
      ui->comboBox->addItem("<no "+type+"s exist in global>");
      ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);
    }

    // list of existing global variable names
    globalElementNames = globalRoot->getChildrenNames();

    // signal-slot connections
    connect(ui->radioButtonExisting, SIGNAL(toggled(bool)),               this, SLOT(radioButtonToggled(bool)));
    connect(ui->lineEditName,        SIGNAL(textEdited(const QString &)), this, SLOT(nameEdited(const QString &)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

SetLoopConditionDialog::~SetLoopConditionDialog()
{
  delete ui;
}

/***********************************************/

void SetLoopConditionDialog::radioButtonToggled(bool checked)
{
  try
  {
    ui->comboBox->setEnabled(checked);
    ui->lineEditName->setEnabled(!checked);
    if(checked) // radioButtonExisting checked
      ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);
    else        // radioButtonNew checked
      ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(!ui->lineEditName->text().isEmpty() && !globalElementNames.contains(ui->lineEditName->text(), Qt::CaseInsensitive));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void SetLoopConditionDialog::nameEdited(const QString &name)
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

QString SetLoopConditionDialog::name()
{
  try
  {
    return ui->radioButtonExisting->isChecked() ? ui->comboBox->currentText() : ui->lineEditName->text();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
