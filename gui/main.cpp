/***********************************************/
/**
* @file main.cpp
*
* @brief Main.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2011-07-12
*/
/***********************************************/

#include <QtDebug>
#include <QtWidgets>
#include <QApplication>
#include <QMessageBox>
#include <QFileInfo>
#include <QMainWindow>
#include <QSettings>
#include "base/importGroops.h"
#include "base/schema.h"
#include "settingsDialog/settingsPathDialog.h"
#include "mainWindow/mainWindow.h"

/***********************************************/

class MyApplication : public QApplication
{
public:
  MyApplication(int &argc, char **argv) : QApplication(argc, argv) {}

  virtual bool notify(QObject *receiver, QEvent *event);
};

/***********************************************/

bool MyApplication::notify(QObject *receiver, QEvent *event)
{
  try
  {
    return QApplication::notify(receiver, event);
  }
  catch(std::exception &e)
  {
    qWarning("%s", e.what());
    QMessageBox::critical(nullptr, "GROOPS", QString("****** Error ******\n")+e.what());
  }
  catch(...)
  {
    qWarning("Unknown ERROR");
    QMessageBox::critical(nullptr, "GROOPS", "Unknown ERROR");
  }
  return false;
}

/***********************************************/

int main(int argc, char *argv[])
{
  try
  {
    MyApplication app(argc, argv);
    app.connect( &app, SIGNAL(lastWindowClosed()), &app, SLOT(quit()) );

    QCoreApplication::setApplicationName("groopsGui");
    QCoreApplication::setOrganizationName("groops-devs");
    QApplication::setFont(QApplication::font("QMenu"));

    QCommandLineParser parser;
    std::stringstream ss;
    ss<<"Graphical User Interface for GROOPS (Gravity Recovery Object Oriented Programming System)"<<std::endl<<std::endl;
    ss<<"Written by: Torsten Mayer-Guerr, Wolfgang Mayer-Guerr, Sebastian Strasser"<<std::endl;
    ss<<"GitHub repository: https://github.com/groops-devs/groops"<<std::endl;
    ss<<"(Compiled: "<<__DATE__<<" "<<__TIME__<<")";
    parser.setApplicationDescription(QString::fromStdString(ss.str()));
    parser.addHelpOption();
    parser.addOptions(
    {
      {{"c", "clean"},  "Do not open files from last session at startup."},
      {{"f", "file"},   "Open <file> at startup.", "file"},
      {{"s", "schema"}, "Set XML schema <file> at startup.", "file"}
    });
    parser.process(app); // Handle command line options

    MainWindow window;
    window.show();

    QSettings settings;
    QString fileNameSchema = settings.value("files/schemaFile").toString();

    // Set XML schema if passed as argument
    // ------------------------------------
    if(parser.isSet("schema"))
    {
      QFileInfo fileFromArg(parser.value("schema"));
      Schema schema;
      if(fileFromArg.isFile() && schema.readFile(fileFromArg.absoluteFilePath()))
      {
        fileNameSchema = fileFromArg.absoluteFilePath();
        settings.setValue("files/schemaFile", fileNameSchema);
        QStringList schemaFiles = settings.value("files/schemaFiles").toStringList();
        if(!schemaFiles.contains(fileNameSchema))
        {
          schemaFiles.append(fileNameSchema);
          schemaFiles.sort();
          settings.setValue("files/schemaFiles", schemaFiles);
        }
      }
      else
        QMessageBox::critical(&window, "GROOPS", QString("File '%1' passed as argument seems not to be a valid XSD schema").arg(fileFromArg.absoluteFilePath()));
    }

    // check schema
    // ------------
    if(fileNameSchema.isEmpty())
    {
      QMessageBox::information(&window , "GROOPS", "GROOPS seems not to be configured yet. You should set at least the XSD schema file.");
      SettingsPathDialog dialog(&window);
      if(!dialog.exec())
      {
        window.close();
        return 0;
      }
    }

    // Open the files from the previous session (from settings)
    // --------------------------------------------------------
    bool isOpen = false;
    if(!parser.isSet("clean"))
      for(const auto &fileName : settings.value("openFiles").toStringList())
        if(QFileInfo::exists(fileName))
          isOpen = window.fileOpen(QFileInfo(fileName).absoluteFilePath()) || isOpen;

    // Open files from command line arguments
    // --------------------------------------
    for(const auto &fileName : parser.values("file"))
      if(QFileInfo::exists(fileName))
        isOpen = window.fileOpen(QFileInfo(fileName).absoluteFilePath()) || isOpen;

    if(!isOpen)
      window.fileNew();

    return app.exec();
  }
  catch(std::exception &e)
  {
    qWarning("%s", e.what());
    QMessageBox::critical(nullptr, "GROOPS", QString("****** Error ******\n")+e.what());
  }
  catch(...)
  {
    qWarning("Unknown ERROR");
    QMessageBox::critical(nullptr, "GROOPS", "Unknown ERROR");
  }

  return 0;
}

/***********************************************/
