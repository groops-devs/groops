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
#include <iostream>
#include "base/importGroops.h"
#include "base/schema.h"
#include "mainWindow/mainWindow.h"

/***********************************************/

class MyApplication : public QApplication
{
public:
  MyApplication(int &argc, char **argv) : QApplication(argc, argv) {}
  virtual ~MyApplication() {}

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

    // Handle command line options
    parser.process(app);
    const Bool cleanStartup = parser.isSet("clean");
    QFileInfoList fileNames;
    for(const auto &file : parser.values("file"))
      fileNames.push_back(file);
    QFileInfo schema = parser.isSet("schema") ? parser.value("schema") : QFileInfo();

    // Set XML schema if passed as argument
    const bool isValidSchema = schema.isFile() ? Schema::validateSchema(schema.absoluteFilePath()) : false;
    if(isValidSchema)
    {
      QSettings *settings = new QSettings();
      settings->setValue("files/schemaFile", schema.absoluteFilePath());
      QStringList schemaFiles = settings->value("files/schemaFiles").toStringList();
      if(!schemaFiles.contains(schema.absoluteFilePath()))
      {
        schemaFiles.append(schema.absoluteFilePath());
        schemaFiles.sort();
        settings->setValue("files/schemaFiles", schemaFiles);
      }
    }

    MainWindow window;
    window.show();

    if(!schema.absoluteFilePath().isEmpty() && !isValidSchema)
      QMessageBox::critical(nullptr, "GROOPS", QString("File '%1' passed as argument seems not to be a valid XSD schema").arg(schema.absoluteFilePath()));

    // Open files from command line arguments
    if(!cleanStartup)
      window.fileOpenInitial();
    else if(fileNames.size() == 0)
      window.fileNew();

    for(int i = 0; i < fileNames.size(); i++)
      if(fileNames.at(i).isFile())
        window.fileOpen(fileNames.at(i).absoluteFilePath());

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
