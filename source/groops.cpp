/***********************************************/
/**
* @file groops.cpp
*
* @brief main.
*
* @author Torsten Mayer-Guerr
* @date 2007
*
*/
/***********************************************/

/**
* @mainpage
*
@verbatim
Gravity Recovery Object Oriented Programming System (GROOPS)
Usage: groops [--log <logfile.txt>] [--settings <groopsDefaults.xml>] [--silent] [--global name=value] <configfile.xml>
       groops --write-settings <groopsDefaults.xml>
       groops --xsd <schemafile.xsd>
       groops --doc <documentation/>

-h, --help           this text
-l, --log            append messages to logfile. If a directory is given, one time-stamped logfile will be created inside for each groops script.
-g, --global         pass a global variable to config files as name=value pair
-c, --settings       read constants from file (default search: groopsDefaults.xml)
-s, --silent         runs silently
-d, --doc            generate documentation files (latex/html/...)
-x, --xsd            write xsd-schema of xml-configfile options
-C, --write-settings write the users current settings to file

GitHub repository: https://github.com/groops-devs/groops
@endverbatim
*/

/** @defgroup base Base
* @brief Basic functions and classes. */

/** @defgroup parserGroup Parser
* @brief XML, mathematical expression and string parser. */

/** @defgroup inputOutputGroup Input/Output
* @brief Low level input/output. */

/** @defgroup parallelGroup Parallel
* @brief Wrapper for Message Passing Interface (MPI). */

/** @defgroup filesGroup File formats
* @brief Definition of groops file formats. */

/** @defgroup configGroup Config
* @brief Config and documentation. */

/** @defgroup classesGroup Classes
* @brief Abstract interfaces. */

/** @defgroup gnssGroup GNSS
* @brief Global Navigation Satellite Systems. */

/** @defgroup plotGroup Plots
* @brief Classes and functions for plotting. */

/** @defgroup miscGroup Misc
* @brief Classes and functions which do not fit into other groups. */

/** @defgroup programsGroup Programs
* @brief The Apps. */

/** @defgroup programsConversionGroup Programs/conversion
* @ingroup programsGroup
* @brief Programs for conversion of strange formats. */

/***********************************************/

#include "programs/program.h"
#include "inputOutput/settings.h"
#include "inputOutput/system.h"
#include "config/generateDocumentation.h"

/***********************************************/

static void groopsHelp(const std::string &progName, Parallel::CommunicatorPtr comm)
{
  if(Parallel::isMaster(comm))
  {
    std::cout<<"Gravity Recovery Object Oriented Programming System (GROOPS)"<<std::endl;
    std::cout<<"Usage: "<<progName<<" [--log <logfile.txt>] [--settings <groopsDefaults.xml>] [--silent] [--global name=value] <configfile.xml>"<<std::endl;
    std::cout<<"       "<<progName<<" --write-settings <groopsDefaults.xml>"<<std::endl;
    std::cout<<"       "<<progName<<" --xsd <schemafile.xsd>"<<std::endl;
    std::cout<<"       "<<progName<<" --doc <documentation/>"<<std::endl;
    std::cout<<std::endl;
    std::cout<<" -h, --help           this text"<<std::endl;
    std::cout<<" -l, --log            append messages to logfile. If a directory is given, one time-stamped logfile will be created inside for each groops script."<<std::endl;
    std::cout<<" -g, --global         pass a global variable to config files as name=value pair"<<std::endl;
    std::cout<<" -c, --settings       read constants from file (default search: groopsDefaults.xml)"<<std::endl;
    std::cout<<" -s, --silent         runs silently"<<std::endl;
    std::cout<<" -d, --doc            generate documentation files (latex/html/...)"<<std::endl;
    std::cout<<" -x, --xsd            write xsd-schema of xml-configfile options"<<std::endl;
    std::cout<<" -C, --write-settings write the users current settings to file"<<std::endl;
    std::cout<<""<<std::endl;
    std::cout<<"GitHub repository: https://github.com/groops-devs/groops"<<std::endl;
    std::cout<<"(Compiled: "<<__DATE__<<" "<<__TIME__<<")"<<std::endl;
  }
  exit(EXIT_FAILURE);
}

/***********************************************/

int main(int argc, char *argv[])
{
  Parallel::CommunicatorPtr comm = Parallel::init(argc, argv);


  try
  {
    Parallel::broadCastExceptions(comm, [&](Parallel::CommunicatorPtr comm)
    {
      // init logging
      // ------------
      Log::setSend(Parallel::addChannel(Log::getRecieve(), comm));
      Log::setRank(Parallel::myRank(comm));
      Log::enableOutput(Parallel::isMaster(comm));

      // handle commandline options
      // --------------------------
      FileName logFileName;
      FileName schemaFileName;
      FileName docFileName;
      FileName settingsFileName;
      FileName writeSettingsFileName;
      Bool     silent   = FALSE;
      Bool     workDone = FALSE;
      std::map<std::string, std::string> commandlineGlobals;
      std::vector<FileName> configFileNames;

      for(int i=1; i<argc; i++)
      {
        auto optArg = [&]()
        {
          if((i+1 >= argc) || (argv[i+1][0] == '-'))
          {
            logWarningOnce<<"Expected argument for: '"<<argv[i]<<"'"<<Log::endl;
            groopsHelp(argv[0], comm);
            return std::string();
          }
          return std::string(argv[++i]);
        };

        const std::string opt(argv[i]);
        if     ((opt == "-l") || (opt == "--log"))            {logFileName           = FileName(optArg());}
        else if((opt == "-x") || (opt == "--xsd"))            {schemaFileName        = FileName(optArg());}
        else if((opt == "-d") || (opt == "--doc"))            {docFileName           = FileName(optArg());}
        else if((opt == "-c") || (opt == "--settings"))       {settingsFileName      = FileName(optArg());}
        else if((opt == "-C") || (opt == "--write-settings")) {writeSettingsFileName = FileName(optArg());}
        else if((opt == "-s") || (opt == "--silent"))         {silent = TRUE;}
        else if((opt == "-h") || (opt == "--help"))           {groopsHelp(argv[0], comm);}
        else if((opt == "-g") || (opt == "--global"))
        {
          std::string keyVal(optArg());
          std::size_t delim = keyVal.find("=");
          if(delim == std::string::npos || delim+1 == keyVal.size())
          {
            logWarningOnce<<"Unable to parse key-value pair <"<<keyVal<<"> for option '-g'."<<Log::endl;
            groopsHelp(argv[0], comm);
          }
          else
            commandlineGlobals[keyVal.substr(0, delim)] = keyVal.substr(delim+1);
        }
        else if(opt[0] == '-')
        {
          logWarningOnce<<"Unknown option: '"<<opt<<"'"<<Log::endl;
          groopsHelp(argv[0], comm);
        }
        else
        {
          configFileNames.push_back(FileName(opt));
        }
      }

      // ============================================

      // start logging
      // -------------
      Log::setSilent(silent);
      if(!System::isDirectory(logFileName))
        Log::setLogFile(logFileName);
      logStatus<<"=== Starting GROOPS ==="<<Log::endl;

      // read default settings and constants
      // -----------------------------------
      if(!settingsFileName.empty())
      {
        logInfo<<"settings: <"<<settingsFileName<<">"<<Log::endl;
        readFileSettings(settingsFileName);
      }
      else // exists groopsDefaults.xml?
      {
        settingsFileName = "groopsDefaults.xml";
        if(System::exists(settingsFileName))
        {
          logInfo<<"settings: <"<<settingsFileName<<">"<<Log::endl;
          readFileSettings(settingsFileName);
        }
      }

      // writing xsd schema file
      // ------------------------
      if(!schemaFileName.empty())
      {
        workDone = TRUE;
        logStatus<<"writing xsd schema file: <"<<schemaFileName<<">"<<Log::endl;
        if(Parallel::isMaster(comm))
          Config::writeSchema(schemaFileName);
      }

      // generate documentation
      // ----------------------
      if(!docFileName.empty())
      {
        workDone = TRUE;
        logStatus<<"generate documentation files in <"<<docFileName<<">"<<Log::endl;
        if(Parallel::isMaster(comm))
          Documentation::write(docFileName);
      }

      // write settings
      // --------------
      if(!writeSettingsFileName.empty())
      {
        workDone = TRUE;
        logStatus<<"writing settings file: <"<<writeSettingsFileName<<">"<<Log::endl;
        if(Parallel::isMaster(comm))
          writeFileSettings(writeSettingsFileName);
      }

      // Starting Programs
      // -----------------
      for(auto &configFileName : configFileNames)
      {
        // If the user specifies a directory as the logging target,
        // a time-stamped log file is created under that directory for each groops script,
        // and log output for that script is redirected there.
        if(System::isDirectory(logFileName))
        {
          FileName thisLogFileName = logFileName.append(configFileName.stripDirectory().str()+"_"+System::now().dateTimeStr()+".log");
          logInfo<<"Future logs are written to file <"<<thisLogFileName<<">"<<Log::endl;
          Log::setLogFile(thisLogFileName);
        }

        logInfo<<"Config file: <"<<configFileName<<">"<<Log::endl;
        Config config(configFileName, commandlineGlobals);
        ProgramConfig programs;
        readConfig(config, "program", programs, Config::OPTIONAL, "", "");
        programs.run(config.getVarList(), comm);
        workDone = TRUE;
      }

      // ============================================

      if(!workDone)
        groopsHelp(argv[0], comm);

      Parallel::barrier(comm);
      logStatus<<"=== Finished GROOPS ==="<<Log::endl;
      Parallel::barrier(comm);
    }); // Parallel::broadCastExceptions()
  }
  catch(std::exception &e)
  {
    if(Parallel::isMaster(comm))
    {
      logError<<"****** Error ******"<<Log::endl;
      logError<<e.what()<<Log::endl;
    }
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/***********************************************/
