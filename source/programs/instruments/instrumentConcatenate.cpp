/***********************************************/
/**
* @file instrumentConcatenate.cpp
*
* @brief Concatenate arcs from several files.
*
* @author Torsten Mayer-Guerr
* @author Norbert Zehentner
* @date 2001-06-08
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program concatenate the arcs from several \file{instrument files}{instrument}
and write it to a new \file{file}{instrument}. Input files must be of the same type.
The arcs are merged to one arc even though there is a gap inbetween.
To split the data into arcs use \program{InstrumentSynchronize}.
Three options are available: \config{sort}, \config{removeDuplicates} and \config{checkForNaNs}.
If \config{sort} is enabled, the program reads all files, no matter if they are sorted correctly in time, and
then sorts the epochs. If \config{removeDuplicates} is enabled, the program checks the whole data set
for epochs that are contained twice. And if \config{checkForNaNs} is enabled the data set is checked for
invalid epochs containing NaNs.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Concatenate arcs from several files.
* @ingroup programsGroup */
class InstrumentConcatenate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentConcatenate, SINGLEPROCESS, "concatenate arcs from several files", Instrument)
GROOPS_RENAMED_PROGRAM(ArcConcatenate, InstrumentConcatenate, date2time(2020, 05, 25))

/***********************************************/

void InstrumentConcatenate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outName;
    std::vector<FileName> inName;
    Bool                  sort, checkNaN;
    std::string           choiceRemoveDuplicates;
    Double                margin;

    readConfig(config, "outputfile", outName, Config::MUSTSET,  "",  "");
    readConfig(config, "inputfile",  inName,  Config::MUSTSET,  "",  "");
    readConfig(config, "sort",       sort,    Config::DEFAULT,  "0", "sort epochs with increasing time");
    if(readConfigChoice(config, "removeDuplicates", choiceRemoveDuplicates, Config::OPTIONAL, "", "remove duplicate epochs"))
    {
      if(readConfigChoiceElement(config, "keepFirst", choiceRemoveDuplicates, "keep first epoch with the same time stamp, remove all others"))
        readConfig(config, "margin", margin, Config::DEFAULT, "1e-5", "margin for identical times [seconds]");
      if(readConfigChoiceElement(config, "keepLast",  choiceRemoveDuplicates, "keep last epoch with the same time stamp, remove all others"))
        readConfig(config, "margin", margin, Config::DEFAULT, "1e-5", "margin for identical times [seconds]");
      endChoice(config);
    }
    readConfig(config, "checkForNaNs", checkNaN, Config::DEFAULT,  "0", "remove epochs with NaN values in one of the data fields");
    if(isCreateSchema(config)) return;

    // read data
    // ---------
    Arc arc;
    for(UInt i=0; i<inName.size(); i++)
    {
      try
      {
        logStatus<<"read instrument file <"<<inName.at(i)<<">"<<Log::endl;
        arc.append(InstrumentFile::read(inName.at(i)));
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<" continue..."<<Log::endl;
        continue;
      }
    }

    // sort data
    // ---------
    if(sort || (choiceRemoveDuplicates == "keepFirst" || choiceRemoveDuplicates == "keepLast"))
    {
      logStatus<<"sort epochs"<<Log::endl;
      arc.sort();
    }

    // eliminate duplicates
    // --------------------
    if(choiceRemoveDuplicates == "keepFirst" || choiceRemoveDuplicates == "keepLast")
    {
      logStatus<<"eliminate duplicates"<<Log::endl;
      UInt oldSize = arc.size();
      arc.removeDuplicateEpochs(choiceRemoveDuplicates == "keepFirst", margin);
      logInfo<<" "<<oldSize-arc.size()<<" duplicates removed!"<<Log::endl;
    }

    // eliminate NaNs
    // --------------------
    if(checkNaN)
    {
      logStatus<<"search for NaNs"<<Log::endl;
      UInt removed=0;
      logTimerStart;
      for(UInt i=0; i<arc.size(); i++)
      {
        logTimerLoop(i,arc.size());
        Vector data=arc.at(i).data();
        for(UInt j=0; j<data.rows(); j++)
          if(std::isnan(data.at(j)))
          {
            arc.remove(i);
            removed++;
            i--;
            break;
          }
      }
      logTimerLoopEnd(arc.size());
      logInfo<<" "<<removed<<" epochs with NaN values removed!"<<Log::endl;
    }

    // save
    // ----
    logStatus<<"write instrument file <"<<outName<<">"<<Log::endl;
    InstrumentFile::write(outName, arc);
    Arc::printStatistics(arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
