/***********************************************/
/**
* @file instrumentMultiplyAdd.cpp
*
* @brief Multiply instrument data with a factor and add them together.
*
* @author Torsten Mayer-Guerr
* @date 2012-06-24
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program multiply \file{instrument data}{instrument} with a factor and add them together.
Afterwards the mean of each arc and data column can be removed with \config{removeArcMean}.
The instrument files must be synchronized (\program{InstrumentSynchronize}).

See also \program{InstrumentArcCalculate}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Multiply instrument data with a factor and add them together.
* @ingroup programsGroup */
class InstrumentMultiplyAdd
{
public:
  class Data
  {
    public:
    FileName fileName;
    Double   factor;
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentMultiplyAdd, PARALLEL, "Multiply instrument data with a factor and add them together", Instrument)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, InstrumentMultiplyAdd::Data &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileInstrument", var.fileName, Config::MUSTSET,  "", "");
  readConfig(config, "factor",              var.factor,   Config::DEFAULT,  "1.0", "");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void InstrumentMultiplyAdd::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOut;
    std::vector<Data> data;
    Bool removeMean;

    readConfig(config, "outputfileInstrument", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "instrument",           data,        Config::MUSTSET,  "", "");
    readConfig(config, "removeArcMean",        removeMean,  Config::DEFAULT,  "0", "remove mean value of each arc");
    if(isCreateSchema(config)) return;

    // open files and check consistency
    std::vector<InstrumentFilePtr> instrumentFile(data.size());
    for(UInt i=0; i<instrumentFile.size(); i++)
    {
      logStatus<<"read instrument data <"<<data.at(i).fileName<<">"<<Log::endl;
      instrumentFile.at(i) = InstrumentFile::newFile(data.at(i).fileName);
      InstrumentFile::checkArcCount({*instrumentFile.at(0), *instrumentFile.at(i)});
      if(instrumentFile.at(i)->getType() != instrumentFile.at(0)->getType())
        throw(Exception("instruments types are different: "+instrumentFile.at(i)->getTypeName()+", "+instrumentFile.at(0)->getTypeName()));
    }

    logStatus<<"combine instrument data"<<Log::endl;
    std::vector<Arc> arcList(instrumentFile.at(0)->arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      Arc arc = instrumentFile.at(0)->readArc(arcNo);
      Matrix A = data.at(0).factor * arc.matrix();
      for(UInt i=1; i<instrumentFile.size(); i++)
      {
        Arc arc2 = instrumentFile.at(i)->readArc(arcNo);
        Arc::checkSynchronized({arc, arc2});
        axpy(data.at(i).factor, arc2.matrix(), A);
      }

      if(removeMean)
        for(UInt k=0; k<A.columns(); k++)
          A.column(k) -= mean(A.column(k));

      return Arc(arc.times(), A, arc.getType());
    }, comm);


    if(Parallel::isMaster(comm))
    {
      logStatus<<"write instrument data to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
