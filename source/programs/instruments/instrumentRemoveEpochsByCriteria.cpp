/***********************************************/
/**
* @file instrumentRemoveEpochsByCriteria.cpp
*
* @brief Remove epochs through evaluating expressions
*
* @author Norbert Zehentner
* @author Beate Klinger
* @author Matthias Ellmer
* @author Andreas Kvas
* @date 2015-05-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program removes epochs from \configFile{inputfileInstrument}{instrument}
by evaluating a set of \config{removalCriteria} expressions. For the data
columns the standard data variables are available,
see~\reference{dataVariables}{general.parser:dataVariables}.

The instrument data can be reduced by data from \configFile{inputfileInstrumentReference}{instrument}
prior to evaluation of the expressions.

To reduce the data by its median, use an expression like \verb|data1-data1mean|.
To remove epochs that deviate by more than 3 sigma use \verb|abs(data1)>3*data1std|
or \verb|abs(data0-data0median)>3*1.4826*data0mad|.

All arcs in the input instrument file are concatenated, meaning expressions
like \verb|data1mean| refer to the complete dataset. The removed epochs can be saved
in a separate \configFile{outputfileInstrumentRemovedEpochs}{instrument}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Remove epochs which meet one of several criteria.
* @ingroup programsGroup */
class InstrumentRemoveEpochsByCriteria
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentRemoveEpochsByCriteria, SINGLEPROCESS, "Remove epochs which criteria defined through expressions", Instrument)

/***********************************************/

void InstrumentRemoveEpochsByCriteria::run(Config &config)
{
  try
  {
    FileName fileNameOut, fileNameRemovedEpochs;
    FileName fileNameIn,  fileNameRef;
    std::vector<ExpressionVariablePtr> expressions;
    Double   margin;

    readConfig(config, "outputfileInstrument",               fileNameOut,           Config::OPTIONAL, "", "all data is stored in one arc");
    readConfig(config, "outputfileInstrumentRemovedEpochs",  fileNameRemovedEpochs, Config::OPTIONAL, "", "all data is stored in one arc");
    readConfig(config, "inputfileInstrument",                fileNameIn,            Config::MUSTSET,  "", "arcs are concatenated for processing");
    readConfig(config, "inputfileInstrumentReference",       fileNameRef,           Config::OPTIONAL, "", "if given, the reference data is reduced prior to the expressions being evaluated");
    readConfig(config, "removalCriteria",                    expressions,           Config::MUSTSET,  "abs(data0-data0median) > 3*1.4826*data0mad",  "epochs are removed if one criterion evaluates true. data0 is the first data field.");
    readConfig(config, "margin",                             margin,                Config::DEFAULT,  "1e-5", "remove data around identified epochs (on both sides) [seconds]");
    if(isCreateSchema(config)) return;

    // ======================================================

    // read instrument data
    // --------------------
    logStatus<<"read instrument data <"<<fileNameIn<<">"<<Log::endl;
    const Arc arc = InstrumentFile::read(fileNameIn);
    Matrix data = arc.matrix();
    if(data.columns()<2)
      throw(Exception("input file <"+fileNameIn.str()+"> does not have enough columns (Found: "+data.columns()%"%i)."s));


    // reduce reference data
    // ---------------------
    if(!fileNameRef.empty())
    {
      logStatus<<"read reference data <"<<fileNameRef<<">"<<Log::endl;
      const Arc arcRef = InstrumentFile::read(fileNameRef);
      Arc::checkSynchronized({arc, arcRef});

      // Reduce reference data
      Matrix reference = arcRef.matrix();
      if(reference.columns() != data.columns())
        throw(Exception("Reference file <"+fileNameRef.str()+"> has incompatible number of data columns with data file <"+fileNameIn.str()+"> ("+ reference.columns() % "%i vs "s + data.columns() % "%i )."s));
      axpy(-1., reference.column(1, reference.columns()-1), data.column(1, data.columns()-1) );
    }
    // Remove time column
    data = data.column(1, data.columns()-1);

    // ======================================================

    // initialize data variables
    // -------------------------
    auto varList = config.getVarList();
    std::set<std::string> usedVariables;
    for(UInt i=0; i<expressions.size(); i++)
      expressions.at(i)->usedVariables(varList, usedVariables);
    addDataVariables(data, varList, usedVariables);
    for(UInt i=0; i<expressions.size(); i++)
      expressions.at(i)->simplify(varList);

    // criteria evaluation
    // -------------------
    std::vector<Time> times;
    logStatus<<"evaluate criteria"<<Log::endl;
    logTimerStart;
    for(UInt idEpoch=0; idEpoch<data.rows(); idEpoch++)
    {
      logTimerLoop(idEpoch, data.rows());
      evaluateDataVariables(data, idEpoch, varList);
      for(UInt i=0; i<expressions.size(); i++)
        if(expressions.at(i)->evaluate(varList))
        {
          times.push_back(arc.at(idEpoch).time);
          break;
        }
    }
    logTimerLoopEnd(data.rows());

    // ======================================================

    // remove epochs within buffer
    // ---------------------------
    Arc arcNew, arcEpochsForRemoval;
    if(times.size())
    {
      logStatus<<"remove epochs (+/- "<<margin<<" sec) from instrument data"<<Log::endl;
      UInt idxTime = 0;

      logTimerStart;
      for(UInt i=0; i<arc.size(); i++)
      {
        logTimerLoop(i, arc.size());

        while((idxTime < times.size()) && ((times.at(idxTime)-arc.at(i).time).seconds() < -margin))
          idxTime++;
        if((idxTime >= times.size()) || ((times.at(idxTime)-arc.at(i).time).seconds() > +margin))
          arcNew.push_back(arc.at(i));
        else
          arcEpochsForRemoval.push_back(arc.at(i));
      }
      logTimerLoopEnd(arc.size());
    }
    else
    {
      arcNew = arc;
    }

    logInfo<<"  "<<arcEpochsForRemoval.size()<<" epochs meet criteria"<<Log::endl;
    logInfo<<"  "<<arc.size()-arcNew.size()<<" epochs removed"<<Log::endl;

    // ======================================================

    // save file
    // ---------
    if((!fileNameOut.empty()) && ((times.size()!=0) || (fileNameIn.str() != fileNameOut.str())))
    {
      logStatus<<"write instrument data to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcNew);
      Arc::printStatistics(arcNew);
    }
    if(!fileNameRemovedEpochs.empty())
    {
      logStatus<<"write removed epochs to file <"<<fileNameRemovedEpochs<<">"<<Log::endl;
      InstrumentFile::write(fileNameRemovedEpochs, arcEpochsForRemoval);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
