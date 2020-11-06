/***********************************************/
/**
* @file instrumentInsertNAN.cpp
*
* @brief insert epochs with data NAN into instrument files
*
* @author Matthias Ellmer
* @date 2017-11-28
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program inserts NAN epochs into \configFile{inputfileInstrument}{instrument} files,
either at specific \configClass{times}{timeSeriesType} or where gaps in the instrument are detected.
)";

/***********************************************/

#include "programs/program.h"
#include "classes/timeSeries/timeSeries.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Insert epochs with data NAN into instrument files
* @ingroup programsGroup */
class InstrumentInsertNAN
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentInsertNAN, SINGLEPROCESS, "Insert epochs with data NAN into instrument files", Instrument, Plot)

/***********************************************/

void InstrumentInsertNAN::run(Config &config)
{
  try
  {
    FileName      inputfileInstrument, outputfileInstrument;
    TimeSeriesPtr timeSeries;
    Bool          atGaps;
    Bool          atArcEnds;

    readConfig(config, "outputfileInstrument", outputfileInstrument, Config::MUSTSET, "",  "");
    readConfig(config, "inputfileInstrument",  inputfileInstrument,  Config::MUSTSET, "",  "");
    readConfig(config, "times",                timeSeries,           Config::DEFAULT, "",  "Insert NAN at specific times.");
    readConfig(config, "atGaps",               atGaps,               Config::DEFAULT, "0", "Insert NAN where epochs are more than 1.5 times the median sampling apart.");
    readConfig(config, "atArcEnds",            atArcEnds,            Config::DEFAULT, "0", "Insert one epoch with data NAN at arc ends");
    if(isCreateSchema(config)) return;

    const auto nanTimes = timeSeries->times();

    InstrumentFile inFile(inputfileInstrument);
    const UInt arcCount = inFile.arcCount();

    std::vector<Arc> arcList(arcCount);
    UInt total = 0;
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      Arc arc = inFile.readArc(arcNo);
      if(!arc.size())
        return arc;

      Epoch *epoch = Epoch::create(arc.getType());
      epoch->setData(arc.at(0).data() * NAN_EXPR);

      const auto times = arc.times();
      const Time epsilonTime = Time(0,std::numeric_limits<Double>::epsilon());

      // ----------

      if(nanTimes.size())
      {
        // First value in nanTimes larger than times.begin()
        UInt nanIndex = std::distance(nanTimes.begin(), std::upper_bound(nanTimes.begin(), nanTimes.end(), times.front()));
        auto it = times.begin();

        std::vector<UInt> indicesToInsert;
        for(; nanIndex<nanTimes.size(); nanIndex++)
        {
          // Find the position of the nan epoch in this arc
          it = std::lower_bound(it, times.end(), nanTimes.at(nanIndex));
          if(it == times.end())
            break;

          indicesToInsert.push_back(std::distance(times.begin(), it));
        }

        // Insert starting from back to not screw up the indices.
        std::reverse(indicesToInsert.begin(), indicesToInsert.end());
        for(const UInt i : indicesToInsert)
        {
          epoch->time = times.at(i-1) + epsilonTime;
          arc.insert(i, *epoch);
          total++;
        }
      }

      // ----------

      if(atGaps && times.size() > 2)
      {
        std::vector<Time> diff(times.size());
        std::adjacent_difference(times.begin(), times.end(), diff.begin());

        std::vector<UInt> indicesToInsert;
        const Double minGap = 1.5 * medianSampling(times).seconds();
        for(UInt i=1; i<diff.size(); i++)
        {
          if(diff.at(i).seconds() > minGap)
            indicesToInsert.push_back(i);
        }

        std::reverse(indicesToInsert.begin(), indicesToInsert.end());
        for(const UInt i : indicesToInsert)
        {
          epoch->time = times.at(i-1) + epsilonTime;
          arc.insert(i, *epoch);
          total++;
        }
      }

      // ----------

      if(atArcEnds)
      {
        epoch->time = times.back() + epsilonTime;
        arc.push_back(*epoch);
        total++;
      }

      return arc;
    });

    Parallel::reduceSum(total);
    logStatus<<"  inserted "<<total<<" NaNs"<<Log::endl;

    if(Parallel::isMaster())
    {
      logStatus<<"write instrument to file <"<<outputfileInstrument<<">"<<Log::endl;
      InstrumentFile::write(outputfileInstrument, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
