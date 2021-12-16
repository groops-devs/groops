/***********************************************/
/**
* @file instrumentSynchronize.cpp
*
* @brief Synchronize instrument data.
*
* @author Torsten Mayer-Guerr
* @date 2001-06-08
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads several \file{instrument files}{instrument} and synchronize the data.
Every epoch with some missing data will be deleted so the remaining epochs
have data from every instrument.

In a second step the epochs are divided into arcs with maximal epochs
(or \config{maxArcLen}) without having a gap inside an arc.
A Gap is defined by a time step with at least \config{minGap} seconds
between consecutive epochs or if not set the 1.5 of the median sampling.
Arc with an epoch count less than \config{minArcLen} will be rejected.

A specific region can be selected with \configClass{border}{borderType}.
In this case one of the instrument data must an orbit.

If \configClass{timeIntervals}{timeSeriesType} is given the data are also divided into time bins.
The assignment of arcs to the bins can be saved in \configFile{outputfileArcList}{arcList}.
This file can be used for the variational equation approach or \program{KalmanBuildNormals}.

Instrument files from \config{irregularData} are not synchronized but
divided into the same number of arcs within the same time intervals.
Data outside the defined arcs will be deleted.
)";


/***********************************************/

#include "programs/program.h"
#include "base/planets.h"
#include "files/fileArcList.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/border/border.h"

/***** CLASS ***********************************/

/** @brief Synchronize instrument data.
* @ingroup programsGroup */
class InstrumentSynchronize
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);

  class Data
  {
    public:
    FileName inName, outName;
  };

  class DataIrregular
  {
    public:
    FileName inName, outName;
    UInt     minArcLen;
  };

private:
  enum ArcType {ALL, ASCENDING, DESCENDING};

  OrbitArc          orbitArc;
  std::vector<Time> times;
  std::vector<Time> timesInterval;
  ArcType           arcType;
  BorderPtr         border;
  UInt              searchInterval(UInt i);
};

GROOPS_REGISTER_PROGRAM(InstrumentSynchronize, SINGLEPROCESS, "Synchronize instrument data", Instrument)
GROOPS_RENAMED_PROGRAM(ArcSynchronize, InstrumentSynchronize, date2time(2020, 05, 25))

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, InstrumentSynchronize::Data &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "outputfileInstrument", var.outName, Config::OPTIONAL, "", "");
  readConfig(config, "inputfileInstrument",  var.inName,  Config::MUSTSET,  "", "");
  endSequence(config);
  return TRUE;
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, InstrumentSynchronize::DataIrregular &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "outputfileInstrument", var.outName,   Config::OPTIONAL, "",  "");
  readConfig(config, "inputfileInstrument",  var.inName,    Config::MUSTSET,  "",  "");
  readConfig(config, "minArcLength",         var.minArcLen, Config::DEFAULT,  "0", "minimal number of epochs in an arc");
  endSequence(config);
  return TRUE;
}

/***********************************************/

UInt InstrumentSynchronize::searchInterval(UInt i)
{
  // in time interval?
  UInt idx = 0;
  if(timesInterval.size())
  {
    if(times.at(i)<timesInterval.at(idx))
      return MAX_UINT;
    while((idx+1 < timesInterval.size()) && (times.at(i) >= timesInterval.at(idx+1)))
      idx++;
    if(idx+1 >= timesInterval.size())
      return MAX_UINT;
  }

  // ascending / descending
  if((arcType!=ALL)&&(i>0))
  {
    if((arcType == ASCENDING) && (orbitArc.at(i).position.z() < orbitArc.at(i-1).position.z()))
      return MAX_UINT;
    if((arcType == DESCENDING) && (orbitArc.at(i).position.z() > orbitArc.at(i-1).position.z()))
      return MAX_UINT;
  }

  // in area?
  if(!border || (border->isInnerPoint(Planets::celestial2TerrestrialFrame(times.at(i)).rotate(orbitArc.at(i).position))))
    return idx;
  return MAX_UINT;
}

/***********************************************/

void InstrumentSynchronize::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    std::vector<Data>          data;
    std::vector<DataIrregular> data2;
    Double            minGap = NAN_EXPR, margin;
    UInt              minArcLen = 1, maxArcLen = MAX_UINT;
    TimeSeriesPtr     timeSeries;
    FileName          outArcList;
    arcType  = ALL;

    readConfig(config, "data",              data,       Config::MUSTSET,  "",     "");
    readConfig(config, "margin",            margin,     Config::DEFAULT,  "1e-5", "margin for identical times [seconds]");
    readConfig(config, "minGap",            minGap,     Config::OPTIONAL, "",     "minimal time to define a gap and to begin a new arc, 0: no dividing [seconds], if not set 1.5*median sampling is used");
    readConfig(config, "minArcLength",      minArcLen,  Config::DEFAULT,  "1",    "minimal number of epochs of an arc");
    readConfig(config, "maxArcLength",      maxArcLen,  Config::OPTIONAL, "",     "maximal number of epochs of an arc");
    std::string choice;
    if(readConfigChoice(config, "arcType", choice, Config::OPTIONAL, "", "all arcs or only ascending or descending arcs are selected"))
    {
      if(readConfigChoiceElement(config, "ascending",  choice, "")) arcType = ASCENDING;
      if(readConfigChoiceElement(config, "descending", choice, "")) arcType = DESCENDING;
      endChoice(config);
    }
    readConfig(config, "border",            border,     Config::OPTIONAL, "", "only data in a specific region is selected");
    readConfig(config, "timeIntervals",     timeSeries, Config::DEFAULT,  "", "divide data into time bins");
    readConfig(config, "outputfileArcList", outArcList, Config::OPTIONAL, "", "arc and time bin mapping");
    readConfig(config, "irregularData",     data2,      Config::OPTIONAL, "", "instrument files with irregular sampling");
    if(isCreateSchema(config)) return;

    // =============================================

    // read data
    // ---------
    std::vector<Arc> arc(data.size());
    for(UInt k=0; k<data.size(); k++)
    {
      logStatus<<"read instrument file <"<<data.at(k).inName<<">"<<Log::endl;
      arc.at(k) = InstrumentFile::read(data.at(k).inName);
      logInfo<<"  epochs = "<<arc.at(k).size()<<Log::endl;
    }

    std::vector<Arc> arcIrregular(data2.size());
    for(UInt k=0; k<data2.size(); k++)
    {
      logStatus<<"read irregular instrument file <"<< data2.at(k).inName<<Log::endl;
      arcIrregular.at(k) = InstrumentFile::read(data2.at(k).inName);
      logInfo<<"  epochs = "<<arcIrregular.at(k).size()<<Log::endl;
    }

    // find orbit data
    // ---------------
    UInt indexOrbit = NULLINDEX;
    if(border || (arcType != ALL))
    {
      for(UInt k=0; k<arc.size(); k++)
        if(arc.at(k).getType() == Epoch::ORBIT)
        {
          indexOrbit = k;
          break;
        }
      if(indexOrbit == NULLINDEX)
        throw(Exception("One instrument file must be an orbit."));
    }

    // synchronize data
    // ----------------
    logStatus<<"synchronize data"<<Log::endl;
    times.resize(arc.at(0).size());
    std::vector<UInt> index(arc.size(),0);
    for(UInt i=0; i<arc.at(0).size(); i++)
    {
      Time time = arc.at(0).at(i).time;

      // this point of time in all files?
      Bool synchron = TRUE;
      Bool eof      = FALSE;
      for(UInt k=1; k<arc.size(); k++)
      {
        while(((arc.at(k).at(index.at(k)).time-time).seconds() < -margin) && (++index.at(k)<arc.at(k).size()));

        if(index.at(k)>=arc.at(k).size())
        {
          eof = TRUE;
          break;
        }
        if(std::fabs((arc.at(k).at(index.at(k)).time-time).seconds()) > margin)
        {
          synchron = FALSE;
          break;
        }
      }
      if(eof)       break;
      if(!synchron) continue;
      times.at(index.at(0)++) = time;
    }
    times.resize(index.at(0));

    if(times.size()==0)
    {
      logWarning<<"found no data with identical time"<<Log::endl;
      return;
    }

    // delete other data
    // -----------------
    logStatus<<"delete asynchronize data"<<Log::endl;
    for(UInt k=0; k<arc.size(); k++)
      arc.at(k).synchronize(times, margin);

    if(indexOrbit != NULLINDEX)
      orbitArc = arc.at(indexOrbit);

    // divide arcs
    // -----------
    logStatus<<"divide arcs"<<Log::endl;
    // determine min gap size
    Time timeGap;
    if(!std::isnan(minGap))
      timeGap = (minGap>0) ? seconds2time(minGap) : seconds2time(100*365*86400.);
    else
    {
      // median sampling
      std::vector<Time> timeDiff(times.size()-1);
      for(UInt i=0; i<times.size()-1; i++)
        timeDiff.at(i) = times.at(i+1)-times.at(i);
      std::sort(timeDiff.begin(),timeDiff.end());
      timeGap = 1.5*timeDiff.at(timeDiff.size()/2);
      logInfo<<"  begin new arc by a time step of at least "<<timeGap.seconds()<<" seconds"<<Log::endl;
    }
    minArcLen = std::max(minArcLen, static_cast<UInt>(1));
    timesInterval = timeSeries->times();
    std::vector<UInt> arcsInterval(timesInterval.size(), 0);

    std::vector<UInt> subArcStart, subArcLen;
    std::vector< std::vector<UInt> > irregularSubArcStart(data2.size()), irregularSubArcLen(data2.size());

    UInt idx = 0;
    while((times.size()-idx) >= minArcLen)
    {
      // search state change
      UInt idxStart    = idx;
      UInt idxInterval = searchInterval(idx++);
      while((idx<times.size()) && (times.at(idx)-times.at(idx-1) < timeGap) && ((idx-idxStart)<maxArcLen) && (idxInterval == searchInterval(idx)))
        idx++;

      // reject arc?
      if((idxInterval == MAX_UINT) || ((idx-idxStart) < minArcLen))
        continue;

      // sort irregular data into arcs
      Bool shortArc = FALSE;
      for(UInt k=0; k<data2.size(); k++)
      {
        UInt idx2 = 0;
        while((idx2<arcIrregular.at(k).size()) && ((arcIrregular.at(k).at(idx2).time-times.at(idxStart)).seconds() < -margin))
          idx2++;
        UInt idx2Start = idx2;
        while((idx2<arcIrregular.at(k).size()) && ((arcIrregular.at(k).at(idx2).time-times.at(idx-1)).seconds() < +margin))
          idx2++;
        irregularSubArcStart.at(k).push_back(idx2Start);
        irregularSubArcLen.at(k).push_back(idx2-idx2Start);
        if(idx2-idx2Start < data2.at(k).minArcLen)
          shortArc = TRUE;
      }

      // reject arc when too few irregular data epochs were found in the arc
      if(shortArc)
      {
        for(UInt k=0; k<data2.size(); k++)
        {
          irregularSubArcStart.at(k).pop_back();
          irregularSubArcLen.at(k).pop_back();
        }
        continue;
      }

      // valid arc -> save
      subArcStart.push_back(idxStart);
      subArcLen.push_back(idx-idxStart);

      if(arcsInterval.size())
        for(UInt i=idxInterval+1; i<arcsInterval.size(); i++)
          arcsInterval.at(i) = subArcLen.size();
    } // end while((times.size()-idx) >= minArcLen)

    if(subArcStart.size()==0)
    {
      logWarning<<"no arcs found"<<Log::endl;
      return;
    }

    // =============================================

    // save files
    // ----------
    if(!outArcList.empty())
    {
      logStatus<<"write arc list <"<<outArcList<<">"<<Log::endl;
      writeFileArcList(outArcList, arcsInterval, timesInterval);
    }

    std::vector<Arc> arcList;
    for(UInt k=0; k<data.size(); k++)
    {
      if(!data.at(k).outName.empty())
      {
        logStatus<<"write instrument file <"<<data.at(k).outName<<">"<<Log::endl;
        arcList.resize(0);
        for(UInt i=0; i<subArcStart.size(); i++)
          arcList.push_back( arc.at(k).subArc(subArcStart.at(i), subArcLen.at(i)) );
        InstrumentFile::write(data.at(k).outName, arcList);
      }
    }

    // write irregular data to file
    for(UInt k=0; k<data2.size(); k++)
    {
      if(!data2.at(k).outName.empty())
      {
        logStatus<<"write irregular instrument file <"<<data2.at(k).outName<<">"<<Log::endl;
        std::vector<Arc> arcIrregularList;
        for(UInt i=0; i<irregularSubArcStart.at(k).size(); i++)
          arcIrregularList.push_back( arcIrregular.at(k).subArc(irregularSubArcStart.at(k).at(i), irregularSubArcLen.at(k).at(i)) );
        InstrumentFile::write(data2.at(k).outName, arcIrregularList);
      }
    }

    Arc::printStatistics(arcList);
    for(UInt i=1; i<arcsInterval.size(); i++)
      if(arcsInterval.at(i-1) == arcsInterval.at(i))
        logWarning << i<<". intervall ("<<timesInterval.at(i-1).dateTimeStr()<<" - "<<timesInterval.at(i).dateTimeStr()<<") is empty"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
