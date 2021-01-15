/***********************************************/
/**
* @file GnssEstimateClockShift.cpp
*
* @brief Estimate clock shift for a constellation of satellites.
*
* @author Sebastian Strasser
* @date 2018-02-06
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring =
R"(
This program estimates an epoch-wise clock shift in a constellation of GNSS satellites.
Each separate \config{data} represents a satellite... (e.g. 32 GPS satellites).
The shift to reference clocks can be estimated by providing \configFile{inputfileInstrumentRef}{instrument}.
Clock shifts are estimated for each epoch given by \configClass{timeSeries}{timeSeriesType}.
)";

/***********************************************/

#include "base/import.h"
#include "classes/timeSeries/timeSeries.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "parallel/parallel.h"
#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief Estimate clock shift for a constellation of satellites.
* @ingroup programsGroup */
class GnssEstimateClockShift
{
public:
  class Data
  {
  public:
    FileName outNameInstrument, outNameInstrumentDiff, inNameInstrument, inNameInstrumentRef;
    std::vector<Time> times;
    std::vector<UInt> index;
    Vector clock, clockDiff;
    Bool   useable;

    void init(const std::vector<Time> &timesGlobal, Double margin);
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssEstimateClockShift, SINGLEPROCESS, "Estimate clock shift for a constellation of satellites.", Gnss)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssEstimateClockShift::Data &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "outputfileInstrument",     var.outNameInstrument,     Config::OPTIONAL, "", "corrected clocks");
  readConfig(config, "outputfileInstrumentDiff", var.outNameInstrumentDiff, Config::OPTIONAL, "", "clock difference after correction");
  readConfig(config, "inputfileInstrument",      var.inNameInstrument,      Config::MUSTSET,  "", "input clocks");
  readConfig(config, "inputfileInstrumentRef",   var.inNameInstrumentRef,   Config::OPTIONAL, "", "reference clocks (subtracted from input clocks)");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void GnssEstimateClockShift::Data::init(const std::vector<Time> &timesGlobal, Double margin)
{
  try
  {
    useable = FALSE;

    MiscValueArc arc, arcRef;
    try
    {
      arc = InstrumentFile::read(inNameInstrument);
      if(!inNameInstrumentRef.empty())
        arcRef = InstrumentFile::read(inNameInstrumentRef);
    }
    catch(std::exception &e)
    {
      logWarning << e.what() << Log::endl;
      return;
    }

    if(!arc.size() || (!inNameInstrumentRef.empty() && arc.times() != arcRef.times()))
    {
      logWarning << "arcs empty or not synchronized: " << arc.size() << " != " << arcRef.size() << ", skipping" << Log::endl;
      return;
    }

    times = arc.times();
    clock = clockDiff = arc.matrix().column(1);
    if(arcRef.size())
      clockDiff -= arcRef.matrix().column(1);

    // find indices of given epochs
    Time timeMargin = seconds2time(margin);
    index.resize(timesGlobal.size(), MAX_UINT);
    UInt idEpochGlobal = 0;
    for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
    {
      while(idEpochGlobal+1 < timesGlobal.size() && timesGlobal.at(idEpochGlobal+1) < times.at(idEpoch)+timeMargin)
        idEpochGlobal++;
      if(timesGlobal.at(idEpochGlobal) >= times.at(idEpoch)-timeMargin && timesGlobal.at(idEpochGlobal) <= times.at(idEpoch)+timeMargin)
        index.at(idEpochGlobal) = idEpoch;
    }

    useable = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssEstimateClockShift::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName          fileNameShiftTimeSeries;
    std::vector<Data> data;
    TimeSeriesPtr     timeSeries;
    Double            margin;

    readConfig(config, "outputfileShiftTimeSeries", fileNameShiftTimeSeries, Config::OPTIONAL, "",    "columns: mjd, clock shift");
    readConfig(config, "data",                      data,                    Config::MUSTSET,  "",    "e.g. satellite");
    readConfig(config, "timeSeries",                timeSeries,              Config::MUSTSET,  "",    "clock epochs");
    readConfig(config, "margin",                    margin,                  Config::DEFAULT,  "0.1", "[s] margin for time comparison");
    if(isCreateSchema(config)) return;

    const std::vector<Time> times = timeSeries->times();

    // initialize clock data
    for(auto &&d : data)
    {
      logStatus << "read instrument file <" << d.inNameInstrument << ">" << Log::endl;
      d.init(times, margin);
    }
    data.erase(std::remove_if(data.begin(), data.end(), [](const Data &d){ return !d.useable; }), data.end());
    if(!data.size())
      throw(Exception("no data found"));

    std::vector<Double> shiftTimeSeries(times.size(), NAN_EXPR);
    Single::forEach(times.size(), [&](UInt idEpoch)
    {
      // build observation vector
      std::vector<Double> l;
      for(const auto &d : data)
        if(d.index.at(idEpoch) != MAX_UINT)
          l.push_back(d.clockDiff(d.index.at(idEpoch)));

      if(l.size() == 0)
      {
        logWarning<<"No data found at epoch "+times.at(idEpoch).dateTimeStr()+". continue with next epoch"<<Log::endl;
        return;
      }

      // compute and remove shift
      shiftTimeSeries.at(idEpoch) = mean(Vector(l));
      for(auto &&d : data)
        if(d.index.at(idEpoch) != MAX_UINT)
        {
          d.clock(d.index.at(idEpoch))     -= shiftTimeSeries.at(idEpoch);
          d.clockDiff(d.index.at(idEpoch)) -= shiftTimeSeries.at(idEpoch);
        }
    });

    // save output files
    for(const auto &d : data)
    {
      if(!d.outNameInstrument.empty())
      {
        logStatus << "write corrected instrument file <" << d.outNameInstrument << ">" << Log::endl;
        Matrix A(d.times.size(), 2);
        copy(d.clock, A.column(1));
        InstrumentFile::write(d.outNameInstrument, Arc(d.times, A, Epoch::MISCVALUE));
      }
      if(!d.outNameInstrumentDiff.empty())
      {
        logStatus << "write clock difference instrument file <" << d.outNameInstrumentDiff << ">" << Log::endl;
        Matrix A(d.times.size(), 2);
        copy(d.clockDiff, A.column(1));
        InstrumentFile::write(d.outNameInstrumentDiff, Arc(d.times, A, Epoch::MISCVALUE));
      }
    }

    if(!fileNameShiftTimeSeries.empty())
    {
      logStatus<<"Write clock shift time series to file <"<<fileNameShiftTimeSeries<<">"<<Log::endl;
      Matrix A(shiftTimeSeries.size(), 2);
      copy(Vector(shiftTimeSeries), A.column(1));
      InstrumentFile::write(fileNameShiftTimeSeries, Arc(times, A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
