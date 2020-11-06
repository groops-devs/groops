/***********************************************/
/**
* @file instrumentEstimateHelmertTransformation.cpp
*
* @brief Estimate Helmert transformation parameters between two networks (frame realizations).
*
* @author Sebastian Strasser
* @date 2018-01-31
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates a 3D Helmert transformation between two networks
(frame realizations, e.g. GNSS satellite or station network).
Each separate \config{data} represents a satellite/station/\ldots (e.g. 32 GPS satellites).
The instrument data (x,y,z position) considered can be set with \config{startData}.
The Helmert parameters are set up according to \configClass{parametrizationTemporal}{parametrizationTemporalType}
for each \configClass{timeIntervals}{timeSeriesType} and are estimated using a
\reference{robust least squares adjustment}{fundamentals.robustLeastSquares}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "misc/varianceComponentEstimation.h"

/***** CLASS ***********************************/

/** @brief Estimate Helmert transformation parameters between two networks (frame realizations).
* @ingroup programsGroup */
class InstrumentEstimateHelmertTransformation
{
public:
  class Data
  {
  public:
    FileName outNameInstrumentTransformed, outNameInstrumentDiff, inNameInstrument, inNameInstrumentRef;
    std::vector<Time> times;
    std::vector<UInt> intervalBoundaries;
    Matrix pos, posRef;
    UInt   startData;
    Bool   useable;

    void init(const std::vector<Time> &intervals);
    void removeEpochsOutsideIntervals(const FileName &fileName, const std::vector<Time> &intervals, Arc &arc) const;
  };

  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentEstimateHelmertTransformation, SINGLEPROCESS, "Compute RMS time series from instrument file(s).", Instrument, Gnss)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, InstrumentEstimateHelmertTransformation::Data &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;

  renameDeprecatedConfig(config, "startData", "startDataFields", date2time(2020, 7, 14));

  readConfig(config, "outputfileInstrument",         var.outNameInstrumentTransformed, Config::OPTIONAL, "", "transformed positions as instrument type Vector3d");
  readConfig(config, "outputfileInstrumentDiff",     var.outNameInstrumentDiff,        Config::OPTIONAL, "", "position difference as instrument type Vector3d");
  readConfig(config, "inputfileInstrument",          var.inNameInstrument,             Config::MUSTSET,  "", "");
  readConfig(config, "inputfileInstrumentReference", var.inNameInstrumentRef,          Config::MUSTSET,  "", "");
  readConfig(config, "startDataFields",              var.startData,                    Config::DEFAULT,  "0", "start index of position (x,y,z) columns");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void InstrumentEstimateHelmertTransformation::run(Config &config)
{
  try
  {
    FileName                   fileNameHelmertTimeSeries;
    std::vector<Data>          data;
    TimeSeriesPtr              intervalsPtr;
    Bool                       estimateShift, estimateScale, estimateRotation;
    ParametrizationTemporalPtr parametrization;
    Double                     huber, huberPower;
    UInt                       iterCount;

    renameDeprecatedConfig(config, "temporalRepresentation", "parametrizationTemporal", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "intervals",              "timeIntervals",           date2time(2020, 7, 14));

    readConfig(config, "outputfileHelmertTimeSeries", fileNameHelmertTimeSeries, Config::OPTIONAL, "",  "columns: mjd, Tx,Ty,Tz,s,Rx,Ry,Rz according to temporal parametrization");
    readConfig(config, "data",                        data,                      Config::MUSTSET,  "",  "e.g. satellite, station");
    readConfig(config, "timeIntervals",               intervalsPtr,              Config::MUSTSET,  "",  "parameters are estimated per interval");
    readConfig(config, "parametrizationTemporal",     parametrization,           Config::MUSTSET,  "",  "temporal parametrization");
    readConfig(config, "estimateShift",               estimateShift,             Config::DEFAULT,  "1", "coordinate center");
    readConfig(config, "estimateScale",               estimateScale,             Config::DEFAULT,  "1", "scale factor of position");
    readConfig(config, "estimateRotation",            estimateRotation,          Config::DEFAULT,  "1", "rotation");
    readConfig(config, "huber",                       huber,                     Config::DEFAULT,  "2.5", "for robust least squares");
    readConfig(config, "huberPower",                  huberPower,                Config::DEFAULT,  "1.5", "for robust least squares");
    readConfig(config, "huberMaxIteration",           iterCount,                 Config::DEFAULT,  "30",  "(maximum) number of iterations for robust estimation");
    if(isCreateSchema(config)) return;

    const std::vector<Time> intervals = intervalsPtr->times();
    if(intervals.size() < 2)
      throw(Exception("invalid intervals"));

    for(auto &&d : data)
    {
      logStatus << "read instrument file <" << d.inNameInstrument << ">" << Log::endl;
      d.init(intervals);
    }
    data.erase(std::remove_if(data.begin(), data.end(), [](const Data &d){ return !d.useable; }), data.end());

    const UInt paramCount    = 3*estimateShift+estimateScale+3*estimateRotation;
    const UInt intervalCount = intervals.size() - 1;
    std::vector<Vector> helmertTimeSeries;
    logTimerStart;
    for(UInt idInterval = 0; idInterval < intervalCount; idInterval++)
    {
      logTimerLoop(idInterval, intervalCount);
      parametrization->setInterval(intervals.at(idInterval), intervals.at(idInterval+1));
      const UInt trendCount = parametrization->parameterCount();

      // check if interval has data
      UInt epochCount = 0;
      for(const auto &d : data)
        epochCount += d.intervalBoundaries.at(idInterval+1) - d.intervalBoundaries.at(idInterval);
      if(epochCount == 0)
      {
        logWarning<<"No data found in interval "+intervals.at(idInterval).dateTimeStr()+" -> "+intervals.at(idInterval+1).dateTimeStr()+". continue with next interval"<<Log::endl;
        continue;
      }

      // check if interval is solvable
      if(3*static_cast<UInt>(std::count_if(data.begin(), data.end(), [&](const Data &d){ return ((d.intervalBoundaries.at(idInterval+1)-d.intervalBoundaries.at(idInterval)) >= trendCount); })) < paramCount)
        throw(Exception("interval "+intervals.at(idInterval).dateTimeStr()+" -> "+intervals.at(idInterval+1).dateTimeStr()+" not solvable"));

      // set up observation vector and design matrix
      UInt idRow = 0;
      Vector l(3*epochCount);
      Matrix A(3*epochCount, paramCount*trendCount);
      for(const auto &d : data)
      {
        for(UInt idEpoch=d.intervalBoundaries.at(idInterval); idEpoch<d.intervalBoundaries.at(idInterval+1); idEpoch++)
        {
          copy((d.pos.row(idEpoch) - d.posRef.row(idEpoch)).trans(), l.row(idRow, 3));

          std::vector<UInt>   index;
          std::vector<Double> factor;
          parametrization->factors(d.times.at(idEpoch), index, factor);

          for(UInt k = 0; k < trendCount; k++)
          {
            UInt idCol = index.at(k)*paramCount;
            if(estimateShift)
            {
              A(idRow+0, idCol++) = factor.at(k);
              A(idRow+1, idCol++) = factor.at(k);
              A(idRow+2, idCol++) = factor.at(k);
            }
            if(estimateScale)
            {
              Vector3d position(d.posRef.row(idEpoch));
              A(idRow+0, idCol)   = position.x()/DEFAULT_R*factor.at(k);
              A(idRow+1, idCol)   = position.y()/DEFAULT_R*factor.at(k);
              A(idRow+2, idCol++) = position.z()/DEFAULT_R*factor.at(k);
            }
            if(estimateRotation)
            {
              Vector3d position(d.posRef.row(idEpoch));
              A(idRow+1, idCol)   =  position.z()/DEFAULT_R*factor.at(k); //rx
              A(idRow+2, idCol++) = -position.y()/DEFAULT_R*factor.at(k); //rx
              A(idRow+0, idCol)   = -position.z()/DEFAULT_R*factor.at(k); //ry
              A(idRow+2, idCol++) =  position.x()/DEFAULT_R*factor.at(k); //ry
              A(idRow+0, idCol)   =  position.y()/DEFAULT_R*factor.at(k); //rz
              A(idRow+1, idCol++) = -position.x()/DEFAULT_R*factor.at(k); //rz
            }
          }
          idRow += 3;
        }
      }

      // solve system and compute residuals
      Vector x;
      if(A.size())
      {
        Vector sigma;
        x = Vce::robustLeastSquares(A, l, 3, huber, huberPower, iterCount, sigma);
        l -= A*x;
      }

      // write residuals into input data
      idRow = 0;
      for(const auto &d : data)
        for(UInt idEpoch=d.intervalBoundaries.at(idInterval); idEpoch<d.intervalBoundaries.at(idInterval+1); idEpoch++)
          copy(l.row(3*idRow++, 3).trans(), d.pos.row(idEpoch));

      // add Helmert parameters to time series
      Vector h(1+paramCount*trendCount);
      h(0) = 0.5*(intervals.at(idInterval)+intervals.at(idInterval+1)).mjd();
      if(x.size())
      {
        copy(x, h.row(1, x.size()));
        for(UInt k = 0; k < trendCount; k++)
          if(estimateScale || estimateRotation)
            h.row(1+k*paramCount+3*estimateShift, paramCount-3*estimateShift) /= DEFAULT_R;
      }
      helmertTimeSeries.push_back(h);
    }
    logTimerLoopEnd(intervalCount);

    // save output files
    for(const auto &d : data)
      if(!d.outNameInstrumentTransformed.empty())
      {
        logStatus << "write transformed instrument file <" << d.outNameInstrumentTransformed << ">" << Log::endl;
        Matrix A(d.times.size(), 4);
        copy(d.pos+d.posRef, A.column(1,3));
        InstrumentFile::write(d.outNameInstrumentTransformed, Arc(d.times, A, Epoch::VECTOR3D));
      }

    for(const auto &d : data)
      if(!d.outNameInstrumentDiff.empty())
      {
        logStatus << "write difference instrument file <" << d.outNameInstrumentDiff << ">" << Log::endl;
        Matrix A(d.times.size(), 4);
        copy(d.pos, A.column(1,3));
        InstrumentFile::write(d.outNameInstrumentDiff, Arc(d.times, A, Epoch::VECTOR3D));
      }

    if(!fileNameHelmertTimeSeries.empty() && (estimateShift || estimateScale || estimateRotation))
    {
      logStatus<<"write Helmert time series to file <"<<fileNameHelmertTimeSeries<<">"<<Log::endl;
      Matrix H(helmertTimeSeries.size(), 1+paramCount*parametrization->parameterCount());
      for(UInt i=0; i<H.rows(); i++)
        copy(helmertTimeSeries.at(i).trans(), H.row(i));
      writeFileMatrix(fileNameHelmertTimeSeries, H);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentEstimateHelmertTransformation::Data::init(const std::vector<Time> &intervals)
{
  try
  {
    useable = FALSE;

    Arc arc, arcRef;
    try
    {
      arc    = InstrumentFile::read(inNameInstrument);
      arcRef = InstrumentFile::read(inNameInstrumentRef);
    }
    catch(std::exception &e)
    {
      logWarning << e.what() << Log::endl;
      return;
    }
    removeEpochsOutsideIntervals(inNameInstrument,    intervals, arc);
    removeEpochsOutsideIntervals(inNameInstrumentRef, intervals, arcRef);

    if(!arc.size() || arc.times() != arcRef.times())
    {
      logWarning << "arcs empty or not synchronized: " << arc.size() << " != " << arcRef.size() << ", skipping" << Log::endl;
      return;
    }

    times  = arc.times();
    pos    = arc.matrix().column(1+startData,3);
    posRef = arcRef.matrix().column(1+startData,3);

    // find epochs IDs of interval boundaries
    UInt idEpoch = 0;
    intervalBoundaries.push_back(idEpoch);
    for(UInt idInterval = 0; idInterval < intervals.size()-1; idInterval++)
    {
      while(idEpoch < times.size() && times.at(idEpoch) < intervals.at(idInterval+1))
        idEpoch++;
      intervalBoundaries.push_back(idEpoch);
    }

    useable = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentEstimateHelmertTransformation::Data::removeEpochsOutsideIntervals(const FileName &fileName, const std::vector<Time> &intervals, Arc &arc) const
{
  UInt removedEpochs = arc.size();
  for(UInt idEpoch = arc.size(); idEpoch --> 0; )
    if(arc.at(idEpoch).time < intervals.front() || arc.at(idEpoch).time >= intervals.back())
      arc.remove(idEpoch);
  removedEpochs -= arc.size();
  if(removedEpochs > 0)
    logWarning << "removed " << removedEpochs << " epochs outside interval boundaries from <" << fileName << ">" << Log::endl;
}

/***********************************************/
