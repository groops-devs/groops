/***********************************************/
/**
* @file psmslOceanBottomPressure2TimeSeries.cpp
*
* @brief Read ocean bottom pressure (OBP) time series in the PSMSL file format.
*
* @author Andreas Kvas
* @date 2019-12-10
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This programs reads ocean bottom pressure time series from the Permanent Service for Mean Sea Level (PSMSL).

In addition to the OBP measurements, the recorder position can be written to a \file{grid file}{griddedData}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "files/fileGriddedData.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Read ocean bottom pressure (OBP) time series in the PSMSL file format
* @ingroup programsConversionGroup */
class PsmslOceanBottomPressure2TimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(PsmslOceanBottomPressure2TimeSeries, SINGLEPROCESS, "read ocean bottom pressure (OBP) time series in the PSMSL file format", Conversion)
GROOPS_RENAMED_PROGRAM(ReadOceanBottomPressurePSMSL, PsmslOceanBottomPressure2TimeSeries, date2time(2020, 9, 10))

/***********************************************/

void PsmslOceanBottomPressure2TimeSeries::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName      outName, posName, timeName;
    FileName      inName;
    Double        a, f;
    TimeSeriesPtr timeSeries;
    Bool          isDaily, ignoreBadData;

    readConfig(config, "outputfileTimeSeries", outName,       Config::MUSTSET,  "",  "");
    readConfig(config, "outputfilePosition",   posName,       Config::OPTIONAL, "",  "recorder position as gridded data");
    readConfig(config, "inputfile",            inName,        Config::MUSTSET,  "",  "");
    readConfig(config, "isDaily",              isDaily,       Config::MUSTSET,  "0", "");
    readConfig(config, "ignoreBadData",        ignoreBadData, Config::MUSTSET,  "1", "");
    readConfig(config, "R",                    a,             Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "");
    readConfig(config, "inverseFlattening",    f,             Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "");
    readConfig(config, "timeSeries",           timeSeries,    Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    std::vector<Time> times = timeSeries->times();

    logStatus<<"read inputfile <"<<inName<<">"<<Log::endl;
    InFile file(inName);

    // Header einlesen
    Double L=0, B=0;
    std::string line;
    while(std::getline(file, line))
    {
      std::stringstream ss(line);
      std::string tag;
      std::string tmp;
      ss>>tag;
      if(tag == "Latitude")
        ss>>tmp>>tmp>>B;
      if(tag == "Longitude")
        ss>>tmp>>tmp>>L;
      if(tag == "Recno")
        break;
    }

    std::vector<Time>   times2;
    std::vector<Double> values2;
    while(std::getline(file, line))
    {
      std::stringstream ss(line);
      std::string tag;
      std::string tmp;
      UInt year, day, flag;
      Double hour;
      ss>>tmp;
      ss>>flag;
      if(ignoreBadData && flag==1) continue;
      ss>>year>>day;

      Time time = date2time(year, 1, 0);
      time += seconds2time(static_cast<Double>(day)*86400.0);
      if(!isDaily)
      {
        ss>>hour;
        time += seconds2time(hour*60.0*60.0);
      }
      Double total, tide, drift, resid;
      if(!isDaily) ss>>total>>tide;
      ss>>drift>>resid;
      if(!ss.good())
        continue;

      times2.push_back(time);
      values2.push_back(resid);
    }

    Double mean = 0.0;
    for(UInt i=0; i<values2.size(); i++)
      mean += values2.at(i);
    mean /= values2.size();
    for(UInt i=0; i<values2.size(); i++)
      values2.at(i) -= mean;

    UInt idx = 0;
    std::vector<UInt>   count(times.size()-1, 0);
    std::vector<Double> value(times.size()-1, 0.0);
    for(UInt i=0; i<times2.size(); i++)
    {
      if(times2.at(i) < times.at(idx))
        continue;
      while((idx<times.size()-1) && (times2.at(i)>=times.at(idx+1))) idx++;
      if((idx>=times.size()-1)||(times2.at(i)>=times.at(times.size()-1)))
        break;

      count.at(idx)++;
      value.at(idx) += values2.at(i);
    }

    logStatus<<"remove mean"<<Log::endl;
    mean = 0.0;
    UInt meanCount = 0;
    for(UInt i=0; i<value.size(); i++)
      if(count.at(i) != 0)
      {
        value.at(i) /= count.at(i);
        mean += value.at(i);
        meanCount++;
      }
    mean /= meanCount;
    meanCount = 0;
    for(UInt i=0; i<value.size(); i++)
      if(count.at(i)!=0)
      {
        value.at(i) -= mean;
        meanCount++;
      }
      else
        count.at(i)=0;

    logInfo<<"  mean  = "<<mean<<Log::endl;
    logInfo<<"  count = "<<meanCount<<Log::endl;
    logInfo<<"  L     = "<<L<<Log::endl;
    logInfo<<"  B     = "<<B<<Log::endl;

    if(meanCount == 0)
      return;

    if(!outName.empty())
    {
      logStatus<<"write time series to file <"<<outName<<">"<<Log::endl;
      Matrix A(meanCount,2);
      UInt idx = 0;
      for(UInt i=0; i<value.size(); i++)
        if(count.at(i) != 0)
        {
          A(idx,0) = times.at(i).mjd();//(0.5*times.at(i)+0.5*times.at(i+1)).mjd();
          A(idx,1) = value.at(i);
          idx++;
        }
      InstrumentFile::write(outName, Arc(A));
    }

    if(!posName.empty())
    {
      logStatus<<"write position to file <"<<posName<<">"<<Log::endl;
      Ellipsoid ellipsoid(a,f);
      writeFileGriddedData(posName, GriddedData(ellipsoid, {ellipsoid(Angle(L*DEG2RAD), Angle(B*DEG2RAD), 0.0)}, {1.}, {}));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
