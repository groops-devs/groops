/***********************************************/
/**
* @file icgem2PotentialCoefficients.cpp
*
* @brief Read spherical harmonics in ICGEM format.
*
* @author Andreas Kvas
* @date 2019-06-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read spherical harmonics in ICGEM format (\url{http://icgem.gfz-potsdam.de/}).
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileSphericalHarmonics.h"

/***** CLASS ***********************************/

/** @brief Write spherical harmonics in ICGEM format.
* @ingroup programsConversionGroup */
class Icgem2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);

  class Coefficient
  {
  public:
    Coefficient(UInt _n, UInt _m, Double _cnm, Double _snm, Double _cnm_e, Double _snm_e) { n = _n; m = _m; cnm=_cnm; snm=_snm; cnm_error=_cnm_e; snm_error=_snm_e; }

    enum Type {STATIC, STATIC_INTERVAL, TREND, OSC_COSINE, OSC_SINE };

    Double cnm, snm, cnm_error, snm_error;
    UInt n, m;

    Time t0, t1;
    Double period;

    Type coefficientType;
  };

  class Field
  {
  public:
    Field() {}

    Field(UInt maxDegree)
    {
      _cnm = Matrix(maxDegree+1, maxDegree+1);
      _snm = Matrix(maxDegree+1, maxDegree+1);
      _cnm_error = Matrix(maxDegree+1, maxDegree+1);
      _snm_error = Matrix(maxDegree+1, maxDegree+1);
    }

    Matrix _cnm, _snm, _cnm_error, _snm_error;
  };

};

GROOPS_REGISTER_PROGRAM(Icgem2PotentialCoefficients, SINGLEPROCESS, "read spherical harmonics in ICGEM format", Conversion, PotentialCoefficients)

/***********************************************/

void Icgem2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameStatic, fileNameTrend, fileNameOscCos, fileNameOscSin, fileNameIntervals;
    FileName fileNameIn;
    Bool     useFormalErrors;

    readConfig(config, "outputfileStaticCoefficients", fileNameStatic,    Config::MUSTSET,  "",  "static potential coefficients in GROOPS gfc format. Available variables (icgem2.0): epochStart, epochEnd, epochMid; (icgem1.0) epochReference");
    readConfig(config, "outputfileTrendCoefficients",  fileNameTrend,     Config::OPTIONAL, "",  "trend potential coefficients in GROOPS gfc format.  Available variables (icgem2.0): epochStart, epochEnd, epochMid; (icgem1.0) epochReference");
    readConfig(config, "outputfileOscillationCosine",  fileNameOscCos,    Config::OPTIONAL, "",  "oscillation cosine coefficients in GROOPS gfc format. Available variables (icgem2.0): epochStart, epochEnd, epochMid, oscillationPeriod; (icgem1.0) epochReference, oscillationPeriod");
    readConfig(config, "outputfileOscillationSine",    fileNameOscSin,    Config::OPTIONAL, "",  "oscillation sine coefficients in GROOPS gfc format. Available variables (icgem2.0): epochStart, epochEnd, epochMid, oscillationPeriod; (icgem1.0) epochReference, oscillationPeriod");
    readConfig(config, "outputfileIntervals",          fileNameIntervals, Config::OPTIONAL, "",  "two column ASCII file with all intervals found (only sensible for icgem2.0). The base name will be extended with .static, .trend, .annualCos, and .annualSin.");
    readConfig(config, "inputfileIcgem",               fileNameIn,        Config::MUSTSET,  "",  "ICGEM GFC file");
    readConfig(config, "useFormalErrors",              useFormalErrors,   Config::DEFAULT,  "0", "use formal errors if both formal and calibrated errors are given");
    if(isCreateSchema(config)) return;

    // convenience functions
    // ---------------------
    auto parseTimeStamp = [](const std::string &timeStamp)
    {
      UInt year  = std::stoul(timeStamp.substr(0, 4));
      UInt month = std::stoul(timeStamp.substr(4, 2));
      UInt day   = std::stoul(timeStamp.substr(6, 2));
      Double fraction = 0.0;
      if(timeStamp.find('.') != std::string::npos)
        fraction = std::stod(timeStamp.substr(8, timeStamp.size()-8));
      return date2time(year, month, day) + mjd2time(fraction);
    };

    auto splitLine = [](const std::string &line)
    {
      std::vector<std::string> tokens;
      std::stringstream ss(line);
      std::string token;
      while(ss >> token)
        tokens.push_back(token);
      return tokens;
    };

    // read header
    // -----------
    logStatus<<"read header of file <"<<fileNameIn<<">"<<Log::endl;
    Double GM, R;
    Bool hasFormalError = FALSE;
    Bool hasCalibratedError = FALSE;
    Bool isVersion2 = FALSE;

    InFile inputFile(fileNameIn);
    std::string line;
    while(!inputFile.eof())
    {
      std::getline(inputFile, line);
      if(String::startsWith(line, "earth_gravity_constant"))
      {
        std::stringstream ss(line);
        std::string dummy; ss>>dummy; ss>>GM;
      }
      if(String::startsWith(line, "radius"))
      {
        std::stringstream ss(line);
        std::string dummy; ss>>dummy; ss>>R;
      }
      if(String::startsWith(line, "format"))
      {
        std::stringstream ss(line);
        std::string dummy; ss>>dummy; ss>>dummy;
        if(dummy == "icgem2.0")
          isVersion2 = TRUE;
      }
      if(String::startsWith(line, "errors"))
      {
        std::stringstream ss(line);
        std::string dummy; ss>>dummy; ss>>dummy;
        if(dummy == "calibrated")
          hasCalibratedError = TRUE;

        if(dummy == "formal")
          hasFormalError = TRUE;

        if(dummy == "calibrated_and_formal")
        {
          hasCalibratedError = TRUE;
          hasFormalError = TRUE;
        }

      }
      if(String::startsWith(line, "end_of_head"))
        break;
    }

    // read data
    // ---------
    logStatus<<"read potential coeffcients"<<Log::endl;
    std::vector<Coefficient> coefficients;
    UInt maxDegree = 0;
    while(!inputFile.eof())
    {
      std::getline(inputFile, line);
      if(line.size() == 0)
        continue;
      std::vector<std::string> tokens = splitLine(line);
      if(tokens.size()<5)
        continue;
      UInt offset = 1;
      UInt n = std::stoul(tokens.at(offset++));
      maxDegree = std::max(n, maxDegree);
      UInt m = std::stoul(tokens.at(offset++));

      Double cnm = std::stod(tokens.at(offset++));
      Double snm = std::stod(tokens.at(offset++));

      Double cnm_error = 0.0;
      Double snm_error = 0.0;

      if(hasFormalError && hasCalibratedError)
      {
        if(useFormalErrors) offset += 2;
        cnm_error = std::stod(tokens.at(offset++));
        snm_error = std::stod(tokens.at(offset++));
      }
      if(hasFormalError || hasCalibratedError)
      {
        cnm_error = std::stod(tokens.at(offset++));
        snm_error = std::stod(tokens.at(offset++));
      }
      Coefficient c(n, m, cnm, snm, cnm_error*cnm_error, snm_error*snm_error);

      Time versionOneRefTime;
      if(tokens.front() == "gfc")
      {
        c.coefficientType = Coefficient::STATIC;
      }
      else if(tokens.front() == "gfct")
      {
        c.coefficientType = Coefficient::STATIC_INTERVAL;
        c.t0 = parseTimeStamp(tokens.at(offset++));
        versionOneRefTime = c.t0;
        if(isVersion2)
          c.t1 = parseTimeStamp(tokens.at(offset++));
      }
      else if(tokens.front() == "trnd")
      {
        c.coefficientType = Coefficient::TREND;
        if(isVersion2)
        {
          c.t0 = parseTimeStamp(tokens.at(offset++));
          c.t1 = parseTimeStamp(tokens.at(offset++));
        }
      }
      else if(tokens.front() == "acos")
      {
        c.coefficientType = Coefficient::OSC_COSINE;
        if(isVersion2)
        {
          c.t0 = parseTimeStamp(tokens.at(offset++));
          c.t1 = parseTimeStamp(tokens.at(offset++));
        }
        c.period = std::stod(tokens.at(offset));
      }
      else if(tokens.front() == "asin")
      {
        c.coefficientType = Coefficient::OSC_SINE;
        if(isVersion2)
        {
          c.t0 = parseTimeStamp(tokens.at(offset++));
          c.t1 = parseTimeStamp(tokens.at(offset++));
        }
        c.period = std::stod(tokens.at(offset));
      }
      if(!isVersion2)
      {
        c.t0 = versionOneRefTime;
        c.t1 = versionOneRefTime;
      }
      coefficients.push_back(c);
    }
    // sort and reorder data
    // ---------------------
    logStatus<<"split coefficients into static and time-variable parts"<<Log::endl;
    if(!isVersion2)
    {
      Time refTime;
      VariableList fileNameVariableList;
      addVariable("epochReference", fileNameVariableList);
      addVariable("oscillationPeriod", fileNameVariableList);
      {
        Field staticField(maxDegree);
        for(const Coefficient &c : coefficients)
        {
          if( (c.coefficientType == Coefficient::STATIC) || (c.coefficientType == Coefficient::STATIC_INTERVAL) )
          {
            staticField._cnm(c.n, c.m) = c.cnm;
            staticField._snm(c.n, c.m) = c.snm;
            staticField._cnm_error(c.n, c.m) = c.cnm_error;
            staticField._snm_error(c.n, c.m) = c.snm_error;
            if(c.coefficientType == Coefficient::STATIC_INTERVAL)
              refTime = c.t0;
          }
        }
        SphericalHarmonics harm(GM, R, staticField._cnm, staticField._snm, staticField._cnm_error, staticField._snm_error);
        fileNameVariableList["epochReference"]->setValue(refTime.mjd());
        logStatus<<"write static potential coefficients to <"<<fileNameStatic(fileNameVariableList)<<">"<<Log::endl;
        writeFileSphericalHarmonics(fileNameStatic(fileNameVariableList), harm);
      }

      if(!fileNameTrend.empty())
      {
        Field trendField(maxDegree);
        UInt maxDegreeTemp = 0;
        for(const Coefficient &c : coefficients)
          if( (c.coefficientType == Coefficient::TREND) )
          {
            maxDegreeTemp = std::max(maxDegreeTemp, c.n);
            trendField._cnm(c.n, c.m) = c.cnm;
            trendField._snm(c.n, c.m) = c.snm;
            trendField._cnm_error(c.n, c.m) = c.cnm_error;
            trendField._snm_error(c.n, c.m) = c.snm_error;
          }
        SphericalHarmonics harm(GM, R, trendField._cnm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
        trendField._snm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
        trendField._cnm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
        trendField._snm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1));
        logStatus<<"write trend potential coefficients to <"<<fileNameTrend(fileNameVariableList)<<">"<<Log::endl;
        writeFileSphericalHarmonics(fileNameTrend(fileNameVariableList), harm);
      }

      if(!fileNameOscCos.empty())
      {
        std::map<Double, Field> oscillationMap;
        UInt maxDegreeTemp = 0;
        for(const Coefficient &c : coefficients)
          if( (c.coefficientType == Coefficient::OSC_COSINE) )
          {
            maxDegreeTemp = std::max(maxDegreeTemp, c.n);
            if(oscillationMap.find(c.period) == oscillationMap.end())
              oscillationMap[c.period] = Field(maxDegree);

            oscillationMap[c.period]._cnm(c.n, c.m) = c.cnm;
            oscillationMap[c.period]._snm(c.n, c.m) = c.snm;
            oscillationMap[c.period]._cnm_error(c.n, c.m) = c.cnm_error;
            oscillationMap[c.period]._snm_error(c.n, c.m) = c.snm_error;
          }

        for(auto &entry : oscillationMap)
        {
          Double period = entry.first;

          SphericalHarmonics harm(GM, R, entry.second._cnm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._cnm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1));
          fileNameVariableList["oscillationPeriod"]->setValue(period);
          logStatus<<"write cosine potential coefficients to <"<<fileNameOscCos(fileNameVariableList)<<">"<<Log::endl;
          writeFileSphericalHarmonics(fileNameOscCos(fileNameVariableList), harm);
        }
      }
      if(!fileNameOscCos.empty())
      {
        std::map<Double, Field> oscillationMap;
        UInt maxDegreeTemp = 0;
        for(const Coefficient &c : coefficients)
          if( (c.coefficientType == Coefficient::OSC_SINE) )
          {
            maxDegreeTemp = std::max(maxDegreeTemp, c.n);
            if(oscillationMap.find(c.period) == oscillationMap.end())
              oscillationMap[c.period] = Field(maxDegree);

            oscillationMap[c.period]._cnm(c.n, c.m) = c.cnm;
            oscillationMap[c.period]._snm(c.n, c.m) = c.snm;
            oscillationMap[c.period]._cnm_error(c.n, c.m) = c.cnm_error;
            oscillationMap[c.period]._snm_error(c.n, c.m) = c.snm_error;
          }

        for(auto &entry : oscillationMap)
        {
          Double period = entry.first;

          SphericalHarmonics harm(GM, R, entry.second._cnm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._cnm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1));
          fileNameVariableList["oscillationPeriod"]->setValue(period);
          logStatus<<"write sine potential coefficients to <"<<fileNameOscSin(fileNameVariableList)<<">"<<Log::endl;
          writeFileSphericalHarmonics(fileNameOscSin(fileNameVariableList), harm);
        }
      }

    }
    else // icgem2.0
    {
      VariableList fileNameVariableList;
      addVariable("epochStart", fileNameVariableList);
      addVariable("epochEnd", fileNameVariableList);
      addVariable("epochMid", fileNameVariableList);
      addVariable("oscillationPeriod", fileNameVariableList);

      Field staticGLobal(maxDegree);
      for(const Coefficient &c : coefficients)
        if(c.coefficientType == Coefficient::STATIC)
        {
          staticGLobal._cnm(c.n, c.m) = c.cnm;
          staticGLobal._snm(c.n, c.m) = c.snm;
          staticGLobal._cnm_error(c.n, c.m) = c.cnm_error;
          staticGLobal._snm_error(c.n, c.m) = c.snm_error;
        }

      {
        std::map< std::pair<Time,Time>, Field> staticInterval;
        UInt maxDegreeTemp = maxDegree;
        for(const Coefficient &c : coefficients)
        {
          if(c.coefficientType == Coefficient::STATIC_INTERVAL)
          {
            std::pair<Time, Time> interval(c.t0, c.t1);
            if(staticInterval.find(interval) == staticInterval.end())
              staticInterval[interval] = Field(staticGLobal);

            staticInterval[interval]._cnm(c.n, c.m) = c.cnm;
            staticInterval[interval]._snm(c.n, c.m) = c.snm;
            staticInterval[interval]._cnm_error(c.n, c.m) = c.cnm_error;
            staticInterval[interval]._snm_error(c.n, c.m) = c.snm_error;

            maxDegreeTemp = std::max(maxDegreeTemp, c.n);
          }
        }
        if(staticInterval.size() == 0)
        {
          SphericalHarmonics harm(GM, R, staticGLobal._cnm, staticGLobal._snm, staticGLobal._cnm_error, staticGLobal._snm_error);
          logStatus<<"write static potential coefficients to <"<<fileNameStatic<<">"<<Log::endl;
          writeFileSphericalHarmonics(fileNameStatic, harm);
        }
        std::vector<std::pair<Double,Double>> intervals;
        for(auto &entry : staticInterval)
        {
          Time timeStart = entry.first.first;
          Time timeEnd = entry.first.second;
          fileNameVariableList["epochStart"]-> setValue(timeStart.mjd());
          fileNameVariableList["epochMid"]-> setValue(timeStart.mjd()*0.5 + timeEnd.mjd()*0.5);
          fileNameVariableList["epochEnd"]-> setValue(timeEnd.mjd());

          intervals.push_back(std::pair<Double,Double>(timeStart.mjd(), timeEnd.mjd()));

           SphericalHarmonics harm(GM, R, entry.second._cnm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._cnm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1));
          logStatus<<"write static potential coefficients for interval ("<<timeStart.dateStr()<<", "<<timeEnd.dateStr()<<") to <"<<fileNameStatic(fileNameVariableList)<<">"<<Log::endl;
          writeFileSphericalHarmonics(fileNameStatic(fileNameVariableList), harm);
        }
        if(!fileNameIntervals.empty())
        {
          OutFile intervalFile(fileNameIntervals.appendBaseName(".static"));
          for(auto &i : intervals)
            intervalFile<<i.first<<" "<<i.second<<std::endl;
        }
      }

      if(!fileNameTrend.empty())
      {
        std::map< std::pair<Time,Time>, Field> trendMap;
        UInt maxDegreeTemp = 0;
        for(const Coefficient &c : coefficients)
        {
          if(c.coefficientType == Coefficient::TREND)
          {
            std::pair<Time, Time> interval(c.t0, c.t1);
            if(trendMap.find(interval) == trendMap.end())
              trendMap[interval] = Field(maxDegree);

            trendMap[interval]._cnm(c.n, c.m) = c.cnm;
            trendMap[interval]._snm(c.n, c.m) = c.snm;
            trendMap[interval]._cnm_error(c.n, c.m) = c.cnm_error;
            trendMap[interval]._snm_error(c.n, c.m) = c.snm_error;

            maxDegreeTemp = std::max(maxDegreeTemp, c.n);
          }
        }
        std::vector<std::pair<Double,Double>> intervals;
        for(auto &entry : trendMap)
        {
          Time timeStart = entry.first.first;
          Time timeEnd = entry.first.second;
          fileNameVariableList["epochStart"]-> setValue(timeStart.mjd());
          fileNameVariableList["epochMid"]-> setValue(timeStart.mjd()*0.5 + timeEnd.mjd()*0.5);
          fileNameVariableList["epochEnd"]-> setValue(timeEnd.mjd());

          intervals.push_back(std::pair<Double,Double>(timeStart.mjd(), timeEnd.mjd()));

          SphericalHarmonics harm(GM, R, entry.second._cnm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._cnm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1));
          logStatus<<"write trend potential coefficients for interval ("<<timeStart.dateStr()<<", "<<timeEnd.dateStr()<<") to <"<<fileNameTrend(fileNameVariableList)<<">"<<Log::endl;
          writeFileSphericalHarmonics(fileNameTrend(fileNameVariableList), harm);
        }
        if(!fileNameIntervals.empty())
        {
          OutFile intervalFile(fileNameIntervals.appendBaseName(".trend"));
          for(auto &i : intervals)
            intervalFile<<i.first<<" "<<i.second<<std::endl;
        }
      }

      if(!fileNameOscCos.empty())
      {
        std::map< std::tuple<Time,Time,Double>, Field> oscMap;
        UInt maxDegreeTemp = 0;
        for(const Coefficient &c : coefficients)
        {
          if(c.coefficientType == Coefficient::OSC_COSINE)
          {
            std::tuple<Time, Time, Double> interval(c.t0, c.t1, c.period);
            if(oscMap.find(interval) == oscMap.end())
              oscMap[interval] = Field(maxDegree);

            oscMap[interval]._cnm(c.n, c.m) = c.cnm;
            oscMap[interval]._snm(c.n, c.m) = c.snm;
            oscMap[interval]._cnm_error(c.n, c.m) = c.cnm_error;
            oscMap[interval]._snm_error(c.n, c.m) = c.snm_error;
            maxDegreeTemp = std::max(maxDegreeTemp, c.n);
          }
        }
        std::vector<std::pair<Double,Double>> intervals;
        for(auto &entry : oscMap)
        {
          Time timeStart = std::get<0>(entry.first);
          Time timeEnd = std::get<1>(entry.first);
          Double period = std::get<2>(entry.first);
          fileNameVariableList["epochStart"]-> setValue(timeStart.mjd());
          fileNameVariableList["epochMid"]-> setValue(timeStart.mjd()*0.5 + timeEnd.mjd()*0.5);
          fileNameVariableList["epochEnd"]-> setValue(timeEnd.mjd());
          fileNameVariableList["oscillationPeriod"]->setValue(period);

          intervals.push_back(std::pair<Double,Double>(timeStart.mjd(), timeEnd.mjd()));

          SphericalHarmonics harm(GM, R, entry.second._cnm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._cnm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1));
          logStatus<<"write cosine potential coefficients for period "<<period<<" in interval ("<<timeStart.dateStr()<<", "<<timeEnd.dateStr()<<") to <"<<fileNameOscCos(fileNameVariableList)<<">"<<Log::endl;
          writeFileSphericalHarmonics(fileNameOscCos(fileNameVariableList), harm);
        }
        if(!fileNameIntervals.empty())
        {
          OutFile intervalFile(fileNameIntervals.appendBaseName(".annualCos"));
          for(auto &i : intervals)
            intervalFile<<i.first<<" "<<i.second<<std::endl;
        }
      }
      if(!fileNameOscSin.empty())
      {
        std::map< std::tuple<Time,Time,Double>, Field> oscMap;
        UInt maxDegreeTemp = 0;
        for(const Coefficient &c : coefficients)
        {
          if(c.coefficientType == Coefficient::OSC_SINE)
          {
            std::tuple<Time, Time, Double> interval(c.t0, c.t1, c.period);
            if(oscMap.find(interval) == oscMap.end())
              oscMap[interval] = Field(maxDegree);

            oscMap[interval]._cnm(c.n, c.m) = c.cnm;
            oscMap[interval]._snm(c.n, c.m) = c.snm;
            oscMap[interval]._cnm_error(c.n, c.m) = c.cnm_error;
            oscMap[interval]._snm_error(c.n, c.m) = c.snm_error;
            maxDegreeTemp = std::max(maxDegreeTemp, c.n);
          }
        }
        std::vector<std::pair<Double,Double>> intervals;
        for(auto &entry : oscMap)
        {
          Time timeStart = std::get<0>(entry.first);
          Time timeEnd = std::get<1>(entry.first);
          Double period = std::get<2>(entry.first);
          fileNameVariableList["epochStart"]-> setValue(timeStart.mjd());
          fileNameVariableList["epochMid"]-> setValue(timeStart.mjd()*0.5 + timeEnd.mjd()*0.5);
          fileNameVariableList["epochEnd"]-> setValue(timeEnd.mjd());
          fileNameVariableList["oscillationPeriod"]->setValue(period);

          intervals.push_back(std::pair<Double,Double>(timeStart.mjd(), timeEnd.mjd()));

          SphericalHarmonics harm(GM, R, entry.second._cnm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._cnm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1),
          entry.second._snm_error.slice(0, 0, maxDegreeTemp+1, maxDegreeTemp+1));
          logStatus<<"write sine potential coefficients for period"<<period<<" in interval ("<<timeStart.dateStr()<<", "<<timeEnd.dateStr()<<") to <"<<fileNameOscSin(fileNameVariableList)<<">"<<Log::endl;
          writeFileSphericalHarmonics(fileNameOscSin(fileNameVariableList), harm);
        }
        if(!fileNameIntervals.empty())
        {
          OutFile intervalFile(fileNameIntervals.appendBaseName(".annualSin"));
          for(auto &i : intervals)
            intervalFile<<i.first<<" "<<i.second<<std::endl;
        }
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
