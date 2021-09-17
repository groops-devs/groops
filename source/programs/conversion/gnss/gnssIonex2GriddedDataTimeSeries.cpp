/***********************************************/
/**
* @file gnssIonex2GriddedDataTimeSeries.cpp
*
* @brief Converts TEC maps from IGS IONEX file format to GROOPS GriddedDataTimeSeries format.
**
* @author Sebastian Strasser
* @date 2021-09-13
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts TEC maps from IGS \href{https://files.igs.org/pub/data/format/ionex1.pdf}{IONEX file format}
to GROOPS \file{gridded data time series}{griddedDataTimeSeries} format.

Currently only supports 2D TEC maps.

See also \program{GnssGriddedDataTimeSeries2Ionex}, \configClass{IonosphereMap}{gnssParametrizationType:ionosphereMap}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileGriddedDataTimeSeries.h"

/***** CLASS ***********************************/

/** @brief Converts TEC maps from IGS IONEX file format to GROOPS GriddedDataTimeSeries format.
* @ingroup programsGroup */
class GnssIonex2GriddedDataTimeSeries
{
  Bool getLine(InFile &file, std::string &line, std::string &label) const;
  Bool testLabel(const std::string &labelInLine, const std::string &label, Bool optional=TRUE) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssIonex2GriddedDataTimeSeries, SINGLEPROCESS, "Converts TEC maps from IGS IONEX file format to GROOPS GriddedDataTimeSeries format.", Conversion, Gnss, Grid, TimeSplines)

/***********************************************/

void GnssIonex2GriddedDataTimeSeries::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut;
    std::vector<FileName> fileNamesIn;

    readConfig(config, "outputfileGriddedDataTimeSeries", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileIonex",                  fileNamesIn, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    GriddedData grid;
    std::vector<Matrix> data;
    std::vector<Time> times;

    for(const auto &fileName : fileNamesIn)
    {
      logStatus<<"read IONEX file <"<<fileName<<">"<<Log::endl;
      InFile file(fileName);
      file.exceptions(std::ios::badbit|std::ios::failbit);

      Double factor = 0.1;
      std::vector<Angle> latitudes, longitudes;
      Double height = 0;
      Double radius = DEFAULT_GRS80_a;

      // read header
      std::string line, label;
      getLine(file, line, label);
      testLabel(label, "IONEX VERSION / TYPE", FALSE);
      while(getLine(file, line, label))
      {
        if(testLabel(label, "END OF HEADER"))
          break;
        else if(testLabel(label, "MAP DIMENSION"))
        {
          if(String::toInt(line.substr(0,6)) > 2)
            throw(Exception("only 2-dimensional maps are supported yet"));
        }
        else if(testLabel(label, "LAT1 / LAT2 / DLAT"))
        {
          Double lat1 = String::toDouble(line.substr(2,6));
          Double lat2 = String::toDouble(line.substr(8,6));
          Double dlat = String::toDouble(line.substr(14,6));
          for(Double lat=lat1; (dlat > 0 ? lat<=lat2 : lat>=lat2); lat+=dlat)
            latitudes.push_back(Angle(lat*DEG2RAD));
        }
        else if(testLabel(label, "LON1 / LON2 / DLON"))
        {
          Double lon1 = String::toDouble(line.substr(2,6));
          Double lon2 = String::toDouble(line.substr(8,6));
          Double dlon = String::toDouble(line.substr(14,6));
          for(Double lon=lon1; (dlon > 0 ? lon<=lon2 : lon>=lon2); lon+=dlon)
            longitudes.push_back(Angle(lon*DEG2RAD));
        }
        else if(testLabel(label, "HGT1 / HGT2 / DHGT"))
        {
          height = 1e3*String::toDouble(line.substr(2,6));
        }
        else if(testLabel(label, "BASE RADIUS"))
        {
          radius = 1e3*String::toDouble(line.substr(2,6));
        }
        else if(testLabel(label, "EXPONENT"))
        {
          factor = std::pow(10., String::toInt(line.substr(0,6)));
        }
      }

      if(!grid.points.size())
      {
        GriddedDataRectangular gridRectangular;
        gridRectangular.longitudes = longitudes;
        gridRectangular.latitudes  = latitudes;
        gridRectangular.heights.resize(latitudes.size(), height);
        gridRectangular.ellipsoid = Ellipsoid(radius);
        gridRectangular.convert(grid);
      }

      // read data
      while(getLine(file, line, label))
        if(testLabel(label, "START OF TEC MAP"))
        {
          getLine(file, line, label);
          testLabel(label, "EPOCH OF CURRENT MAP", FALSE);
          Time time = date2time(String::toInt(line.substr(0,6)), String::toInt(line.substr(6,6)), String::toInt(line.substr(12,6)),
                                String::toInt(line.substr(18,6)), String::toInt(line.substr(24,6)), String::toInt(line.substr(30,6)));

          if(times.size() && time < times.back())
          {
            logWarning << "skipping epoch " << time.dateTimeStr() << " (unsorted)" << Log::endl;
            continue;
          }

          // replace possible duplicate epoch by newer one (e.g. 00:00 from file 2 instead of 24:00 from file 1)
          if(times.size() && time == times.back())
            data.pop_back();
          else
            times.push_back(time);

          Vector epochData(latitudes.size()*longitudes.size());
          while(getLine(file, line, label))
          {
            if(testLabel(label, "END OF TEC MAP"))
              break;

            testLabel(label, "LAT/LON1/LON2/DLON/H", FALSE);
            auto iter = std::find(latitudes.begin(), latitudes.end(), Angle(String::toDouble(line.substr(2,6))*DEG2RAD));
            UInt index = std::distance(latitudes.begin(), iter);
            for(UInt i = 0; i < std::ceil(longitudes.size()/16.); i++)
            {
              getLine(file, line, label);
              for(UInt j = 0; j < std::min(longitudes.size()-i*16, UInt(16)); j++)
                epochData.at(index*longitudes.size() + i*16 + j) = String::toDouble(line.substr(j*5,5)) * factor;
            }
          }
          data.push_back(epochData);
        }
    }

    logStatus<<"write gridded data time series to <"<<fileNameOut<<">"<<Log::endl;
    writeFileGriddedDataTimeSeries(fileNameOut, 1, times, grid, data);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssIonex2GriddedDataTimeSeries::getLine(InFile &file, std::string &line, std::string &label) const
{
  try
  {
    std::getline(file, line);
    if(line.size() && line.back() == '\r')
      line.pop_back();
    if(line.size()<80)
      line.resize(80,' ');
    label = line.substr(60,20);
    return TRUE;
  }
  catch(...)
  {
    line.clear();
    line.resize(80,' ');
    label = line.substr(60,20);
    return FALSE;
  }
}

/***********************************************/

Bool GnssIonex2GriddedDataTimeSeries::testLabel(const std::string &labelInLine, const std::string &label, Bool optional) const
{
  if(labelInLine.find(label)!=std::string::npos)
    return TRUE;
  if(optional)
    return FALSE;
  throw(Exception(std::string("In line '")+labelInLine+"' label '"+label+"' expected\n"));
}

/***********************************************/
