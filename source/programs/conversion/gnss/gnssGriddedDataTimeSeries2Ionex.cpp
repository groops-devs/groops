/***********************************************/
/**
* @file programTemplate.cpp
*
* @brief Converts TEC maps from GROOPS GriddedDataTimeSeries format to IGS IONEX file format.
*
* @author Sebastian Strasser
* @date 2021-09-15
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts TEC maps from GROOPS \file{gridded data time series}{griddedDataTimeSeries} format
to IGS \href{https://files.igs.org/pub/data/format/ionex1.pdf}{IONEX file format}.

Currently only supports 2D TEC maps.

See also \program{GnssIonex2GriddedDataTimeSeries}, \configClass{IonosphereMap}{gnssParametrizationType:ionosphereMap}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/system.h"
#include "classes/timeSeries/timeSeries.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "parser/dataVariables.h"

/***** CLASS ***********************************/

/** @brief Converts TEC maps from GROOPS GriddedDataTimeSeries format to IGS IONEX file format.
* @ingroup programsGroup */
class GnssGriddedDataTimeSeries2Ionex
{
  static std::string resize(std::string str, UInt length) {str.resize(length, ' '); return str;}

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssGriddedDataTimeSeries2Ionex, SINGLEPROCESS, "Converts TEC maps from GROOPS GriddedDataTimeSeries format to IGS IONEX file format.", Conversion, Gnss, Grid, TimeSplines)

/***********************************************/

void GnssGriddedDataTimeSeries2Ionex::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameIn, fileNameOut;
    TimeSeriesPtr timeSeries;
    std::string program, institution, mappingFunction, observablesUsed;
    std::vector<std::string> comments, descriptions;
    Double elevationCutoff;
    Int exponent;
    ExpressionVariablePtr expression;

    readConfig(config, "outputfileIonex",                fileNameOut,     Config::MUSTSET,  "",       "");
    readConfig(config, "inputfileGriddedDataTimeSeries", fileNameIn,      Config::MUSTSET,  "",       "must contain regular grid");
    readConfig(config, "value",                          expression,      Config::MUSTSET,  "data0",  "expression (e.g. data column)");
    readConfig(config, "timeSeries",                     timeSeries,      Config::OPTIONAL, "",       "(empty = use input file time series)");
    readConfig(config, "program",                        program,         Config::MUSTSET,  "GROOPS", "name of program (for first line)");
    readConfig(config, "institution",                    institution,     Config::MUSTSET,  "TUG (TU Graz)", "name of agency (for first line)");
    readConfig(config, "description",                    descriptions,    Config::OPTIONAL, "",       "description in header");
    readConfig(config, "comment",                        comments,        Config::OPTIONAL, R"(["TEC values in 0.1 TECU, 9999 if no value available"])", "comment in header");
    readConfig(config, "mappingFunction",                mappingFunction, Config::OPTIONAL, "MSLM",   "see IONEX documentation");
    readConfig(config, "elevationCutoff",                elevationCutoff, Config::DEFAULT,  "0",      "see IONEX documentation (0 if unknown)");
    readConfig(config, "observablesUsed",                observablesUsed, Config::OPTIONAL, "",       "see IONEX documentation");
    readConfig(config, "exponent",                       exponent,        Config::DEFAULT,  "-1",     "factor 10^exponent is applied to all values");
    if(isCreateSchema(config)) return;

    logStatus<<"read gridded data time series file <"<<fileNameIn<<">"<<Log::endl;
    InFileGriddedDataTimeSeries fileGrid(fileNameIn);
    GriddedDataRectangular grid;
    if(!grid.init(fileGrid.grid()))
      throw(Exception("grid must be rectangular"));

    std::vector<Time> times = fileGrid.times();
    if(timeSeries)
      times = timeSeries->times();

    logStatus<<"write IONEX file <"<<fileNameOut<<">"<<Log::endl;
    OutFile file(fileNameOut);

    file<<"     1.0            IONOSPHERE MAPS     GNSS                IONEX VERSION / TYPE"<<std::endl;
    file<<resize(program, 20)<<resize(institution, 20)<<System::now()%"%D %H:%M:%S "s<<"PGM / RUN BY / DATE"<<std::endl;
    for(auto &description : descriptions)
      file<<resize(description, 20)<<"DESCRIPTION"<<std::endl;
    file<<times.front()%"  %y    %m    %d    %H    %M    %S"s<<std::string(24, ' ')<<"EPOCH OF FIRST MAP"<<std::endl;
    file<<times.back() %"  %y    %m    %d    %H    %M    %S"s<<std::string(24, ' ')<<"EPOCH OF LAST MAP"<<std::endl;
    file<<medianSampling(times).seconds()%"% 6i"s<<std::string(54, ' ')<<"INTERVAL"<<std::endl;
    file<<times.size()%"% 6i"s<<std::string(54, ' ')<<"# OF MAPS IN FILE"<<std::endl;
    if(!mappingFunction.empty()) file<<"  "<<resize(mappingFunction, 4)<<std::string(54, ' ')<<"MAPPING FUNCTION"<<std::endl;
    file<<elevationCutoff%"%8.1f"s<<std::string(52, ' ')<<"ELEVATION CUTOFF"<<std::endl;
    if(!observablesUsed.empty()) file<<resize(observablesUsed, 60)<<"OBSERVABLES USED"<<std::endl;
    file<<(fileGrid.grid().ellipsoid.a()/1000)%"%8.1f"s<<std::string(52, ' ')<<"BASE RADIUS"<<std::endl;
    file<<"     2"<<std::string(54, ' ')<<"MAP DIMENSION"<<std::endl;
    // Header: HGT1 / HGT 2 / DHGT, LAT1 / LAT2 / DLAT, LON1 / LON2 / DLON
    const Double height = mean(Vector(grid.heights));
    UInt columns = (grid.longitudes.front()+2*PI == grid.longitudes.back() ? grid.longitudes.size()-1 : grid.longitudes.size());
    file<<"  "<<(height/1000)%"%6.1f%6.1f"s<<0%"%6.1f"s<<std::string(40, ' ')<<"HGT1 / HGT2 / DHGT"<<std::endl;
    file<<"  "<<(grid.latitudes.front()*RAD2DEG)%"%6.1f"s<<(grid.latitudes.back()*RAD2DEG)%"%6.1f"s
        <<((grid.latitudes.back()-grid.latitudes.front())/grid.latitudes.size()*RAD2DEG)%"%6.1f"s<<std::string(40, ' ')<<"LAT1 / LAT2 / DLAT"<<std::endl;
    file<<"  "<<(grid.longitudes.front()*RAD2DEG)%"%6.1f"s<<(grid.longitudes.back()*RAD2DEG)%"%6.1f"s
        <<((grid.longitudes.back()-grid.longitudes.front())/columns*RAD2DEG)%"%6.1f"s<<std::string(40, ' ')<<"LON1 / LON2 / DLON"<<std::endl;
    file<<exponent%"% 6i"s<<std::string(54, ' ')<<"EXPONENT"<<std::endl;
    for(auto &comment : comments)
      file<<resize(comment ,60)<<"COMMENT"<<std::endl;
    file<<std::string(60, ' ')<<"END OF HEADER"<<std::endl;

    // Data records
    const Double factor = std::pow(10, exponent);
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      file<<(idEpoch+1)%"% 6i"s<<std::string(54, ' ')<<"START OF TEC MAP"<<std::endl;
      file<<times.at(idEpoch)%"  %y    %m    %d    %H    %M    %S"s<<std::string(24, ' ')<<"EPOCH OF CURRENT MAP"<<std::endl;
      for(UInt i=0; i<grid.latitudes.size(); i++)
      {
        file<<"  "<<(grid.latitudes.at(i)*RAD2DEG)%"% 6.1f"s<< (grid.longitudes.front()*RAD2DEG)%"%6.1f"s<<(grid.longitudes.back()*RAD2DEG)%"%6.1f"s
            <<((grid.longitudes.back()-grid.longitudes.front())/columns*RAD2DEG)%"%6.1f"s<<(height/1000)%"%6.1f"s<<std::string(28, ' ')<<"LAT/LON1/LON2/DLON/H"<<std::endl;

        Matrix data = fileGrid.data(times.at(idEpoch));
        VariableList varList;
        addDataVariables(data, varList);

        for(UInt k=0; k<grid.longitudes.size(); k++)
        {
          evaluateDataVariables(data, i*grid.longitudes.size()+k, varList);
          Double value = expression->evaluate(varList);
          if(k > 0 && k%16 == 0) // start new line after 16 values
            file<<std::endl;
          file<<(std::isnan(value) ? " 9999" : (value/factor)%"% 5i"s);
        }
        file<<std::endl;
      }
      file<<(idEpoch+1)%"% 6i"s<<std::string(54, ' ')<<"END OF TEC MAP"<<std::endl;
    }
    file<<std::string(60, ' ')<<"END OF FILE"<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
