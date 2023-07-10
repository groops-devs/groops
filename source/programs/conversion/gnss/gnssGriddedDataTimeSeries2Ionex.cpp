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
    std::vector<std::string> comment, description;
    Double elevationCutoff;
    Int exponent;
    ExpressionVariablePtr expression;

    readConfig(config, "outputfileIonex",                fileNameOut,     Config::MUSTSET,  "", "");
    readConfig(config, "inputfileGriddedDataTimeSeries", fileNameIn,      Config::MUSTSET,  "", "must contain regular grid");
    readConfig(config, "value",                          expression,      Config::MUSTSET,  "data0", "expression (e.g. data column)");
    readConfig(config, "timeSeries",                     timeSeries,      Config::OPTIONAL, "", "(empty = use input file time series)");
    readConfig(config, "program",                        program,         Config::MUSTSET,  "GROOPS", "name of program (for first line)");
    readConfig(config, "institution",                    institution,     Config::MUSTSET,  "TUG (TU Graz)", "name of agency (for first line)");
    readConfig(config, "description",                    description,     Config::OPTIONAL, "", "description in header");
    readConfig(config, "comment",                        comment,         Config::OPTIONAL, R"(["TEC values in 0.1 TECU, 9999 if no value available"])", "comment in header");
    readConfig(config, "mappingFunction",                mappingFunction, Config::DEFAULT,  "NONE", "see IONEX documentation");
    readConfig(config, "elevationCutoff",                elevationCutoff, Config::DEFAULT,  "0", "see IONEX documentation (0 if unknown)");
    readConfig(config, "observablesUsed",                observablesUsed, Config::OPTIONAL, "", "see IONEX documentation");
    readConfig(config, "exponent",                       exponent,        Config::DEFAULT,  "-1", "factor 10^exponent is applied to all values");
    if(isCreateSchema(config)) return;

    logStatus<<"read gridded data time series file <"<<fileNameIn<<">"<<Log::endl;
    InFileGriddedDataTimeSeries infile(fileNameIn);
    GriddedDataRectangular grid;
    if(!grid.init(infile.grid()))
      throw(Exception("grid must be rectangular"));

    std::vector<Time> times = infile.times();
    if(timeSeries)
      times = timeSeries->times();

    logStatus << "write IONEX file <" << fileNameOut << ">" << Log::endl;
    OutFile outfile(fileNameOut);

    // Header: IONEX VERSION / TYPE
    outfile << "     1.0            IONOSPHERE MAPS     GNSS                IONEX VERSION / TYPE" << std::endl;

    // Header: PGM / RUN BY / DATE
    program.resize(20, ' ');
    institution.resize(20, ' ');
    outfile << program << institution << System::now()%"%D %H:%M:%S "s << "PGM / RUN BY / DATE" << std::endl;

    // Header: DESCRIPTION
    for(auto &d : description)
    {
      d.resize(60, ' ');
      outfile << d << "DESCRIPTION" << std::endl;
    }

    // Header: EPOCH OF FIRST/LAST MAP, INTERVAL, # OF MAPS IN FILE
    outfile << times.front()%"  %y    %m    %d    %H    %M    %S"s << std::string(24, ' ') << "EPOCH OF FIRST MAP" << std::endl;
    outfile << times.back()%"  %y    %m    %d    %H    %M    %S"s << std::string(24, ' ') << "EPOCH OF LAST MAP" << std::endl;
    outfile << medianSampling(times).seconds()%"% 6i"s << std::string(54, ' ') << "INTERVAL" << std::endl;
    outfile << times.size()%"% 6i"s << std::string(54, ' ') << "# OF MAPS IN FILE" << std::endl;

    // Header: MAPPING FUNCTION, ELEVATION CUTOFF, OBSERVABLES USED, BASE RADIUS, MAP DIMENSION
    mappingFunction.resize(4, ' ');
    observablesUsed.resize(60, ' ');
    outfile << "  " << mappingFunction << std::string(54, ' ') << "MAPPING FUNCTION" << std::endl;
    outfile << elevationCutoff%"%8.1f"s << std::string(52, ' ') << "ELEVATION CUTOFF" << std::endl;
    outfile << observablesUsed << "OBSERVABLES USED" << std::endl;
    outfile << (infile.grid().ellipsoid.a()/1000)%"%8.1f"s << std::string(52, ' ') << "BASE RADIUS" << std::endl;
    outfile << "     2" << std::string(54, ' ') << "MAP DIMENSION" << std::endl;

    // Header: HGT1 / HGT 2 / DHGT, LAT1 / LAT2 / DLAT, LON1 / LON2 / DLON
    const Double height = mean(Vector(grid.heights));
    UInt columns = (grid.longitudes.front()+2*PI == grid.longitudes.back() ? grid.longitudes.size()-1 : grid.longitudes.size());
    outfile << "  " << (height/1000)%"%6.1f%6.1f"s << 0%"%6.1f"s << std::string(40, ' ') << "HGT1 / HGT2 / DHGT" << std::endl;
    outfile << "  " << (grid.latitudes.front()*RAD2DEG)%"%6.1f"s << (grid.latitudes.back()*RAD2DEG)%"%6.1f"s
            << ((grid.latitudes.back()-grid.latitudes.front())/grid.latitudes.size()*RAD2DEG)%"%6.1f"s << std::string(40, ' ') << "LAT1 / LAT2 / DLAT" << std::endl;
    outfile << "  " << (grid.longitudes.front()*RAD2DEG)%"%6.1f"s << (grid.longitudes.back()*RAD2DEG)%"%6.1f"s
            << ((grid.longitudes.back()-grid.longitudes.front())/columns*RAD2DEG)%"%6.1f"s << std::string(40, ' ') << "LON1 / LON2 / DLON" << std::endl;

    // Header: EXPONENT
    outfile << exponent%"% 6i"s << std::string(54, ' ') << "EXPONENT" << std::endl;

    // Header: COMMENT
    for(auto &c : comment)
    {
      c.resize(60, ' ');
      outfile << c << "COMMENT" << std::endl;
    }

    // Header: END OF HEADER
    outfile << std::string(60, ' ') << "END OF HEADER" << std::endl;

    // Data records
    const Double factor = std::pow(10, exponent);
    for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
    {
      outfile << (idEpoch+1)%"% 6i"s << std::string(54, ' ') << "START OF TEC MAP" << std::endl;
      outfile << times.at(idEpoch)%"  %y    %m    %d    %H    %M    %S"s << std::string(24, ' ') << "EPOCH OF CURRENT MAP" << std::endl;
      for(UInt i = 0; i < grid.latitudes.size(); i++)
      {
        outfile << "  " << (grid.latitudes.at(i)*RAD2DEG)%"% 6.1f"s <<  (grid.longitudes.front()*RAD2DEG)%"%6.1f"s << (grid.longitudes.back()*RAD2DEG)%"%6.1f"s
                << ((grid.longitudes.back()-grid.longitudes.front())/columns*RAD2DEG)%"%6.1f"s << (height/1000)%"%6.1f"s << std::string(28, ' ') << "LAT/LON1/LON2/DLON/H" << std::endl;

        Matrix data = infile.data(times.at(idEpoch));
        VariableList varList;
        addDataVariables(data, varList);

        for(UInt j = 0; j < grid.longitudes.size(); j++)
        {
          evaluateDataVariables(data, i*grid.longitudes.size()+j, varList);
          Double value = expression->evaluate(varList);

          if(j > 0 && j%16 == 0) // start new line after 16 values
            outfile << std::endl;
          outfile << (std::isnan(value) ? " 9999" : (value/factor)%"% 5i"s);
        }
        outfile << std::endl;
      }
      outfile << (idEpoch+1)%"% 6i"s << std::string(54, ' ') << "END OF TEC MAP" << std::endl;
    }

    outfile << std::string(60, ' ') << "END OF FILE" << std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
