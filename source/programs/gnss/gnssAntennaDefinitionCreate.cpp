/***********************************************/
/**
* @file gnssAntennaDefinitionCreate.cpp
*
* @brief Create GNSS antenna definition file.
*
* @author Torsten Mayer-Guerr
* @date 2019-08-29
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create a \file{GNSS antenna definition file}{gnssAntennaDefinition} (Antenna Center Variations, ACV) consisting of multiple antennas.
The antennas can be created from scratch or can be selected from existing files.
This program can also be used to modify existing files.

Furthermore it can be used to create accuracy definition files containing azimuth and elevation dependent accuracy values for antennas.
To create an accuracy pattern for phase observations with \verb|1 mm| accuracy at zenith and no azimuth dependency, define a
pattern with \config{type}=\verb|L|, \config{values}=\verb|0.001/cos(zenith/rho)|.

The antennas in \configFile{outputfileAntennaDefinition}{gnssAntennaDefinition}
are sorted by names and duplicates are removed (first one is kept).
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "config/configRegister.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileGnssStationInfo.h"

/***** CLASS ***********************************/

/** @brief Create GNSS antenna definition file.
* @ingroup programsGroup */
class GnssAntennaDefinitionCreate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssAntennaDefinitionCreate, SINGLEPROCESS, "Create GNSS antenna definition file.", Gnss)

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, GnssAntennaPattern &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    Angle maxZenith, dAzimuth;
    ExpressionVariablePtr valueExpr;

    readConfig(config, "type",         var.type,       Config::MUSTSET,  "",   "pattern matching of observation types");
    readConfig(config, "offsetX",      var.offset.x(), Config::DEFAULT,  "0",  "[m] antenna center offset");
    readConfig(config, "offsetY",      var.offset.y(), Config::DEFAULT,  "0",  "[m] antenna center offset");
    readConfig(config, "offsetZ",      var.offset.z(), Config::DEFAULT,  "0",  "[m] antenna center offset");
    readConfig(config, "deltaAzimuth", dAzimuth,       Config::MUSTSET,  "5",  "[degree] step size");
    readConfig(config, "deltaZenith",  var.dZenit,     Config::MUSTSET,  "5",  "[degree] step size");
    readConfig(config, "maxZenith",    maxZenith,      Config::DEFAULT,  "90", "[degree]");
    readConfig(config, "values",       valueExpr,      Config::OPTIONAL, "0",  "[m] expression (zenith, azimuth: variables)");
    endSequence(config);
    if(isCreateSchema(config))
      return TRUE;

    UInt rows   = static_cast<UInt>(std::round(2*PI/dAzimuth));
    UInt cols   = static_cast<UInt>(std::round(maxZenith/var.dZenit))+1;
    var.pattern = Matrix(rows,cols);

    if(valueExpr)
    {
      auto varList = config.getVarList();
      addVariable("zenith",  varList);
      addVariable("azmiuth", varList);
      valueExpr->simplify(varList);

      for(UInt i=0; i<var.pattern.rows(); i++)
        for(UInt k=0; k<var.pattern.columns(); k++)
        {
          varList["azmiuth"]->setValue(i * 360./var.pattern.rows());
          varList["zenith"] ->setValue(k * RAD2DEG * Double(var.dZenit));
          var.pattern(i, k) = valueExpr->evaluate(varList);
        }
    }

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

// Latex documentation
static const char *docstringGnssAntennaDefintionList = R"(
\section{GnssAntennaDefintionList}\label{gnssAntennaDefintionListType}
Provides a list of GnssAntennaDefinitions as used in \program{GnssAntennaDefinitionCreate}.
)";

class GnssAntennaDefintionList;
typedef std::shared_ptr<GnssAntennaDefintionList> GnssAntennaDefintionListPtr;

class GnssAntennaDefintionList
{
protected:
  static Bool match(GnssType type1, GnssType type2)
  {
    return (type1 == type2) &&
           (!type1.hasWildcard(GnssType::PRN)       || type2.hasWildcard(GnssType::PRN))       &&
           (!type1.hasWildcard(GnssType::SYSTEM)    || type2.hasWildcard(GnssType::SYSTEM))    &&
           (!type1.hasWildcard(GnssType::FREQUENCY) || type2.hasWildcard(GnssType::FREQUENCY)) &&
           (!type1.hasWildcard(GnssType::TYPE)      || type2.hasWildcard(GnssType::TYPE))      &&
           (!type1.hasWildcard(GnssType::ATTRIBUTE) || type2.hasWildcard(GnssType::ATTRIBUTE)) &&
           (!type1.hasWildcard(GnssType::FREQ_NO)   || type2.hasWildcard(GnssType::FREQ_NO));
  }

public:
  std::vector<GnssAntennaDefinitionPtr> antennas;

  /// Creates an derived instance of this class.
  static GnssAntennaDefintionListPtr create(Config &config, const std::string &name);
};

/***********************************************/

static const char *docstringGnssAntennaDefintionListNew = R"(
\subsection{New}
Creates a new antenna.
)";

class GnssAntennaDefintionListNew : public GnssAntennaDefintionList
{
public:
  GnssAntennaDefintionListNew(Config &config)
  {
    try
    {
      antennas.push_back(GnssAntennaDefinitionPtr(new GnssAntennaDefinition()));
      readConfig(config, "name",    antennas.back()->name,    Config::OPTIONAL, "", "");
      readConfig(config, "serial",  antennas.back()->serial,  Config::OPTIONAL, "", "");
      readConfig(config, "radome",  antennas.back()->radome,  Config::OPTIONAL, "", "");
      readConfig(config, "comment", antennas.back()->comment, Config::OPTIONAL, "", "");
      readConfig(config, "pattern", antennas.back()->pattern, Config::MUSTSET,  "", "");
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssAntennaDefintionListFromFile = R"(
\subsection{FromFile}
Select all or the first antenna from an \file{antenna definition file}{gnssAntennaDefinition}
which matches the wildcards.
)";

class GnssAntennaDefintionListFromFile : public GnssAntennaDefintionList
{
public:
  GnssAntennaDefintionListFromFile(Config &config)
  {
    try
    {
      FileName    fileNameAntenna;
      std::string name, serial, radome;
      Bool        onlyFirstMatch;

      readConfig(config, "inputfileAntennaDefinition", fileNameAntenna, Config::MUSTSET,  "",  "");
      readConfig(config, "name",                       name,            Config::OPTIONAL, "*", "");
      readConfig(config, "serial",                     serial,          Config::OPTIONAL, "*", "");
      readConfig(config, "radome",                     radome,          Config::OPTIONAL, "*", "");
      readConfig(config, "onlyFirstMatch",             onlyFirstMatch,  Config::DEFAULT,  "0", "otherwise all machting antennas included");
      if(isCreateSchema(config)) return;

      std::vector<GnssAntennaDefinitionPtr> antennasFile;
      readFileGnssAntennaDefinition(fileNameAntenna, antennasFile);

      std::regex pattern = String::wildcard2regex(GnssAntennaDefinition::str(name, serial, radome));
      for(auto &antenna : antennasFile)
        if(std::regex_match(antenna->str(), pattern))
        {
          antennas.push_back(antenna);
          if(onlyFirstMatch)
            break;
        }
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssAntennaDefintionListFromStationInfo = R"(
\subsection{FromStationInfo}
Select all antennas from an \file{antenna definition file}{gnssAntennaDefinition}
which are used by a station within a defined time interval.
With \config{specializeAntenna} an individual antenna is created for each different serial number
using the general type specific values from file.
)";

class GnssAntennaDefintionListFromStationInfo : public GnssAntennaDefintionList
{
public:
  GnssAntennaDefintionListFromStationInfo(Config &config)
  {
    try
    {
      FileName fileNameStationInfo, fileNameAntenna;
      Time     timeStart, timeEnd = date2time(2500,1,1);
      Bool     specialize;

      readConfig(config, "inputfileStationInfo",       fileNameStationInfo, Config::MUSTSET,  "",  "");
      readConfig(config, "inputfileAntennaDefinition", fileNameAntenna,     Config::MUSTSET,  "",  "");
      readConfig(config, "timeStart",                  timeStart,           Config::OPTIONAL, "",  "only antennas used in this time interval");
      readConfig(config, "timeEnd",                    timeEnd,             Config::OPTIONAL, "",  "only antennas used in this time interval");
      readConfig(config, "specializeAntenna",          specialize,          Config::DEFAULT,  "0", "e.g. separate different serial numbers from stationInfo");
      if(isCreateSchema(config)) return;

      GnssStationInfo stationInfo;
      std::vector<GnssAntennaDefinitionPtr> antennasFile;
      readFileGnssStationInfo(fileNameStationInfo, stationInfo);
      readFileGnssAntennaDefinition(fileNameAntenna, antennasFile);
      stationInfo.fillAntennaPattern(antennasFile);

      for(const auto &antennaInfo : stationInfo.antenna)
        if(antennaInfo.timeEnd >= timeStart && antennaInfo.timeStart < timeEnd)
        {
          if(antennaInfo.antennaDef)
          {
            antennas.push_back(std::make_shared<GnssAntennaDefinition>(*antennaInfo.antennaDef));
            if(specialize)
            {
              if(antennas.back()->name.empty())   antennas.back()->name   = antennaInfo.name;
              if(antennas.back()->serial.empty()) antennas.back()->serial = antennaInfo.serial;
              if(antennas.back()->radome.empty()) antennas.back()->radome = antennaInfo.radome;
            }
          }
          else
            logWarning<<stationInfo.markerName<<"."<<stationInfo.markerNumber<<": no antenna definition found for "<<antennaInfo.str()<<Log::endl;
        }
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssAntennaDefintionListResample = R"(
\subsection{Resample}
The azimuth and elevation dependend antenna center variations (patterns) of all \config{antenna}s
are resampled to a new resolution.
)";

class GnssAntennaDefintionListResample : public GnssAntennaDefintionList
{
public:
  GnssAntennaDefintionListResample(Config &config)
  {
    try
    {
      std::vector<GnssAntennaDefintionListPtr> antennaLists;
      Angle dAzimuth (NAN_EXPR);
      Angle dZenith  (NAN_EXPR);
      Angle maxZenith(NAN_EXPR);

      readConfig(config, "antenna",      antennaLists, Config::MUSTSET,   "", "");
      readConfig(config, "deltaAzimuth", dAzimuth,     Config::OPTIONAL, "", "[degree] step size, empty: no change");
      readConfig(config, "deltaZenith",  dZenith,      Config::OPTIONAL, "", "[degree] step size, empty: no change");
      readConfig(config, "maxZenith",    maxZenith,    Config::OPTIONAL, "", "[degree], empty: no change");
      if(isCreateSchema(config)) return;

      for(auto &antennaList : antennaLists)
        for(auto &antenna : antennaList->antennas)
          for(auto &pattern : antenna->pattern)
          {
            Double dz   = (dZenith   > 0) ? dZenith   : pattern.dZenit;
            Double maxz = (maxZenith >=0) ? maxZenith : ((pattern.pattern.columns()-1) * Double(pattern.dZenit));
            UInt cols   = static_cast<UInt>(std::round(maxz/dz))+1;
            UInt rows   = (dAzimuth  > 0) ? static_cast<UInt>(std::round(2*PI/dAzimuth)) : pattern.pattern.rows();

            if((pattern.dZenit != dz) || (pattern.pattern.rows() != rows) || (pattern.pattern.columns() != cols))
            {
              Matrix acv(rows, cols);
              for(UInt i=0; i<rows; i++)
                for(UInt k=0; k<cols; k++)
                  acv(i, k) += pattern.antennaVariations(Angle(2*PI*i/rows), Angle(PI/2-k*dz), FALSE);
              pattern.dZenit  = dz;
              pattern.pattern = acv;
            }
          }

      for(auto &antennaList : antennaLists)
        antennas.insert(antennas.end(), antennaList->antennas.begin(), antennaList->antennas.end());
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssAntennaDefintionListTransform = R"(
\subsection{Transform}
This class can be used to separate general antenna patterns for different \configClass{gnssType}{gnssType}s.
If the \config{antenna}s contain only one pattern for all GPS observations on the L1 frequency (\verb|*1*G**|),
the \config{patternTypes}=\verb|C1*G**| and \verb|L1*G**| create two patterns with the \verb|*1*G**| patterm as template.
The first matching pattern in the \config{antenna} is used as template.
Also new \config{additionalPattern} can be added (e.g. for \verb|*5*G**|).
With \config{addExistingPatterns} all already existing patterns that don't match completely to any of the above are added.
)";

class GnssAntennaDefintionListTransform : public GnssAntennaDefintionList
{
public:
  GnssAntennaDefintionListTransform(Config &config)
  {
    try
    {
      std::vector<GnssAntennaDefintionListPtr> antennaLists;
      std::vector<GnssType>                    types;
      std::vector<GnssAntennaPattern>          patternsAdd;
      Bool                                     addExistingPatterns;

      readConfig(config, "antenna",             antennaLists,        Config::MUSTSET,  "", "");
      readConfig(config, "patternTypes",        types,               Config::OPTIONAL, "",  "gnssType for each pattern (first match is used)");
      readConfig(config, "additionalPattern",   patternsAdd,         Config::OPTIONAL, "",  "additional new patterns");
      readConfig(config, "addExistingPatterns", addExistingPatterns, Config::DEFAULT,  "1", "add existing patterns that don't match completely any of the above");
      if(isCreateSchema(config)) return;

      for(auto &antennaList : antennaLists)
        for(auto &antenna : antennaList->antennas)
        {
          std::vector<GnssAntennaPattern> patternsOld = antenna->pattern;
          antenna->pattern.clear();
          for(GnssType type : types)
          {
            Bool found = FALSE;
            for(auto &patternOld : patternsOld)
              if(match(type, patternOld.type))
              {
                antenna->pattern.push_back(patternOld);
                antenna->pattern.back().type = type;
                found = TRUE;
                break;
              }
            if(!found)
              throw(Exception(antenna->str()+": no pattern found which fits to all possible types: "+type.str()));
          }

          for(auto &patternAdd : patternsAdd)
          {
            Bool found = FALSE;
            for(auto &patternNew : antenna->pattern)
              if(match(patternAdd.type, patternNew.type))
              {
                found = TRUE;
                break;
              }
            if(!found)
              antenna->pattern.push_back(patternAdd);
          }

          if(addExistingPatterns)
            for(auto &patternOld : patternsOld)
            {
              Bool found = FALSE;
              for(auto &patternNew : antenna->pattern)
                if(match(patternOld.type, patternNew.type))
                {
                  found = TRUE;
                  break;
                }
              if(!found)
                antenna->pattern.push_back(patternOld);
            }

          std::sort(antenna->pattern.begin(), antenna->pattern.end(), [](GnssAntennaPattern &a, GnssAntennaPattern &b){return a.type < b.type;});
        } // for(antennas)

      for(auto &antennaList : antennaLists)
        antennas.insert(antennas.end(), antennaList->antennas.begin(), antennaList->antennas.end());
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssAntennaDefintionListRename = R"(
\subsection{Rename}
Replaces parts of the descrption of \config{antenna}s.
The star "\verb|*|" left this part untouched.
)";

class GnssAntennaDefintionListRename : public GnssAntennaDefintionList
{
public:
  GnssAntennaDefintionListRename(Config &config)
  {
    try
    {
      std::vector<GnssAntennaDefintionListPtr> antennaLists;
      std::string name, serial, radome, comment;

      readConfig(config, "antenna", antennaLists, Config::MUSTSET,  "",  "");
      readConfig(config, "name",    name,         Config::OPTIONAL, "*", "*: left this part untouched");
      readConfig(config, "serial",  serial,       Config::OPTIONAL, "*", "*: left this part untouched");
      readConfig(config, "radome",  radome,       Config::OPTIONAL, "*", "*: left this part untouched");
      readConfig(config, "comment", comment,      Config::OPTIONAL, "*", "*: left this part untouched");
      if(isCreateSchema(config)) return;

      for(auto &antennaList : antennaLists)
        for(auto &antenna : antennaList->antennas)
        {
          if(name    != "*") antenna->name    = name;
          if(serial  != "*") antenna->serial  = serial;
          if(radome  != "*") antenna->radome  = radome;
          if(comment != "*") antenna->comment = comment;
          antennas.push_back(antenna);
        } // for(antennas)
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

GROOPS_REGISTER_CLASS(GnssAntennaDefintionList, "gnssAntennaDefintionListType",
                      GnssAntennaDefintionListNew,
                      GnssAntennaDefintionListFromFile,
                      GnssAntennaDefintionListFromStationInfo,
                      GnssAntennaDefintionListResample,
                      GnssAntennaDefintionListTransform,
                      GnssAntennaDefintionListRename)

GROOPS_READCONFIG_CLASS(GnssAntennaDefintionList, "gnssAntennaDefintionListType")

/***********************************************/

GnssAntennaDefintionListPtr GnssAntennaDefintionList::create(Config &config, const std::string &name)
{
  try
  {
    GnssAntennaDefintionListPtr ptr;
    std::string           choice;
    readConfigChoice(config, name, choice, Config::MUSTSET, "", "");
    if(readConfigChoiceElement(config, "new",             choice, ""))
      ptr = GnssAntennaDefintionListPtr(new GnssAntennaDefintionListNew(config));
    if(readConfigChoiceElement(config, "fromFile",        choice, ""))
      ptr = GnssAntennaDefintionListPtr(new GnssAntennaDefintionListFromFile(config));
    if(readConfigChoiceElement(config, "fromStationInfo", choice, ""))
      ptr = GnssAntennaDefintionListPtr(new GnssAntennaDefintionListFromStationInfo(config));
    if(readConfigChoiceElement(config, "resample",        choice, ""))
      ptr = GnssAntennaDefintionListPtr(new GnssAntennaDefintionListResample(config));
    if(readConfigChoiceElement(config, "transform",       choice, ""))
      ptr = GnssAntennaDefintionListPtr(new GnssAntennaDefintionListTransform(config));
    if(readConfigChoiceElement(config, "rename",          choice, ""))
      ptr = GnssAntennaDefintionListPtr(new GnssAntennaDefintionListRename(config));
    endChoice(config);
    return ptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssAntennaDefinitionCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutAntenna;
    std::vector<GnssAntennaDefintionListPtr> antennaLists;

    readConfig(config, "outputfileAntennaDefinition", fileNameOutAntenna, Config::MUSTSET, "", "");
    readConfig(config, "antenna",                     antennaLists,       Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    // ==================================

    std::vector<GnssAntennaDefinitionPtr> antennas;
    for(const auto &antennaList : antennaLists)
      antennas.insert(antennas.end(), antennaList->antennas.begin(), antennaList->antennas.end());

    std::stable_sort(antennas.begin(), antennas.end(), [](GnssAntennaDefinitionPtr a, GnssAntennaDefinitionPtr b)
                     {return (a->name != b->name) ? (a->name < b->name) : ((a->serial != b->serial) ? (a->serial < b->serial) : (a->radome < b->radome));});
    antennas.erase(std::unique(antennas.begin(), antennas.end(), [](GnssAntennaDefinitionPtr a, GnssAntennaDefinitionPtr b) {return (a->name == b->name) && (a->serial == b->serial) && (a->radome == b->radome);}), antennas.end());


    logStatus<<"write antenna definition <"<<fileNameOutAntenna<<">"<<Log::endl;
    writeFileGnssAntennaDefinition(fileNameOutAntenna, antennas);

    for(const auto &antenna : antennas)
    {
      std::string str;
      for(const auto &pattern : antenna->pattern)
        str += " "+pattern.type.str();
      logInfo<<" "<<antenna->str()<<str<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
