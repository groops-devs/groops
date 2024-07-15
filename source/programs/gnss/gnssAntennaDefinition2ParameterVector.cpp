/***********************************************/
/**
* @file gnssAntennaDefinition2ParameterVector.cpp
*
* @brief Estimate parameters of an antenna parametrization to fit an antenna definition.
*
* @author Torsten Mayer-Guerr
* @date 2019-09-08
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Estimates parameters of a parametrization of \configClass{antennaCenterVariations}{parametrizationGnssAntennaType},
which represents all antennas from \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition}
matching the wildcard patterns of \config{name}, \config{serial}, \config{radome}.

The provided values at the area weighted grid points of the pattern of each gnssType are used as pseudo-observations.
A subset of patterns can be selected with \configClass{types}{gnssType}.

The \file{GnssAntennaDefinition file}{gnssAntennaDefinition} can be modified to the demands before with
\program{GnssAntennaDefinitionCreate}.

See also \program{ParameterVector2GnssAntennaDefinition}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileMatrix.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileParameterName.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"

/***** CLASS ***********************************/

/** @brief Estimate parameters of an antenna parametrization to fit an antenna definition.
* @ingroup programsGroup */
class GnssAntennaDefinition2ParameterVector
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssAntennaDefinition2ParameterVector, SINGLEPROCESS, "Estimate parameters of an antenna parametrization to fit an antenna definition", Gnss)
GROOPS_RENAMED_PROGRAM(GnssAntennaDefinition2Representation, GnssAntennaDefinition2ParameterVector, date2time(2020, 6, 26))

/***********************************************/

void GnssAntennaDefinition2ParameterVector::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                      fileNameSolution, fileNameParameterNames;
    ParametrizationGnssAntennaPtr parametrization;
    FileName                      fileNameAntenna;
    std::string                   name, serial, radome;
    std::vector<GnssType>         types;

    readConfig(config, "outputfileSolution",         fileNameSolution,       Config::MUSTSET,  "",  "");
    readConfig(config, "outputfileParameterNames",   fileNameParameterNames, Config::OPTIONAL, "",  "");
    readConfig(config, "antennaCenterVariations",    parametrization,        Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileAntennaDefinition", fileNameAntenna,        Config::MUSTSET,  "",  "");
    readConfig(config, "name",                       name,                   Config::OPTIONAL, "*", "");
    readConfig(config, "serial",                     serial,                 Config::OPTIONAL, "*", "");
    readConfig(config, "radome",                     radome,                 Config::OPTIONAL, "*", "");
    readConfig(config, "types",                      types,                  Config::OPTIONAL, "",  "if not set, all types in the file are used");
    if(isCreateSchema(config)) return;

    // ============================

    logStatus<<"read antenna center variations <"<<fileNameAntenna<<">"<<Log::endl;
    std::vector<GnssAntennaDefinitionPtr> antennaList;
    readFileGnssAntennaDefinition(fileNameAntenna, antennaList);

    std::regex                 pattern = String::wildcard2regex(GnssAntennaDefinition::str(name, serial, radome));
    std::vector<Vector>        solutions;
    std::vector<ParameterName> parameterNames;
    std::vector<ParameterName> baseNames;
    parametrization->parameterName(baseNames);

    for(GnssAntennaDefinitionPtr antenna : antennaList)
      if(std::regex_match(antenna->str(), pattern))
      {
        const std::string antennaName = antenna->str();
        logStatus<<"estimate parameters for "<<antennaName<<Log::endl;

        auto typesAnt = types;
        if(typesAnt.size() == 0)
          for(auto &pattern : antenna->patterns)
            typesAnt.push_back(pattern.type);

        Vector x(typesAnt.size()*parametrization->parameterCount());
        for(UInt idType=0; idType<typesAnt.size(); idType++)
        {
          logInfo<<typesAnt.at(idType).str()<<Log::endl;
          const UInt idPattern = antenna->findAntennaPattern(typesAnt.at(idType), GnssAntennaDefinition::THROW_EXCEPTION);
          GnssAntennaPattern &pattern = antenna->patterns.at(idPattern);

          Matrix A(pattern.pattern.rows()*pattern.pattern.columns(), parametrization->parameterCount());
          Vector l(pattern.pattern.rows()*pattern.pattern.columns());
          UInt idx = 0;
          for(UInt i=0; i<pattern.pattern.rows(); i++)
            for(UInt k=0; k<pattern.pattern.columns(); k++)
            {
              // area weight
              const Double p = std::sqrt((k == 0) ? (1-std::cos(0.5*Double(pattern.dZenit))) : std::sin(k*Double(pattern.dZenit)));
              l(idx) = p * pattern.pattern(i,k);
              axpy(p, parametrization->designMatrix(Angle(i*2*PI/pattern.pattern.rows()), Angle(PI/2-k*Double(pattern.dZenit))), A.row(idx));
              idx++;
            }
          solutions.push_back(leastSquares(A, l));

          std::string typeStr = "." + typesAnt.at(idType).str();
          for(const auto &base : baseNames)
            parameterNames.push_back( ParameterName(antennaName, base.type + typeStr, base.temporal, base.interval) );
        }
      } // for(antenna)

    // ============================

    if(!fileNameSolution.empty())
    {
      logStatus<<"write solution to <"<<fileNameSolution<<">"<<Log::endl;
      const UInt count = std::accumulate(solutions.begin(), solutions.end(), UInt(0), [](UInt count, const Vector &x) {return count + x.size();});
      Vector x(count);
      UInt idx = 0;
      for(const Vector &solution : solutions)
      {
        copy(solution, x.row(idx, solution.rows()));
        idx += solution.rows();
      }
      writeFileMatrix(fileNameSolution, x);
    }

    if(!fileNameParameterNames.empty())
    {
      logStatus<<"write parameter names to <"<<fileNameParameterNames<<">"<<Log::endl;
      writeFileParameterName(fileNameParameterNames, parameterNames);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
