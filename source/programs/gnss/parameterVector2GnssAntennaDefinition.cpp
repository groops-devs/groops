/***********************************************/
/**
* @file parameterVector2GnssAntennaDefinition.cpp
*
* @brief Update antenna definition from parametrization.
*
* @author Torsten Mayer-Guerr
* @date 2012-12-05
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Updates an \file{GnssAntennaDefinition file}{gnssAntennaDefinition} with estimated parameters which belongs
to the parametrization \configClass{antennaCenterVariations}{parametrizationGnssAntennaType}.
The \configFile{outfileAntennaDefinition}{gnssAntennaDefinition} contains all antennas
from \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition}.
The antenna center variations representend by the \configFile{inputfileSolution}{matrix} are added
to the matching antennas.

The \file{GnssAntennaDefinition file}{gnssAntennaDefinition} can be modified to the demands before with
\program{GnssAntennaDefinitionCreate}.

The following steps are used to estimate antenna center variations:
\begin{itemize}
\item \program{GnssAntennaDefinitionCreate} or \program{GnssAntex2AntennaDefinition}
\item \program{GnssProcessing} with \config{inputfileAntennaDefinition} as apriori
      and writing \file{normal equations}{normalEquation} with
      parametrization of \configClass{antennaCenterVariations}{parametrizationGnssAntennaType}
\item \program{NormalsEliminate}: eliminate all other than antenna parameters
\item \program{NormalsAccumulate}: accumulate normals over a sufficient long period
\item \program{GnssAntennaNormalsConstraint}: constrain unsolvable parameter linear combinations
\item \program{NormalsSolverVCE}: estimate the parameter vector
\item \program{ParameterVector2GnssAntennaDefinition}: update \config{inputfileAntennaDefinition}
\end{itemize}

See also \program{ParameterVector2GnssAntennaDefinition}, \program{GnssAntennaNormalsConstraint}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileMatrix.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileParameterName.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"

/***** CLASS ***********************************/

/** @brief Update antenna definition from parametrization.
* @ingroup programsGroup */
class ParameterVector2GnssAntennaDefinition
{
public:

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(ParameterVector2GnssAntennaDefinition, SINGLEPROCESS, "Update antenna definition from parametrization", Gnss)
GROOPS_RENAMED_PROGRAM(GnssRepresentation2AntennaDefinition, ParameterVector2GnssAntennaDefinition, date2time(2020, 6, 26))

/***********************************************/

void ParameterVector2GnssAntennaDefinition::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameAntennaOut;
    FileName fileNameAntennaIn, fileNameSolution, fileNameParameterNames;
    ParametrizationGnssAntennaPtr parametrization;

    readConfig(config, "outfileAntennaDefinition",   fileNameAntennaOut,     Config::MUSTSET, "", "all apriori antennas");
    readConfig(config, "inputfileAntennaDefinition", fileNameAntennaIn,      Config::MUSTSET, "", "apriori antennas");
    readConfig(config, "antennaCenterVariations",    parametrization,        Config::MUSTSET, "", "");
    readConfig(config, "inputfileSolution",          fileNameSolution,       Config::MUSTSET, "", "");
    readConfig(config, "inputfileParameterNames",    fileNameParameterNames, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    // ============================

    logStatus<<"read antenna definitions <"<<fileNameAntennaIn<<">"<<Log::endl;
    std::vector<GnssAntennaDefinitionPtr> antennas;
    readFileGnssAntennaDefinition(fileNameAntennaIn, antennas);

    logStatus<<"read solution vector <"<<fileNameSolution<<">"<<Log::endl;
    Vector solution;
    readFileMatrix(fileNameSolution, solution);

    logStatus<<"read parameter names <"<<fileNameParameterNames<<">"<<Log::endl;
    std::vector<ParameterName> parameterNames;
    readFileParameterName(fileNameParameterNames, parameterNames);

    std::vector<ParameterName> parametrizationNames;
    parametrization->parameterName(parametrizationNames);

    // ============================

    for(UInt i=0; i<parameterNames.size(); i++)
      if(String::startsWith(parameterNames.at(i).type, parametrizationNames.at(0).type))
      {
        // check antenna name
        // ------------------
        std::vector<std::string> parts = String::split(parameterNames.at(i).object, GnssAntennaDefinition::sep);
        auto antenna = GnssAntennaDefinition::find(antennas, parts.at(0), parts.at(1), parts.at(2));
        if(!antenna)
          throw(Exception("antenna not found in list: "+parameterNames.at(i).object));

        // extract GnssType
        // ----------------
        auto pos = parameterNames.at(i).type.rfind('.');
        if((pos == std::string::npos) || (pos+1 == parameterNames.at(i).type.size()))
          throw(Exception(parameterNames.at(i).str() + ": GnssType not found in parameter name"));
        GnssType type(parameterNames.at(i).type.substr(pos+1));

        // add parametrization to all machting patterns
        // -------------------------------------------
        Bool found = FALSE;
        for(auto &pattern : antenna->patterns)
          if(type == pattern.type)
          {
            logInfo<<antenna->str()<<": add parametrization of "<<type.str()<<" to "<<pattern.type.str()<<Log::endl;
            found = TRUE;

            Vector x = solution.row(i, parametrization->parameterCount());
            for(UInt k=0; k<x.rows(); k++)
            {
              if(parameterNames.at(i+k).type.find("antennaCenter.x") == 0) {pattern.offset.x() += x(k); x(k) = 0;}
              if(parameterNames.at(i+k).type.find("antennaCenter.y") == 0) {pattern.offset.y() += x(k); x(k) = 0;}
              if(parameterNames.at(i+k).type.find("antennaCenter.z") == 0) {pattern.offset.z() += x(k); x(k) = 0;}
            }
            for(UInt i=0; i<pattern.pattern.rows(); i++)
              for(UInt k=0; k<pattern.pattern.columns(); k++)
                pattern.pattern(i,k) += inner(parametrization->designMatrix(Angle(2*PI*i/pattern.pattern.rows()), Angle(PI/2-k*Double(pattern.dZenit))).trans(), x);

          // type cannot match further patterns (all wildcards in type are already catched)?
          if((!type.hasWildcard(GnssType::PRN)       || pattern.type.hasWildcard(GnssType::PRN)      ) &&
             (!type.hasWildcard(GnssType::SYSTEM)    || pattern.type.hasWildcard(GnssType::SYSTEM)   ) &&
             (!type.hasWildcard(GnssType::FREQUENCY) || pattern.type.hasWildcard(GnssType::FREQUENCY)) &&
             (!type.hasWildcard(GnssType::TYPE)      || pattern.type.hasWildcard(GnssType::TYPE)     ) &&
             (!type.hasWildcard(GnssType::ATTRIBUTE) || pattern.type.hasWildcard(GnssType::ATTRIBUTE)) &&
             (!type.hasWildcard(GnssType::FREQ_NO)   || pattern.type.hasWildcard(GnssType::FREQ_NO)  ))
            break;
          }

        if(!found)
          logWarning<<antenna->str()<<" has no pattern for "<<type.str()<<Log::endl;
      }

    // ============================

    if(!fileNameAntennaOut.empty())
    {
      logStatus<<"write antenna definition <"<<fileNameAntennaOut<<">"<<Log::endl;
      writeFileGnssAntennaDefinition(fileNameAntennaOut, antennas);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
