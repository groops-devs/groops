/***********************************************/
/**
* @file gnssAntennaNormalsConstraint.cpp
*
* @brief Apply constraints to normals of antenna parametrization.
*
* @author Torsten Mayer-Guerr
* @date 2019-09-14
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Apply constraints to \file{normal equations}{normalEquation}
containing \configClass{antennaCenterVariations}{parametrizationGnssAntennaType}.
Usually the antenna center variations are estimated together with other parameters
like station coordinates, signal biases and slant TEC in \program{GnssProcessing}.
This results in a rank deficient matrix as not all parameters can be separated.
The deficient can be solved by adding pseudo observation equations as constraints.

To separate antenna center variations and signal biases
apply \config{constraint:mean} for each GNSS \configClass{type}{gnssType}.
The observation equation for the integral mean of antenna center variations (ACV)
in all azimuth~$A$ and elevation~$E$ dependent directions
\begin{equation}
  0 = \iint ACV(A,E)\, d\Phi \approx \sum_i ACV(A_i,E_i)\, \Delta\Phi_i
\end{equation}
is approximated by a grid defined by
\config{deltaAzimuth}, \config{deltaZenith}, and \config{maxZenith}.

To separate from station coordinates use \config{constraint:centerMean}
and from slant TEC parameters use \config{constraint:TEC}.

The constraints are applied separately to all antennas matching
the wildcard patterns of \config{name}, \config{serial}, \config{radome}.

See also \program{ParameterVector2GnssAntennaDefinition}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileMatrix.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileParameterName.h"
#include "files/fileNormalEquation.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"

/***** CLASS ***********************************/

/** @brief Apply constraints to normals of antenna parametrization.
* @ingroup programsGroup */
class GnssAntennaNormalsConstraint
{
public:
  class Constraint
  {
  public:
    enum Type {CENTER, CENTERMEAN, CONSTANT, CONSTANTMEAN, TEC};
    Type     type;
    Bool     applyWeight;
    GnssType gnssType;
    Double   sigma;
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssAntennaNormalsConstraint, SINGLEPROCESS, "Apply constraints to normals of antenna parametrization", Gnss)
GROOPS_RENAMED_PROGRAM(GnssAntennaRepresentationConstraint, GnssAntennaNormalsConstraint, date2time(2020, 6, 26))

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssAntennaNormalsConstraint::Constraint &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string choice;
  if(!readConfigChoice(config, name, choice, mustSet, defaultValue, annotation))
    return FALSE;
  if(readConfigChoiceElement(config, "center", choice, "zero center (x,y,z) of a single pattern"))
  {
    var.type = GnssAntennaNormalsConstraint::Constraint::CENTER;
    readConfig(config, "type",         var.gnssType,    Config::MUSTSET,  "",     "applied for each matching types");
    readConfig(config, "applyWeight",  var.applyWeight, Config::DEFAULT,  "1",    "from normal equations");
    readConfig(config, "sigma",        var.sigma,       Config::DEFAULT,  "1e-5", "[m]");
  }
  if(readConfigChoiceElement(config, "centerMean", choice, "zero center (x,y,z) as (weighted) mean of all patterns"))
  {
    var.type = GnssAntennaNormalsConstraint::Constraint::CENTERMEAN;
    readConfig(config, "applyWeight",  var.applyWeight, Config::DEFAULT,  "1",    "from normal equations");
    readConfig(config, "sigma",        var.sigma,       Config::DEFAULT,  "1e-5", "[m]");
  }
  if(readConfigChoiceElement(config, "constant", choice, "zero constant (mean of all directions) of a single pattern"))
  {
    var.type = GnssAntennaNormalsConstraint::Constraint::CONSTANT;
    readConfig(config, "type",         var.gnssType,    Config::MUSTSET,  "",     "applied for each matching types");
    readConfig(config, "applyWeight",  var.applyWeight, Config::DEFAULT,  "1",    "from normal equations");
    readConfig(config, "sigma",        var.sigma,       Config::DEFAULT,  "1e-5", "[m]");
  }
  if(readConfigChoiceElement(config, "constantMean", choice, "zero constant (mean of all directions) as (weighted) mean of all patterns"))
  {
    var.type = GnssAntennaNormalsConstraint::Constraint::CONSTANTMEAN;
    readConfig(config, "applyWeight",  var.applyWeight, Config::DEFAULT,  "1",    "from normal equations");
    readConfig(config, "sigma",        var.sigma,       Config::DEFAULT,  "1e-5", "[m]");
  }
  if(readConfigChoiceElement(config, "TEC", choice, "zero TEC computed as (weighetd) least squares from all types"))
  {
    var.type = GnssAntennaNormalsConstraint::Constraint::TEC;
    readConfig(config, "applyWeight",  var.applyWeight, Config::DEFAULT,  "1",    "from normal equations");
    readConfig(config, "sigma",        var.sigma,       Config::DEFAULT,  "1e-5", "[TECU]");
  }
  endChoice(config);
  return TRUE;
}

/***********************************************/

void GnssAntennaNormalsConstraint::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                      fileNameNormalsOut;
    FileName                      fileNameNormalsIn;
    std::vector<Constraint>       constraints;
    ParametrizationGnssAntennaPtr parametrization;
    std::string                   name, serial, radome;
    Angle                         dAzimuth, dZenith, maxZenith;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "inputfileNormalequation",  "inputfileNormalEquation",  date2time(2020, 6, 3));

    readConfig(config, "outputfileNormalEquation", fileNameNormalsOut, Config::MUSTSET,  "",   "with applied constraints");
    readConfig(config, "inputfileNormalEquation",  fileNameNormalsIn,  Config::MUSTSET,  "",   "");
    readConfig(config, "constraint",               constraints,        Config::MUSTSET,  "",   "");
    readConfig(config, "antennaCenterVariations",  parametrization,    Config::MUSTSET,  "",   "");
    readConfig(config, "antennaName",              name,               Config::OPTIONAL, "*",  "apply constraints to all machting antennas");
    readConfig(config, "antennaSerial",            serial,             Config::OPTIONAL, "*",  "apply constraints to all machting antennas");
    readConfig(config, "antennaRadome",            radome,             Config::OPTIONAL, "*",  "apply constraints to all machting antennas");
    readConfig(config, "deltaAzimuth",             dAzimuth,           Config::DEFAULT,  "1",  "[degree] sampling of pattern to estimate center/constant");
    readConfig(config, "deltaZenith",              dZenith,            Config::DEFAULT,  "1",  "[degree] sampling of pattern to estimate center/constant");
    readConfig(config, "maxZenith",                maxZenith,          Config::DEFAULT,  "90", "[degree] sampling of pattern to estimate center/constant");
    if(isCreateSchema(config)) return;

    // ============================

    logStatus<<"read normal equations <"<<fileNameNormalsIn<<">"<<Log::endl;
    Matrix             normals, rhs;
    NormalEquationInfo info;
    readFileNormalEquation(fileNameNormalsIn, info, normals, rhs);

    // ============================

    // find all antennas in normals/parameterNames
    // -------------------------------------------
    std::vector<std::string>           antennaName;
    std::vector<std::vector<GnssType>> types; // for each antenna
    std::vector<std::vector<UInt>>     index; // for each antenna, type: index in normals

    std::vector<ParameterName> parametrizationNames;
    parametrization->parameterName(parametrizationNames);
    std::regex wildcard = String::wildcard2regex(GnssAntennaDefinition::str(name, serial, radome));

    for(UInt i=0; i<info.parameterName.size(); i++)
      if((info.parameterName.at(i).type.find(parametrizationNames.at(0).type) == 0) && std::regex_match(info.parameterName.at(i).object, wildcard))
      {
        // check antenna name
        // ------------------
        const UInt idAnt = std::distance(antennaName.begin(), std::find(antennaName.begin(), antennaName.end(), info.parameterName.at(i).object));
        if(idAnt >= antennaName.size())
        {
          antennaName.push_back(info.parameterName.at(i).object);
          types.push_back(std::vector<GnssType>());
          index.push_back(std::vector<UInt>());
        }

        // extract GnssType
        auto pos = info.parameterName.at(i).type.rfind('.');
        if((pos == std::string::npos) || (pos+1 == info.parameterName.at(i).type.size()))
          throw(Exception(info.parameterName.at(i).str() + ": GnssType not found in parameter name"));
        types.at(idAnt).push_back(GnssType(info.parameterName.at(i).type.substr(pos+1)));
        index.at(idAnt).push_back(i);
      }

    // ============================

    // partial derivatives of parametrization coeffcients with respect to center (x,y,z)
    // performed by numerical integration and least squares adjustment
    Matrix dacv_dcenter;
    {
      logStatus<<"compute center constraint"<<Log::endl;
      const UInt rows = static_cast<UInt>(std::round(2*PI/dAzimuth));
      const UInt cols = static_cast<UInt>(std::round(maxZenith/dZenith));

      Matrix A(rows*cols, parametrization->parameterCount());
      Matrix L(rows*cols, 3);
      UInt idx = 0;
      for(UInt i=0; i<rows; i++)
        for(UInt k=0; k<cols; k++)
        {
          const Angle azimuth  (i*2*PI/rows);
          const Angle elevation(PI/2-(k+0.5)*Double(dZenith));
          copy(parametrization->designMatrix(azimuth, elevation), A.row(idx));
          L(idx, 0) = -std::cos(elevation) * std::cos(azimuth);  // x
          L(idx, 1) = -std::cos(elevation) * std::sin(azimuth);  // y
          L(idx, 2) = -std::sin(elevation);                      // z
          A.row(idx) *= std::sqrt(std::cos(elevation));          // area weights
          L.row(idx) *= std::sqrt(std::cos(elevation));
          idx++;
        }
      dacv_dcenter = leastSquares(A, L);
    }

    // ============================

    // partial derivatives of parametrization coeffcients with respect to a constant
    // performed by numerical integration and least squares adjustment
    Matrix dacv_dconst;
    {
      logStatus<<"compute constant constraint"<<Log::endl;
      const UInt rows = static_cast<UInt>(std::round(2*PI/dAzimuth));
      const UInt cols = static_cast<UInt>(std::round(maxZenith/dZenith));

      Matrix A(rows*cols, parametrization->parameterCount());
      Vector l(rows*cols);
      UInt idx = 0;
      for(UInt i=0; i<rows; i++)
        for(UInt k=0; k<cols; k++)
        {
          const Angle azimuth  (i*2*PI/rows);
          const Angle elevation(PI/2-(k+0.5)*Double(dZenith));
          copy(parametrization->designMatrix(azimuth, elevation), A.row(idx));
          l(idx) = 1;                                   // constant
          A.row(idx) *= std::sqrt(std::cos(elevation)); // area weights
          l.row(idx) *= std::sqrt(std::cos(elevation));
          idx++;
        }
      dacv_dconst = leastSquares(A, l);
    }

    // ============================

    logStatus<<"apply constraints to antennas"<<Log::endl;
    const UInt parameterCount = parametrization->parameterCount();
    for(UInt idAnt=0; idAnt<antennaName.size(); idAnt++)
    {
      // antenna name
      std::string str;
      for(GnssType type : types.at(idAnt))
        str += " "+type.str();
      logInfo<<" "<<antennaName.at(idAnt)<<str<<Log::endl;

      // weight matrices of each pattern: taken from normal equations
      std::vector<Matrix> P(types.at(idAnt).size());
      for(UInt idType=0; idType<types.at(idAnt).size(); idType++)
      {
        P.at(idType) = normals.slice(index.at(idAnt).at(idType), index.at(idAnt).at(idType), parameterCount, parameterCount);
        fillSymmetric(P.at(idType));
        P.at(idType).setType(Matrix::GENERAL);
      }
      Matrix I = identityMatrix(parameterCount);

      // --- lambda ---------------------
      auto accumulateNormals = [&](Double weight, const std::vector<Matrix> &A)
      {
        // accumulate constraint normals
        for(UInt idType1=0; idType1<types.at(idAnt).size(); idType1++)
        {
          rankKUpdate(weight, A.at(idType1).trans(), normals.slice(index.at(idAnt).at(idType1), index.at(idAnt).at(idType1), parameterCount, parameterCount));
          for(UInt idType2=idType1+1; idType2<types.at(idAnt).size(); idType2++)
            matMult(weight, A.at(idType1).trans(), A.at(idType2), normals.slice(index.at(idAnt).at(idType1), index.at(idAnt).at(idType2), parameterCount, parameterCount));
        }
        info.observationCount += A.at(0).rows();
      };
      // --------------------------------

      for(const auto &constraint : constraints)
      {
        switch(constraint.type)
        {
          case Constraint::CENTER:
          {
            for(UInt idType=0; idType<types.at(idAnt).size(); idType++)
              if(types.at(idAnt).at(idType) == constraint.gnssType)
              {
                logInfo<<"  - center   "<<types.at(idAnt).at(idType).str()<<Log::endl;
                Matrix N = dacv_dcenter.trans() * ((constraint.applyWeight ? P.at(idType) : I) * dacv_dcenter);
                N.setType(Matrix::SYMMETRIC);
                cholesky(N);
                Matrix A = dacv_dcenter.trans() * (constraint.applyWeight ? P.at(idType) : I);
                triangularSolve(1., N.trans(), A);
                triangularSolve(1., N,         A);
                rankKUpdate(1./std::pow(constraint.sigma, 2), A, normals.slice(index.at(idAnt).at(idType), index.at(idAnt).at(idType), parameterCount, parameterCount));
                info.observationCount += 3;
              }
            break;
          }
          case Constraint::CENTERMEAN:
          {
            logInfo<<"  - zero mean center"<<Log::endl;
            // normals: estimate antenna center
            Matrix N(3, Matrix::SYMMETRIC);
            for(UInt idType=0; idType<types.at(idAnt).size(); idType++)
              N += dacv_dcenter.trans() * ((constraint.applyWeight ? P.at(idType) : I) * dacv_dcenter);
            cholesky(N);

            std::vector<Matrix> A(types.at(idAnt).size());
            for(UInt idType=0; idType<types.at(idAnt).size(); idType++)
            {
              A.at(idType) = (constraint.applyWeight ? (dacv_dcenter.trans()*P.at(idType)) : dacv_dcenter.trans());
              triangularSolve(1., N.trans(), A.at(idType));
              triangularSolve(1., N,         A.at(idType));
            }

            accumulateNormals(1./std::pow(constraint.sigma, 2), A);
            break;
          }
          case Constraint::CONSTANT:
          {
            for(UInt idType=0; idType<types.at(idAnt).size(); idType++)
              if(types.at(idAnt).at(idType) == constraint.gnssType)
              {
                logInfo<<"  - constant "<<types.at(idAnt).at(idType).str()<<Log::endl;
                // normals: estimate antenna constant
                Double N = inner(dacv_dconst, (constraint.applyWeight ? (P.at(idType)*dacv_dconst) : dacv_dconst));
                Matrix A = (1./N) * (constraint.applyWeight ? (dacv_dconst.trans()*P.at(idType)) : dacv_dconst.trans());
                // accumulate constraint normals
                rankKUpdate(1./std::pow(constraint.sigma, 2), A, normals.slice(index.at(idAnt).at(idType), index.at(idAnt).at(idType), parameterCount, parameterCount));
                info.observationCount += 1;
              }
            break;
          }
          case Constraint::CONSTANTMEAN:
          {
            // normals: estimate antenna center
            Double N = 0;
            for(UInt idType=0; idType<types.at(idAnt).size(); idType++)
              N += inner(dacv_dconst, (constraint.applyWeight ? (P.at(idType)*dacv_dconst) : dacv_dconst));
            std::vector<Matrix> A(types.at(idAnt).size());
            for(UInt idType=0; idType<types.at(idAnt).size(); idType++)
              A.at(idType) = (1./N) * (constraint.applyWeight ? (dacv_dconst.trans()*P.at(idType)) : dacv_dconst.trans());
            accumulateNormals(1./std::pow(constraint.sigma, 2), A);
            break;
          }
          case Constraint::TEC:
          {
            logInfo<<"  - zero mean TEC"<<Log::endl;
            Matrix N(parametrization->parameterCount(), Matrix::SYMMETRIC);
            for(UInt idType=0; idType<types.at(idAnt).size(); idType++)
              axpy(std::pow(types.at(idAnt).at(idType).ionosphericFactor(), 2), (constraint.applyWeight ? P.at(idType) : I), N);
            cholesky(N);

            // constraint as pseudo observation equations
            std::vector<Matrix> A(types.at(idAnt).size());
            for(UInt idType=0; idType<types.at(idAnt).size(); idType++)
            {
              A.at(idType) = types.at(idAnt).at(idType).ionosphericFactor() * (constraint.applyWeight ? P.at(idType) : I);
              triangularSolve(1., N.trans(), A.at(idType));
              triangularSolve(1., N,         A.at(idType));
            }

            accumulateNormals(1./std::pow(constraint.sigma, 2), A);
            break;
          }
        }
      }
    } // for(idAnt)

    // ============================

    logStatus<<"write constrained normal equations to <"<<fileNameNormalsOut<<">"<<Log::endl;
    writeFileNormalEquation(fileNameNormalsOut, info, normals, rhs);

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
