/***********************************************/
/**
* @file griddedData2PotentialCoefficients.cpp
*
* @brief Estimate potential coefficients from gridded gravity field functionals.
* With quadrature formular or
* least squares adjustment with block diagonal normals (order by order).
*
* @author Annette Eicker
* @author Torsten Mayer-Guerr
* @date 2005-05-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimate potential coefficients from \configFile{inputfileGriddedData}{griddedData}
gravity field functionals. It used a simple quadrature formular
\begin{equation}
  c_{nm} = \frac{1}{4\pi}\frac{R}{GM} \sum_i f_i \left(\frac{r_i}{R}\right)^{n+1} k_n C_{nm}(\lambda_i,\vartheta_i)\,\Delta\Phi_i
\end{equation}
or a \config{leastSquares} adjustment with block diagonal normal matrix (order by order).
For the latter one the data must be regular distributed.

The \config{value}s $f_i$ and the \config{weight}s $\Delta\Phi_i$ are expressions
using the common data variables for grids, see \reference{dataVariables}{general.parser:dataVariables}.
Multiple \configFile{outputfilePotentialCoefficients}{potentialCoefficients} can be estimated in one step.
For each an indivdual \config{value} must be specified.
The type of the gridded data (e.g gravity anomalies or geoid heights)
must be set with \configClass{kernel}{kernelType} $k_n$.

The expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusively. The coefficients are related
to the reference radius~\config{R} and the Earth gravitational constant \config{GM}.

For irregular distributed data and using the full variance covariance matrix use
\program{NormalsSolverVCE} together with \configClass{oberservation:terrestrial}{observationType:terrestrial}
and \configClass{parametrizationGravity:sphericalHarmonics}{parametrizationGravityType:sphericalHarmonics}.

See also \program{GriddedDataTimeSeries2PotentialCoefficients}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileSphericalHarmonics.h"
#include "files/fileGriddedData.h"
#include "classes/kernel/kernel.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Estimate potential coefficients from gridded gravity field functionals.
* With quadrature formular or
* least squares adjustment with block diagonal normals (order by order).
* @ingroup programsGroup */
class GriddedData2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedData2PotentialCoefficients, PARALLEL, "Estimate potential coefficients from gridded gravity field functionals", Grid, PotentialCoefficients)

/***********************************************/

void GriddedData2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    std::vector<FileName> fileNameOut;
    FileName              fileNameGrid;
    std::vector<ExpressionVariablePtr> exprValue;
    ExpressionVariablePtr exprArea;
    KernelPtr             kernel;
    Double                GM, R;
    UInt                  minDegree, maxDegree;
    Bool                  useLeastSquares;

    readConfig(config, "outputfilePotentialCoefficients", fileNameOut,     Config::MUSTSET,  "",      "one file for each value expression");
    readConfig(config, "inputfileGriddedData",            fileNameGrid,    Config::MUSTSET,  "",      "");
    readConfig(config, "value",                           exprValue,       Config::MUSTSET,  "data0", "expression to compute values (input columns are named data0, data1, ...)");
    readConfig(config, "weight",                          exprArea,        Config::MUSTSET,  "area",  "expression to compute values (input columns are named data0, data1, ...)");
    readConfig(config, "kernel",                          kernel,          Config::MUSTSET,  "",      "data type of input values");
    readConfig(config, "minDegree",                       minDegree,       Config::DEFAULT,  "0",     "");
    readConfig(config, "maxDegree",                       maxDegree,       Config::MUSTSET,  "",      "");
    readConfig(config, "GM",                              GM,              Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                               R,               Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "leastSquares",                    useLeastSquares, Config::DEFAULT,  "0",     "false: quadrature formular, true: least squares adjustment order by order");
    if(isCreateSchema(config)) return;

    if(fileNameOut.size() != exprValue.size())
      throw(Exception("number of outputfilePotentialCoefficients must agree with number of value expressions"));

    // reading grids
    // -------------
    logStatus<<"read grid from file <"<<fileNameGrid<<">"<<Log::endl;
    GriddedData grid;
    readFileGriddedData(fileNameGrid, grid);
    MiscGriddedData::printStatistics(grid);

    // evaluate expressions
    // --------------------
    VariableList varList;
    addDataVariables(grid, varList);
    std::for_each(exprValue.begin(), exprValue.end(), [&](auto expr) {expr->simplify(varList);});
    exprArea ->simplify(varList);
    std::vector<std::vector<Double>> values(exprValue.size(), std::vector<Double>(grid.points.size()));
    for(UInt i=0; i<grid.points.size(); i++)
    {
      evaluateDataVariables(grid, i, varList);
      for(UInt k=0; k<exprValue.size(); k++)
        values.at(k).at(i) = exprValue.at(k)->evaluate(varList);
      grid.areas.at(i) = exprArea->evaluate(varList);
    }
    grid.values = values;

    // spherical harmonic analysis
    // ---------------------------
    std::vector<SphericalHarmonics> harmonics = MiscGriddedData::analysisSphericalHarmonics(grid, kernel, minDegree, maxDegree, GM, R, useLeastSquares, comm);

    // write potential coefficients
    // ----------------------------
    if(Parallel::isMaster(comm))
      for(UInt k=0; k<harmonics.size(); k++)
      {
        logStatus<<"write potential coefficients to file <"<<fileNameOut.at(k)<<">"<<Log::endl;
        writeFileSphericalHarmonics(fileNameOut.at(k), harmonics.at(k));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
