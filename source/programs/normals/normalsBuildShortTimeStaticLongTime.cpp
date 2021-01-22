/***********************************************/
/**
* @file normalsBuildShortTimeStaticLongTime.cpp
*
* @brief Normal equations with static, short time, and long time gravity field parameters.
*
* @author Torsten Mayer-Guerr
* @date 2012-08-13
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program sets up normal equations based on \configClass{observation}{observationType}.
Additionally short time and long time variations can be parametrized based on the static parameters
in \configClass{observation}{observationType} in an efficient way. The observation equations
are divided into time intervals $i \in \{1, ..., N\}$ (e.g. daily) as defined in
\configFile{inputfileArcList}{arcList}.

With \config{estimateLongTimeVariations} additional temporal variations can be co-estimated
for a subset of the parameters selected by \configClass{parameterSelection}{parameterSelectorType}.
These parameters might be spherical harmonic coefficients with a limited maximum degree.
The temporal variations are represented by base functions $\Phi_k(t_i)$ (e.g. trend and annual oscillation)
given by \configClass{parametrizationTemporal}{parametrizationTemporalType}.
The temporal base functions are evaluated at the mid time~$t_i$ of each interval~$i$, multiplicated
with the design matrix $\M A_i$ of the selected parameters, and the design matrix is extended
accordingly.

\fig{!hb}{0.8}{normalsBuildShortTimeStaticLongTime}{fig:normalsBuildShortTimeStaticLongTime}{Schema of the extended design matrix.}

With \config{estimateShortTimeVariations} short time variations of the gravity field can be co-estimated.
Their purpose is to mitigate temporal aliasing.
The short time parameters selected by \configClass{parameterSelection}{parameterSelectorType}
(e.g. daily constant or linear splines every 6 hour) are constrained by an
\configClass{autoregressiveModelSequence}{autoregressiveModelSequenceType}. If only a static parameter
set is selected the coressponding part of the design matrix is copied and modeled as a constant value
per interval in \configFile{inputfileArcList}{arcList} additionally so the corresponding temporal factor can be expressed as
\begin{equation}
  \Phi_i(t)  =
  \begin{cases}
    1 &\text{if} \hspace{5pt} t \in [t_i, t_{i+1}) \\
    0 & \text{otherwise}
  \end{cases}.
\end{equation}

Before writing the normal equations to \configFile{outputfileNormalEquation}{normalEquation}
short time gravity and satellite specific parameters can be eliminated with \config{eliminateParameter}.

Example: For the computation of the mean gravity field ITSG-Grace2018s with additional trend and annual signal
the normal equations are computed month by month and accumulated afterwards (see \program{NormalsAccumulate}).
The observations were divided into daily intervals with \configFile{inputfileArcList}{arcList}.
The static gravity field has been parametrized as spherical harmonics
up to degree $n=200$ in \configClass{observation:parametrizationGravity}{parametrizationGravityType}.
The trend and annual signals defined by
\configClass{estimateLongTimeVariations:parametrizationTemporal}{parametrizationTemporalType}
were estimated for selected parameters up to degree $n=120$.
To mitigate temporal aliasing daily gravity fields up to degree $n=40$ were setup and constrained
with an \configClass{autoregressiveModelSequence}{autoregressiveModelSequenceType} up to order three.

A detailed description of the approach is given in:
Kvas, A., Mayer-Gürr, T. GRACE gravity field recovery with background model uncertainties.
J Geod 93, 2543–2552 (2019). \url{https://doi.org/10.1007/s00190-019-01314-1}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileArcList.h"
#include "files/fileNormalEquation.h"
#include "classes/observation/observation.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "misc/kalmanProcessing.h"
#include "misc/normalsShortTimeStaticLongTime.h"

/***** CLASS ***********************************/

/** @brief Normal equations with static, short time, and long time gravity field parameters.
* @ingroup programsGroup */
class NormalsBuildShortTimeStaticLongTime
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsBuildShortTimeStaticLongTime, PARALLEL, "Normal equations with static, short time, and long time gravity field parameters", NormalEquation)
GROOPS_RENAMED_PROGRAM(KalmanStaticTemporalNormals, NormalsBuildShortTimeStaticLongTime, date2time(2020, 12, 7))

/***********************************************/

void NormalsBuildShortTimeStaticLongTime::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName                       fileNameNormals;
    ObservationPtr                 observation;
    AutoregressiveModelSequencePtr arSequence;
    ParameterSelectorPtr           parameterShortTime;
    ParametrizationTemporalPtr     temporalLongTime;
    ParameterSelectorPtr           parameterLongTime;
    FileName                       fileNameArcList;
    UInt                           defaultBlockSize;
    Bool                           eliminateParameter;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "arcList",                  "inputfileArcList",         date2time(2020, 7, 7));

    readConfig(config, "outputfileNormalEquation", fileNameNormals, Config::MUSTSET, "", "outputfile for normal equations");
    readConfig(config, "observation",              observation,     Config::MUSTSET, "", "");
    if(readConfigSequence(config, "estimateShortTimeVariations", Config::OPTIONAL, "", "co-estimate short time gravity field variations"))
    {
      readConfig(config, "autoregressiveModelSequence", arSequence,         Config::MUSTSET, "", "AR model sequence for constraining short time gravity variations");
      readConfig(config, "parameterSelection",          parameterShortTime, Config::MUSTSET, "", "parameters describing the short time gravity field");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateLongTimeVariations", Config::OPTIONAL, "", "co-estimate long time gravity field variations"))
    {
      readConfig(config, "parametrizationTemporal", temporalLongTime,  Config::MUSTSET, "", "parametrization of time variations (trend, annual, ...)");
      readConfig(config, "parameterSelection",      parameterLongTime, Config::MUSTSET, "", "parameters describing the long time gravity field");
      endSequence(config);
    }
    readConfig(config, "inputfileArcList",   fileNameArcList,    Config::MUSTSET, "",     "list to correspond points of time to arc numbers");
    readConfig(config, "defaultBlockSize",   defaultBlockSize,   Config::DEFAULT, "2048", "block size for distributing the normal equations, 0: one block");
    readConfig(config, "eliminateParameter", eliminateParameter, Config::DEFAULT, "1",    "eliminate short time and state parameter");
    if(isCreateSchema(config)) return;

    // =======================

    logStatus<<"read arc list <"<<fileNameArcList<<">"<<Log::endl;
    std::vector<UInt> arcsInterval;
    std::vector<Time> timesInterval;
    readFileArcList(fileNameArcList, arcsInterval, timesInterval);

    // init normal equations
    // ---------------------
    logStatus<<"initialize normal equations"<<Log::endl;
    NormalsShortTimeStaticLongTime normals;
    normals.init(observation, timesInterval, defaultBlockSize, comm, TRUE/*sortStateBeforeGravityParameter*/,
                 (arSequence) ? arSequence->dimension() : 0, parameterShortTime,
                 temporalLongTime, parameterLongTime);
    normals.setBlocks(arcsInterval);

    // setup observation equations
    // ---------------------------
    logStatus<<"accumulate normals from observation equations"<<Log::endl;
    Parallel::forEachInterval(observation->arcCount(), arcsInterval, [&](UInt arcNo)
    {
      // search time interval
      UInt idInterval = 0;
      while(arcsInterval.at(idInterval+1) <= arcNo)
        idInterval++;
      observation->setInterval(timesInterval.at(idInterval), timesInterval.at(idInterval+1));

      // observation equations
      Matrix l, A, B;
      observation->observation(arcNo, l, A, B);
      if(l.rows() == 0)
        return;

      normals.accumulate(idInterval, l, A, B);
    }, comm);
    observation = nullptr;

    // collect system of normal equations
    // ----------------------------------
    logStatus<<"collect system of normal equations"<<Log::endl;
    normals.reduceSum();

    // add normals of short time model
    // -------------------------------
    if(arSequence)
    {
      logStatus<<"add normals of short time model"<<Log::endl;
      normals.addShortTimeNormals(1., arSequence->normalEquationSequence());
    }

    // eliminate interval & state parameters
    // -------------------------------------
    if(eliminateParameter)
    {
      logStatus<<"eliminate interval parameters from normal equations"<<Log::endl;
      normals.regularizeUnusedParameters(normals.blockIndexStatic);
      normals.cholesky(TRUE/*timing*/, 0, normals.blockIndexStatic, TRUE/*collect*/);
      normals.triangularTransSolve(normals.n, 0, normals.blockIndexStatic);
      if(Parallel::isMaster(comm))
      {
        normals.obsCount -= normals.blockIndex(normals.blockIndexStatic);
        // lPl = lPl - n2^T N2^(-1) n2
        for(UInt i=0; i<normals.lPl.rows(); i++)
          normals.lPl(i) -= quadsum(normals.n.slice(0, i, normals.blockIndex(normals.blockIndexStatic),1));
        // remove additional parameters
        normals.n = normals.n.row(normals.blockIndex(normals.blockIndexStatic), normals.parameterCount() - normals.blockIndex(normals.blockIndexStatic));
        normals.parameterNames.erase(normals.parameterNames.begin(), normals.parameterNames.begin()+normals.blockIndex(normals.blockIndexStatic));
      }
      normals.eraseBlocks(0, normals.blockIndexStatic);
    } // if(eliminateParameter)

    // Write normal equations
    // ----------------------
    logStatus<<"write normal equations to <"<<fileNameNormals<<">"<<Log::endl;
    writeFileNormalEquation(fileNameNormals, NormalEquationInfo(normals.parameterNames, normals.lPl, normals.obsCount), normals, normals.n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
