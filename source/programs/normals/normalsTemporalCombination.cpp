/***********************************************/
/**
* @file normalsTemporalCombination.cpp
*
* @brief Normal equations of trend, annual from a time series of normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2016-03-03
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads a times series of \configFile{inputfileNormalequation}{normalEquation}
with asscociated \configClass{timeSeries}{timeSeriesType} and setup a new combined normal equation system.
For each parameter a \configClass{parametrizationTemporal}{parametrizationTemporalType} is used.

It can be used to estimate trend and annual spherical harmonic coefficients from monthly GRACE normal equations.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileNormalEquation.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Normal equations of trend, annual from a time series of normal equations.
* @ingroup programsGroup */
class NormalsTemporalCombination
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsTemporalCombination, PARALLEL, "normal equations of trend, annual from a time series of normal equations", NormalEquation)

/***********************************************/

void NormalsTemporalCombination::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName                   fileNameOut;
    std::vector<FileName>      fileNameIn;
    ParametrizationTemporalPtr temporal;
    TimeSeriesPtr              timeSeries;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "inputfileNormalequation",  "inputfileNormalEquation",  date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "temporalRepresentation", "parametrizationTemporal", date2time(2020, 6, 3));

    readConfig(config, "outputfileNormalEquation", fileNameOut, Config::MUSTSET, "", "");
    readConfig(config, "inputfileNormalEquation",  fileNameIn,  Config::MUSTSET, "", "normal equations for each point in time");
    readConfig(config, "timeSeries",               timeSeries,  Config::MUSTSET, "", "times of each normal equations");
    readConfig(config, "parametrizationTemporal",  temporal,    Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    // =============================================

    // Create time series
    // ------------------
    std::vector<Time> times  = timeSeries->times();
    if(times.size() != fileNameIn.size())
      throw(Exception("fileCount("+fileNameIn.size()%"%i) != timeCount("s+times.size()%"%i)-1"s));

    // =============================================

    // read normal equations
    // ---------------------
    logStatus<<"read normal equations"<<Log::endl;
    MatrixDistributed  normal;
    Matrix             nTotal;
    NormalEquationInfo infoTotal;

    Single::forEach(times.size(), [&](UInt idEpoch)
    {
      // read normals
      // ------------
      Matrix N, n;
      NormalEquationInfo info;
      try
      {
        readFileNormalEquation(fileNameIn.at(idEpoch), info, N, n);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<" continue..."<<Log::endl;
        return;
      }

      // init normals
      // ------------
      const Vector factorTemporal = temporal->factors(times.at(idEpoch));
      if(normal.parameterCount()==0)
      {
        std::vector<UInt> blockIndex(factorTemporal.size()+1);
        for(UInt i=0; i<=factorTemporal.rows(); i++)
          blockIndex.at(i) = i*N.rows();
        normal.initEmpty(blockIndex, comm);
        nTotal = Matrix(normal.parameterCount(), n.columns());
        infoTotal.lPl = Vector(n.columns());
      }

      // test dimension
      // --------------
      if((n.rows() != normal.blockSize(0)) || (n.columns() != nTotal.columns()))
        throw(Exception("Normal equation dimension mismatch"));

      // parameter names
      // ---------------
      if(!infoTotal.parameterName.size())
        infoTotal.parameterName = info.parameterName;
      else
        for(UInt i=0; i<infoTotal.parameterName.size(); i++)
          if(!infoTotal.parameterName.at(i).combine(info.parameterName.at(i)))
            logWarning<<"Parameter names do not match at index "<<i<<": '"<<infoTotal.parameterName.at(i).str()<<"' != '"<< info.parameterName.at(i).str()<<"'"<< Log::endl;

      // setup normal matrix
      // -------------------
      fillSymmetric(N);
      for(UInt i=0; i<normal.blockCount(); i++)
        for(UInt k=i; k<normal.blockCount(); k++)
          if(factorTemporal(i) != 0 && factorTemporal(k) != 0)
          {
            normal.setBlock(i, k);

            if(normal.isMyRank(i,k))
              axpy(factorTemporal(i)*factorTemporal(k), N, normal.N(i,k));
          }

      if(Parallel::isMaster(comm))
      {
        for(UInt i=0; i<normal.blockCount(); i++)
          axpy(factorTemporal(i), n, nTotal.row(normal.blockIndex(i), normal.blockSize(i)));
        infoTotal.observationCount += info.observationCount;
        infoTotal.lPl += info.lPl;
      }
    });

    // =============================================

    // adjust parameter names based on temporal representation
    // -------------------------------------------------------
    std::vector<ParameterName> baseName;
    std::swap(infoTotal.parameterName, baseName);
    temporal->parameterName(baseName, infoTotal.parameterName);

    // write normal equations
    // ----------------------
    logStatus<<"write normal equations to <"<<fileNameOut<<">"<<Log::endl;
    writeFileNormalEquation(fileNameOut, infoTotal, normal, nTotal);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
