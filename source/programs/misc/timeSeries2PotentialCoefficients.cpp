/***********************************************/
/**
* @file timeSeries2PotentialCoefficients.cpp
*
* @brief Write time series as potential coefficients for each epoch.
*
* @author Torsten Mayer-Guerr
* @date 2020-07-20
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Interpret the data columns of \configFile{inputfileTimeSeries}{instrument}
as potential coefficients. The sequence of coefficients is given by
\configClass{numbering}{sphericalHarmonicsNumberingType} starting from data column \config{startDataFields}.

For each epoch a \configFile{outputfilesPotentialCoefficients}{potentialCoefficients}
is written where the \config{variableLoopTime} and \config{variableLoopIndex} are expanded for
each point of the given time series to create the file name for this epoch,
see \reference{text parser}{general.parser:text}.

See also \program{Gravityfield2PotentialCoefficientsTimeSeries}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/sphericalHarmonics.h"
#include "files/fileInstrument.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Write time series as potential coefficients for each epoch.
* @ingroup programsGroup */
class TimeSeries2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(TimeSeries2PotentialCoefficients, SINGLEPROCESS, "write time series as potential coefficients for each epoch.", Misc, TimeSeries, PotentialCoefficients)

/***********************************************/

void TimeSeries2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameOut, fileNameInstrument;
    std::string nameTime, nameIndex, nameCount;
    UInt        startData;
    UInt        minDegree, maxDegree;
    Double      GM, R;
    SphericalHarmonicsNumberingPtr numbering;

    readConfig(config, "outputfilesPotentialCoefficients", fileNameOut,        Config::MUSTSET,  "coeff_{loopTime:%y-%m}.gfc", "for each epoch");
    readConfig(config, "variableLoopTime",                 nameTime,           Config::OPTIONAL, "loopTime",        "variable with time of each epoch");
    readConfig(config, "variableLoopIndex",                nameIndex,          Config::OPTIONAL, "",                "variable with index of current epoch (starts with zero)");
    readConfig(config, "variableLoopCount",                nameCount,          Config::OPTIONAL, "",                "variable with total number of epochs");
    readConfig(config, "inputfileTimeSeries",              fileNameInstrument, Config::MUSTSET,  "",                "each epoch: multiple data for points (MISCVALUES)");
    readConfig(config, "startDataFields",                  startData,          Config::DEFAULT,  "0",               "first data column");
    readConfig(config, "minDegree",                        minDegree,          Config::MUSTSET,  "0",               "minimal degree");
    readConfig(config, "maxDegree",                        maxDegree,          Config::MUSTSET,  "",                "maximal degree");
    readConfig(config, "GM",                               GM,                 Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                                R,                  Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "numbering",                        numbering,          Config::MUSTSET,  "",                "numbering scheme");
    if(isCreateSchema(config)) return;

    std::vector<std::vector<UInt>> idxC, idxS;
    numbering->numbering(maxDegree, minDegree, idxC, idxS);
    const UInt parameterCount = numbering->parameterCount(maxDegree, minDegree);

    logStatus<<"read time series <"<<fileNameInstrument<<">"<<Log::endl;
    MiscValuesArc arc = InstrumentFile::read(fileNameInstrument);
    Arc::printStatistics(arc);
    const UInt dataCount = Epoch::dataCount(arc.getType(), TRUE/*mustDefined*/);
    if(dataCount < startData + parameterCount)
      throw(Exception("file contains not enough data columns: "+dataCount%"%i < "s+startData%"%i + "s+parameterCount%"%i"s));

    // write data
    // ----------
    const std::vector<Time> times = arc.times();
    VariableList varList;
    if(!nameTime.empty())  varList.undefineVariable(nameTime);
    if(!nameIndex.empty()) varList.undefineVariable(nameIndex);
    if(!nameCount.empty()) varList.setVariable(nameCount, times.size());

    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      Matrix cnm(maxDegree+1, Matrix::SYMMETRIC, Matrix::LOWER);
      Matrix snm(maxDegree+1, Matrix::SYMMETRIC, Matrix::LOWER);
      for(UInt n=minDegree; n<=maxDegree; n++)
      {
        if(idxC[n][0] != NULLINDEX) cnm(n,0) = arc.at(idEpoch).values(startData+idxC[n][0]);
        for(UInt m=1; m<=n; m++)
        {
          if(idxC[n][m] != NULLINDEX) cnm(n,m) = arc.at(idEpoch).values(startData+idxC[n][m]);
          if(idxS[n][m] != NULLINDEX) snm(n,m) = arc.at(idEpoch).values(startData+idxS[n][m]);
        }
      }

      if(!nameTime.empty())  varList.setVariable(nameTime, times.at(idEpoch).mjd());
      if(!nameIndex.empty()) varList.setVariable(nameIndex, idEpoch);
      logStatus<<"write potential coefficients <"<<fileNameOut(varList)<<">"<<Log::endl;
      writeFileSphericalHarmonics(fileNameOut(varList), SphericalHarmonics(GM, R, cnm, snm));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
