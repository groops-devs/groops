/***********************************************/
/**
* @file instrumentCovarianceCheck.cpp
*
* @brief Remove non invertible covriance matrices.
*
* @author Norbert Zehentner
* @date 2016-01-20
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program checks \configFile{inputfileCovariance3d}{instrument}
3x3 covariance matrices if they are invertible or not and removes the invalid epochs.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Remove non invertible covariance matrices.
* @ingroup programsGroup */
class InstrumentCovarianceCheck
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentCovarianceCheck, SINGLEPROCESS, "Remove non invertible covariance matrices", Instrument, Covariance)

/***********************************************/

void InstrumentCovarianceCheck::run(Config &config)
{
  try
  {
    FileName fileNameOut, fileNameIn;

    readConfig(config, "outputfileCovariance3d", fileNameOut, Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileCovariance3d",  fileNameIn,  Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    logStatus<<"read and check covariance <"<<fileNameIn<<">"<<Log::endl;
    InstrumentFile file(fileNameIn);
    std::vector<Covariance3dArc> arcCov(file.arcCount());
    UInt removed = 0;
    Parallel::forEach(arcCov, [&](UInt arcNo)
    {
      Covariance3dArc arc = file.readArc(arcNo);
      Covariance3dArc arcNew;
      for(UInt i=0; i<arc.size(); i++)
      {
        try
        {
          Matrix W = arc.at(i).covariance.matrix();
          cholesky(W);
          arcNew.push_back(arc.at(i));
        }
        catch(std::exception &/*e*/)
        {
          removed++;
        }
      }
      return arcNew;
    });
    Parallel::reduceSum(removed);
    logInfo<<"  "<<removed<<" epochs removed"<<Log::endl;

    if(Parallel::isMaster())
    {
      logStatus<<"write covariance to file <"<<fileNameOut<<"> "<<Log::endl;
      InstrumentFile::write(fileNameOut, arcCov);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
