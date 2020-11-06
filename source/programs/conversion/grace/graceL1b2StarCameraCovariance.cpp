/***********************************************/
/**
* @file graceL1b2StarCameraCovariance.cpp
*
* @brief Covariance matrix from star camera flags.
*
* @author Torsten Mayer-Guerr
* @date 2017-08-22
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes star camera covariance matrices (\file{instrument file, COVARIANE3D}{instrument})
for a GRACE satellite under consideration of the active camera heads and an a priori variance factor.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Covariance matrix from star camera flags.
* @ingroup programsConversionGroup */
class GraceL1b2StarCameraCovariance
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2StarCameraCovariance, PARALLEL, "Covariance matrix from star camera flags", Conversion, Grace, Covariance, Instrument)

/***********************************************/

void GraceL1b2StarCameraCovariance::run(Config &config)
{
  try
  {
    FileName fileNameCovariance;
    FileName fileNameStarCamera, fileNameFlags, fileNameSoeQSA;
    Double   sigma0;

    readConfig(config, "outputfileStarCameraCovariance", fileNameCovariance, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCameraFlags",       fileNameFlags,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileSequenceOfEventsQSA",   fileNameSoeQSA,     Config::OPTIONAL, "", "");
    readConfig(config, "sigma0",                         sigma0,             Config::DEFAULT,  "6", "[seconds of arc]");
    if(isCreateSchema(config)) return;

    // Covariance of one head
    Tensor3d cov0;
    cov0.xx() = std::pow(  sigma0*DEG2RAD/3600, 2);
    cov0.yy() = std::pow(  sigma0*DEG2RAD/3600, 2);
    cov0.zz() = std::pow(8*sigma0*DEG2RAD/3600, 2);
    Tensor3d cov1  = rotaryX(Angle(-135*DEG2RAD)).rotate(cov0);
    Tensor3d cov2  = rotaryX(Angle(+135*DEG2RAD)).rotate(cov0);
    Matrix   inv1  = cov1.matrix(); inverse(inv1);
    Matrix   inv2  = cov2.matrix(); inverse(inv2);
    Matrix   inv12 = inv1+inv2;     inverse(inv12);
    Tensor3d cov12(inv12);

    // star camera orientations
    MiscValuesArc qsa = InstrumentFile::read(fileNameSoeQSA);

    logStatus<<"computing covariance matrix"<<Log::endl;
    InstrumentFile flagsFile(fileNameFlags);

    std::vector<Arc> arcList(flagsFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      MiscValuesArc flags = flagsFile.readArc(arcNo);
      Covariance3dArc arc;

      for(UInt i=0; i<flags.size(); i++)
      {
        // update rotations
        while(qsa.size() && (qsa.front().time <= flags.at(i).time)) // Move to correct interval
        {
          const Rotary3d rot1(qsa.front().values.row(0,4)); // qsa is rotation SCA frame to SRF
          const Rotary3d rot2(qsa.front().values.row(4,4));
          cov1 = rot1.rotate(cov0);
          cov2 = rot2.rotate(cov0);
          Matrix inv1  = cov1.matrix(); inverse(inv1);
          Matrix inv2  = cov2.matrix(); inverse(inv2);
          Matrix inv12 = inv1+inv2;     inverse(inv12);
          cov12 = inv12;
          qsa.remove(0);
        }

        const Double factor = 1.;
        // factor *= std::pow(epochFlags.values(1), 2); // = q_res;

        Covariance3dEpoch epoch;
        epoch.time = flags.at(i).time;
        if(flags.at(i).values(0) == 1)      // SCA 1
          epoch.covariance = factor * cov1;
        else if(flags.at(i).values(0) == 2) // SCA 2
          epoch.covariance = factor * cov2;
        else if(flags.at(i).values(0) == 4) // SCA 1 + SCA 2
          epoch.covariance = factor * cov12;
        else
        {
//           logWarning<<epoch.time.dateTimeStr()<<" Unknown SCA identification number: "<<flags.at(i).values(0)<<Log::endl;
          epoch.covariance = factor * cov12;
        }
        arc.push_back(epoch);
      }
      return arc;
    });

    if(Parallel::isMaster())
    {
      logStatus<<"write covariance data to file <"<<fileNameCovariance<<">"<<Log::endl;
      InstrumentFile::write(fileNameCovariance, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
