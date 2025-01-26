/***********************************************/
/**
* @file gnssResiduals2TransmitterAccuracyDefinition.cpp
*
* @brief Compute antenna accuracies from observation residuals.
*
* @author Torsten Mayer-Guerr
* @date 2024-08-08
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Compute antenna accuracies from observation \configFile{inputfileResiduals}{instrument}.
The \configFile{inputfileTransmitterInfo}{platform} is needed to assign
the residuals to the equipped antenna at observation times.

The \configFile{outputfileAccuracyDefinition}{gnssAntennaDefinition} contains
at first step the same accuracy information for all antennas as the input file.
Only the azimuth~$A$ and elevation~$E$ dependent grid points of the patterns
where enough residuals are available ($>$ \config{minRedundancy})
are replaced by estimated accuracy
\begin{equation}
 \sigma(A,E) = \sqrt{\frac{\sum_i e_i^2(A,E)}{\sum_i r_i(A,E)}},
\end{equation}
where $e_i$ are the azimuth and elevation dependent residuals and $r_i$ the
corresponding redundancies (number of observations minus the contribution to
the estimated parameters).

The \configFile{inputfileAccuracyDefinition}{gnssAntennaDefinition} can be modified
to the demands before with \program{GnssAntennaDefinitionCreate}
(e.g. with \config{antenna:resample}).

To verify the results the \configFile{outputfileAntennaMean}{gnssAntennaDefinition}
and the accumulated \configFile{outputfileAntennaRedundancy}{gnssAntennaDefinition}
of the computed pattern grid points can be written.

See also \program{GnssResiduals2AccuracyDefinition}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/filePlatform.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief Compute antenna definition from observation residuals.
* @ingroup programsGroup */
class GnssResiduals2TransmitterAccuracyDefinition
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssResiduals2TransmitterAccuracyDefinition, SINGLEPROCESS, "Compute accuracy definition from observation residuals", Gnss)

/***********************************************/

void GnssResiduals2TransmitterAccuracyDefinition::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameAntennaMean, fileNameAntennaAccuracy, fileNameAntennaRedundancy;
    FileName              fileNameTransmitterInfo, fileNameAntenna;
    std::vector<FileName> fileNameResiduals;
    Double                minRedundancy;

    readConfig(config, "outputfileAccuracyDefinition", fileNameAntennaAccuracy,   Config::OPTIONAL, "",   "elevation and azimuth dependent accuracy");
    readConfig(config, "outputfileAntennaMean",        fileNameAntennaMean,       Config::OPTIONAL, "",   "weighted mean of the residuals");
    readConfig(config, "outputfileAntennaRedundancy",  fileNameAntennaRedundancy, Config::OPTIONAL, "",   "redundancy of adjustment");
    readConfig(config, "inputfileAccuracyDefinition",  fileNameAntenna,           Config::MUSTSET,  "",   "apriori accuracies");
    readConfig(config, "inputfileTransmitterInfo",     fileNameTransmitterInfo,   Config::MUSTSET,  "",   "to assign residuals to antennas");
    readConfig(config, "minRedundancy",                minRedundancy,             Config::DEFAULT,  "3",  "min number of residuals. to estimate sigma");
    readConfig(config, "inputfileResiduals",           fileNameResiduals,         Config::MUSTSET,  "",   "GNSS receiver residuals");
    if(isCreateSchema(config)) return;

    std::vector<GnssAntennaDefinitionPtr> antennaList;
    readFileGnssAntennaDefinition(fileNameAntenna, antennaList);

    std::map<GnssType, Platform> platforms;
    for(const FileName &fileName : fileNameResiduals)
    {
      logStatus<<"read GNSS residuals <"<<fileName<<">"<<Log::endl;
      if(!System::exists(fileName))
      {
        logWarning<<"file not exist -> continue"<<Log::endl;
        continue;
      }

      InstrumentFile fileReceiver(fileName);
      for(UInt arcNo=0; arcNo<fileReceiver.arcCount(); arcNo++)
      {
        GnssReceiverArc arc = fileReceiver.readArc(arcNo);
        for(auto &epoch : arc)
        {
          UInt idObs = 0;
          for(GnssType satType : epoch.satellite)
          {
            if(platforms.find(satType) == platforms.end())
            {
              VariableList fileNameVariableList;
              fileNameVariableList.setVariable("prn", satType.prnStr());
              readFilePlatform(fileNameTransmitterInfo(fileNameVariableList), platforms[satType]);
              platforms[satType].fillGnssAntennaDefinition(antennaList);
            }
            // find antenna for epoch
            auto ant = platforms[satType].findEquipment<PlatformGnssAntenna>(epoch.time);
            if(!ant)
              continue;
            GnssAntennaDefinitionPtr antenna = ant->antennaDef;
            if(!antenna)
              continue;
              // throw(Exception(epoch.time.dateTimeStr()+": antenna not found: "+ant->str()));

            // find type for the satellite system
            UInt idType = 0;
            while(epoch.obsType.at(idType) != satType)
              idType++;

            // azimuth and elevation
            if((epoch.obsType.at(idType+0) != (GnssType::AZIMUT    + GnssType::L1)) ||
               (epoch.obsType.at(idType+1) != (GnssType::ELEVATION + GnssType::L1)) ||
               (epoch.obsType.at(idType+2) != (GnssType::AZIMUT    + GnssType::L2)) ||
               (epoch.obsType.at(idType+3) != (GnssType::ELEVATION + GnssType::L2)))
              throw(Exception("azimuth and elevation expected"));

            const Double azimuth   = epoch.observation.at(idObs+2); // transmitter
            const Double elevation = epoch.observation.at(idObs+3); // transmitter

            idObs  += 4;  // skip azimuth and elevation
            idType += 4;

            // resiudals, redundancy, sigma/sigma0
            while((idType<epoch.obsType.size()) && (idObs<epoch.observation.size()) && (epoch.obsType.at(idType) == satType))
            {
              GnssType type  = epoch.obsType.at(idType++);
              Double   value = epoch.observation.at(idObs++);

              Double redundancy=0, sigma=0;
              if((idType < epoch.obsType.size()) && (type == epoch.obsType.at(idType))) // next redundancy?
              {
                type       = epoch.obsType.at(idType++);
                redundancy = epoch.observation.at(idObs++);
              }
              if((idType < epoch.obsType.size()) && (type == epoch.obsType.at(idType))) // next sigma?
              {
                type  = epoch.obsType.at(idType++);
                sigma = epoch.observation.at(idObs++);
              }
              while((idType < epoch.obsType.size()) && (type == epoch.obsType.at(idType))) // other additional information?
                idObs++, idType++;

              if(!value)
                continue;


              for(GnssAntennaPattern &pattern : antenna->patterns)
                if(type+satType == pattern.type)
                {
                  const UInt idxL = static_cast<UInt>(std::round((Double(azimuth)+2*PI)/(2*PI)*pattern.pattern.rows()))%pattern.pattern.rows();
                  const UInt idxB = static_cast<UInt>(std::round((PI/2-Double(elevation))/Double(pattern.dZenit)));
                  if(idxB >= pattern.pattern.columns())
                    break;

                  if(!pattern.count.size())
                  {
                    pattern.sum        = Matrix(pattern.pattern.rows(), pattern.pattern.columns());
                    pattern.ePe        = Matrix(pattern.pattern.rows(), pattern.pattern.columns());
                    pattern.redundancy = Matrix(pattern.pattern.rows(), pattern.pattern.columns());
                    pattern.count      = Matrix(pattern.pattern.rows(), pattern.pattern.columns());
                  }

                  // residuals?
                  if((redundancy > 0) && (sigma > 0))
                  {
                    const Double p = 1./std::pow(sigma, 2); // weight
                    pattern.ePe(idxL,idxB)        += p * std::pow(value, 2);
                    pattern.redundancy(idxL,idxB) += redundancy;
                    pattern.sum(idxL,idxB)        += p * value;
                    pattern.count(idxL,idxB)      += p;
                  }
                }
            } // while()
          } // for(satType)
        } // for(epoch)
      } // for(arcNo)
    } // for(idFile)

    // ============================

    // only one value at zenith
    for(auto &antenna : antennaList)
      for(auto &pattern : antenna->patterns)
        if(pattern.count.size())
        {
          copy(Vector(pattern.pattern.rows(), sum(pattern.sum       .column(0))), pattern.sum       .column(0));
          copy(Vector(pattern.pattern.rows(), sum(pattern.ePe       .column(0))), pattern.ePe       .column(0));
          copy(Vector(pattern.pattern.rows(), sum(pattern.redundancy.column(0))), pattern.redundancy.column(0));
          copy(Vector(pattern.pattern.rows(), sum(pattern.count     .column(0))), pattern.count     .column(0));
        }

    // ============================


    if(!fileNameAntennaAccuracy.empty())
    {
      logStatus<<"write accuracy definition <"<<fileNameAntennaAccuracy<<">"<<Log::endl;
      for(auto &antenna : antennaList)
        for(auto &pattern : antenna->patterns)
          if(pattern.count.size())
          {
            pattern.offset = Vector3d();
            for(UInt i=0; i<pattern.pattern.rows(); i++)
              for(UInt k=0; k<pattern.pattern.columns(); k++)
                if(pattern.redundancy(i, k) >= minRedundancy)
                  pattern.pattern(i, k) = std::sqrt(pattern.ePe(i, k)/pattern.redundancy(i, k)/pattern.count(i, k));
          }
      writeFileGnssAntennaDefinition(fileNameAntennaAccuracy, antennaList);
    }

    for(auto &antenna : antennaList)
      for(auto &pattern : antenna->patterns)
      {
        pattern.offset   = Vector3d();
        pattern.pattern *= NAN_EXPR;
      }

    if(!fileNameAntennaMean.empty())
    {
      logStatus<<"write antenna definition <"<<fileNameAntennaMean<<">"<<Log::endl;
      for(auto &antenna : antennaList)
        for(auto &pattern : antenna->patterns)
          if(pattern.count.size())
            for(UInt i=0; i<pattern.pattern.rows(); i++)
              for(UInt k=0; k<pattern.pattern.columns(); k++)
                pattern.pattern(i, k) = pattern.sum(i, k)/pattern.count(i, k);
      writeFileGnssAntennaDefinition(fileNameAntennaMean, antennaList);
    }

    if(!fileNameAntennaRedundancy.empty())
    {
      logStatus<<"write redundancy <"<<fileNameAntennaRedundancy<<">"<<Log::endl;
      for(auto &antenna : antennaList)
        for(auto &pattern : antenna->patterns)
          if(pattern.count.size())
            for(UInt i=0; i<pattern.pattern.rows(); i++)
              for(UInt k=0; k<pattern.pattern.columns(); k++)
                if(pattern.redundancy(i, k) >= minRedundancy)
                  pattern.pattern(i, k) = pattern.redundancy(i, k);
      writeFileGnssAntennaDefinition(fileNameAntennaRedundancy, antennaList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
