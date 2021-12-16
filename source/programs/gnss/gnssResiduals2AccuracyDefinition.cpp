/***********************************************/
/**
* @file gnssResiduals2AccuracyDefinition.cpp
*
* @brief Compute antenna accuracies from observation residuals.
*
* @author Torsten Mayer-Guerr
* @date 2012-11-19
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Compute antenna accuracies from observation \configFile{inputfileResiduals}{instrument}.
The \configFile{inputfileStationInfo}{gnssStationInfo} is needed to assign
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

Example: Analysis of TerraSAR-X residuals of one month shows that low elevation
GPS satellites are not tracked by the onboard receiver. An estimation of accuracies
for these directions is not possible from the residuals and the apriori accuracies
are left untouched. The other directions show very low phase noise hardly elevation
and azimuth dependent for L2W. A nearly zero mean indicates the use of adequate antennca
center variations in the processing.

\fig{!hb}{0.8}{gnssResiduals2AccuracyDefinition}{fig:gnssResiduals2AccuracyDefinition}{L2W accuracies of TerraSAR-X determined from residuals of one month}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileGnssStationInfo.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief Compute antenna definition from observation residuals.
* @ingroup programsGroup */
class GnssResiduals2AccuracyDefinition
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssResiduals2AccuracyDefinition, SINGLEPROCESS, "Compute accuracy definition from observation residuals", Gnss)
GROOPS_RENAMED_PROGRAM(GnssResiduals2AntennaDefinition, GnssResiduals2AccuracyDefinition, date2time(2020, 6, 26))

/***********************************************/

void GnssResiduals2AccuracyDefinition::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameAntennaMean, fileNameAntennaAccuracy, fileNameAntennaRedundancy;
    FileName              fileNameStationInfo, fileNameAntenna;
    std::vector<FileName> fileNameResiduals;
    Bool                  isTransmitter;
    Double                thresholdOutlier, minRedundancy;

    renameDeprecatedConfig(config, "outputfileAntennaDefinition", "outputfileAntennaMean", date2time(2020, 7, 4));

    readConfig(config, "outputfileAccuracyDefinition", fileNameAntennaAccuracy,   Config::OPTIONAL, "",   "elevation and azimuth dependent accuracy");
    readConfig(config, "outputfileAntennaMean",        fileNameAntennaMean,       Config::OPTIONAL, "",   "weighted mean of the residuals");
    readConfig(config, "outputfileAntennaRedundancy",  fileNameAntennaRedundancy, Config::OPTIONAL, "",   "redundancy of adjustment");
    readConfig(config, "inputfileAccuracyDefinition",  fileNameAntenna,           Config::MUSTSET,  "",   "apriori accuracies");
    readConfig(config, "inputfileStationInfo",         fileNameStationInfo,       Config::MUSTSET,  "",   "to assign residuals to antennas");
    readConfig(config, "isTransmitter",                isTransmitter,             Config::DEFAULT,  "0",  "stationInfo is of a transmitter");
    readConfig(config, "thresholdOutlier",             thresholdOutlier,          Config::DEFAULT,  "10", "ignore residuals with sigma/sigma0 greater than threshold");
    readConfig(config, "minRedundancy",                minRedundancy,             Config::DEFAULT,  "3",  "min number of residuals. to estimate sigma");
    readConfig(config, "inputfileResiduals",           fileNameResiduals,         Config::MUSTSET,  "",   "GNSS receiver residuals");
    if(isCreateSchema(config)) return;

    // ============================

    GnssStationInfo stationInfo;
    if(!fileNameStationInfo.empty())
      readFileGnssStationInfo(fileNameStationInfo, stationInfo);

    GnssType typePRN;
    if(isTransmitter)
      typePRN = GnssType("***"+stationInfo.markerNumber);

    std::vector<GnssAntennaDefinitionPtr> antennaList;
    if(!fileNameAntenna.empty())
    {
      readFileGnssAntennaDefinition(fileNameAntenna, antennaList);
      stationInfo.fillAntennaPattern(antennaList);
    }

    // ============================

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
          // find antenna for epoch
          const UInt idAnt = stationInfo.findAntenna(epoch.time);
          if(idAnt == NULLINDEX)
            continue;
          GnssAntennaDefinitionPtr antenna = stationInfo.antenna.at(idAnt).antennaDef;
          if(!antenna)
            throw(Exception(epoch.time.dateTimeStr()+": antenna not found: "+stationInfo.antenna.at(idAnt).str()));

          UInt idObs = 0;
          for(GnssType satType : epoch.satellite)
          {
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

            Double azimuth   = epoch.observation.at(idObs+0);
            Double elevation = epoch.observation.at(idObs+1);
            if(isTransmitter)
            {
              azimuth   = epoch.observation.at(idObs+2);
              elevation = epoch.observation.at(idObs+3);
            }

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

              if(value)
                for(GnssAntennaPattern &pattern : antenna->pattern)
                  if(type+satType == pattern.type+typePRN)
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
                    if((redundancy > 0) && (sigma > 0) && (sigma <= thresholdOutlier))
                    {
                      const Double p = 1./std::pow(sigma, 2); // weight
                      pattern.ePe(idxL,idxB)        += std::pow(value, 2);
                      pattern.redundancy(idxL,idxB) += redundancy;
                      pattern.sum(idxL,idxB)        += p * value;
                      pattern.count(idxL,idxB)      += p;
                    }

                    // observations?
                    if(std::isnan(redundancy) || (redundancy <= 0))
                    {
                      pattern.ePe(idxL,idxB)        += std::pow(value, 2);
                      pattern.redundancy(idxL,idxB) += 1;
                      pattern.sum(idxL,idxB)        += value;
                      pattern.count(idxL,idxB)      += 1;
                    }

                    break;
                  }
            } // while()
          } // for(satType)
        } // for(epoch)
      } // for(arcNo)
    } // for(idFile)

    // ============================

    // only one value at zenith
    for(auto &antenna : antennaList)
      for(auto &pattern : antenna->pattern)
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
        for(auto &pattern : antenna->pattern)
          if(pattern.count.size())
          {
            pattern.offset = Vector3d();
            for(UInt i=0; i<pattern.pattern.rows(); i++)
              for(UInt k=0; k<pattern.pattern.columns(); k++)
                if(pattern.redundancy(i, k) >= minRedundancy)
                  pattern.pattern(i, k) = std::sqrt(pattern.ePe(i, k)/pattern.redundancy(i, k));
          }
      writeFileGnssAntennaDefinition(fileNameAntennaAccuracy, antennaList);
    }

    for(auto &antenna : antennaList)
      for(auto &pattern : antenna->pattern)
      {
        pattern.offset   = Vector3d();
        pattern.pattern *= NAN_EXPR;
      }

    if(!fileNameAntennaMean.empty())
    {
      logStatus<<"write antenna definition <"<<fileNameAntennaMean<<">"<<Log::endl;
      for(auto &antenna : antennaList)
        for(auto &pattern : antenna->pattern)
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
        for(auto &pattern : antenna->pattern)
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
