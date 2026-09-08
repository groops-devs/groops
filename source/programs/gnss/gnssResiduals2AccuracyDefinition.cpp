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
The \configFile{inputfileStationInfo}{platform} is needed to assign
the residuals to the equipped antenna at observation times.

The \configFile{outputfileAccuracyDefinition}{gnssAntennaDefinition} contains
at first step the same accuracy information for all antennas as the input file.
Only the azimuth~$A$ and elevation~$E$ dependent grid points of the patterns
where enough residuals are available ($>$ \config{minRedundancy})
are replaced by the estimated robust weighted accuracy
\begin{equation}
  \sigma(A,E) = \sqrt{\alpha\frac{\sum_i e_i^2 p_i}{\sum_i r_i \bar{p}}},
\end{equation}
where $e_i$ are the azimuth and elevation dependent residuals and $r_i$ the
corresponding redundancies (number of observations minus the contribution to
the estimated parameters). The weights $p_i=1/\sigma_i^2$ are computed from the
estimated accuracies in the least squares adjustment. Since large outliers are downweighted
using \config{huber} and \config{huberPower}, a correction factor $\alpha$ must be applied to
obtain an unbiased estimate of $\sigma$.

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
and azimuth dependent for L2W. A nearly zero mean indicates the use of adequate antenna
center variations in the processing.

\fig{!hb}{0.8}{gnssResiduals2AccuracyDefinition}{fig:gnssResiduals2AccuracyDefinition}{L2W accuracies of TerraSAR-X determined from residuals of one month}

See also \program{GnssResiduals2TransmitterAccuracyDefinition}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/filePlatform.h"
#include "inputOutput/system.h"
#include "misc/varianceComponentEstimation.h"

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
    Double                minRedundancy;
    Double                huber, huberPower;

    renameDeprecatedConfig(config, "outputfileAntennaDefinition", "outputfileAntennaMean", date2time(2020, 7, 4));

    readConfig(config, "outputfileAccuracyDefinition", fileNameAntennaAccuracy,   Config::OPTIONAL, "",    "elevation and azimuth dependent accuracy");
    readConfig(config, "outputfileAntennaMean",        fileNameAntennaMean,       Config::OPTIONAL, "",    "weighted mean of the residuals");
    readConfig(config, "outputfileAntennaRedundancy",  fileNameAntennaRedundancy, Config::OPTIONAL, "",    "redundancy of adjustment");
    readConfig(config, "inputfileAccuracyDefinition",  fileNameAntenna,           Config::MUSTSET,  "",    "apriori accuracies");
    readConfig(config, "inputfileStationInfo",         fileNameStationInfo,       Config::MUSTSET,  "",    "to assign residuals to antennas");
    readConfig(config, "isTransmitter",                isTransmitter,             Config::DEFAULT,  "0",   "stationInfo is of a transmitter");
    readConfig(config, "minRedundancy",                minRedundancy,             Config::DEFAULT,  "3",   "min. redundancy of residuals to estimate sigma");
    readConfig(config, "huber",                        huber,                     Config::DEFAULT,  "2.5", "residuals > huber*sigma0 have been downweighted");
    readConfig(config, "huberPower",                   huberPower,                Config::DEFAULT,  "1.5", "residuals > huber: sigma=(e/huber)^huberPower*sigma0");
    readConfig(config, "inputfileResiduals",           fileNameResiduals,         Config::MUSTSET,  "",    "GNSS receiver residuals");
    if(isCreateSchema(config)) return;

    // ============================

    Platform stationInfo;
    if(!fileNameStationInfo.empty())
      readFilePlatform(fileNameStationInfo, stationInfo);

    GnssType typePRN;
    if(isTransmitter)
      typePRN = GnssType("***"+stationInfo.markerNumber);

    std::vector<GnssAntennaDefinitionPtr> antennaList;
    if(!fileNameAntenna.empty())
    {
      readFileGnssAntennaDefinition(fileNameAntenna, antennaList);
      stationInfo.fillGnssAntennaDefinition(antennaList);
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

      for(auto &epoch : GnssReceiverArc(InstrumentFile::read(fileName)))
      {
        // find antenna for epoch
        auto ant= stationInfo.findEquipment<PlatformGnssAntenna>(epoch.time);
        if(!ant)
          continue;
        GnssAntennaDefinitionPtr antenna = ant->antennaDef;
        if(!antenna)
          throw(Exception(epoch.time.dateTimeStr()+": antenna not found: "+ant->str()));

        UInt idObs = 0;
        for(GnssType satType : epoch.satellite)
        {
          Double azimuth=NAN_EXPR, elevation=NAN_EXPR;

          // find type for the satellite system, loop over all obs for this satellite
          UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), satType));
          while((idType<epoch.obsType.size()) && (idObs<epoch.observation.size()) && (epoch.obsType.at(idType) == satType))
          {
            GnssType type  = epoch.obsType.at(idType++) + satType;
            Double   value = epoch.observation.at(idObs++);

            if(isTransmitter) // isTransmitter
            {
              if(type == (GnssType::AZIMUT    + GnssType::L2)) {azimuth   = value; continue;}
              if(type == (GnssType::ELEVATION + GnssType::L2)) {elevation = value; continue;}
            }
            else
            {
              if(type == (GnssType::AZIMUT    + GnssType::L1)) {azimuth   = value; continue;}
              if(type == (GnssType::ELEVATION + GnssType::L1)) {elevation = value; continue;}
            }

            Double redundancy=NAN_EXPR, sigma=NAN_EXPR;
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

            if(!value || std::isnan(value))
              continue;

            for(GnssAntennaPattern &pattern : antenna->patterns)
              if(type+satType == pattern.type+typePRN)
              {
                if(std::isnan(azimuth) || std::isnan(elevation))
                  throw(Exception("file must contain azimuth and elevation."));
                const UInt idxL = static_cast<UInt>(std::round((Double(azimuth)+2*PI)/(2*PI)*pattern.pattern.rows()))%pattern.pattern.rows();
                const UInt idxB = static_cast<UInt>(std::round((PI/2-Double(elevation))/Double(pattern.dZenit)));
                if(idxB >= pattern.pattern.columns())
                  break;

                if(!pattern.count.size())
                {
                  pattern.ePe        = Matrix(pattern.pattern.rows(), pattern.pattern.columns());
                  pattern.redundancy = Matrix(pattern.pattern.rows(), pattern.pattern.columns());
                  pattern.sum        = Matrix(pattern.pattern.rows(), pattern.pattern.columns());
                  pattern.weight     = Matrix(pattern.pattern.rows(), pattern.pattern.columns());
                  pattern.count      = Matrix(pattern.pattern.rows(), pattern.pattern.columns());
                }

                // residuals?
                if((redundancy > 0) && (sigma > 0))
                {
                  const Double p = 1./std::pow(sigma, 2); // weight
                  pattern.ePe(idxL,idxB)        += p * std::pow(value, 2);
                  pattern.redundancy(idxL,idxB) += redundancy;
                  pattern.sum(idxL,idxB)        += p * value;
                  pattern.weight(idxL,idxB)     += p;
                  pattern.count(idxL,idxB)      += 1;
                }

                break;
              }
          } // while()
        } // for(satType)
      } // for(epoch)
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
          copy(Vector(pattern.pattern.rows(), sum(pattern.weight    .column(0))), pattern.weight    .column(0));
        }

    // ============================

    for(auto &antenna : antennaList)
      for(auto &pattern : antenna->patterns)
        pattern.offset = Vector3d();

    if(!fileNameAntennaAccuracy.empty())
    {
      // unbiased estimation of sigmas needs to consider huber downweighting
      // -------------------------------------------------------------------
      // numerical integration of normal distribution
      constexpr Double dx = 1e-4;
      Double x      = dx/2;
      Double factor = 0;
      for(; x<std::min(huber, 10.); x+=dx)
        factor += std::exp(-0.5*x*x) * dx;
      // variance of downweighted normal distribution
      for(; x<10.; x+=dx)
        factor += std::pow(x/huber, -2*huberPower) * std::exp(-0.5*x*x) * dx;
      factor *= 2./std::sqrt(2*PI);

      logStatus<<"write accuracy definition <"<<fileNameAntennaAccuracy<<">"<<Log::endl;
      for(auto &antenna : antennaList)
        for(auto &pattern : antenna->patterns)
          if(pattern.count.size())
            for(UInt i=0; i<pattern.pattern.rows(); i++)
              for(UInt k=0; k<pattern.pattern.columns(); k++)
                if(pattern.redundancy(i, k) >= minRedundancy)
                  pattern.pattern(i, k) = std::sqrt(factor*pattern.count(i, k)/pattern.weight(i, k))
                                        * Vce::standardDeviation(pattern.ePe(i, k), pattern.redundancy(i, k), huber, huberPower);
      writeFileGnssAntennaDefinition(fileNameAntennaAccuracy, antennaList);
    }

    if(!fileNameAntennaMean.empty())
    {
      logStatus<<"write antenna definition <"<<fileNameAntennaMean<<">"<<Log::endl;
      for(auto &antenna : antennaList)
        for(auto &pattern : antenna->patterns)
        {
          pattern.pattern *= NAN_EXPR;
          if(pattern.count.size())
            for(UInt i=0; i<pattern.pattern.rows(); i++)
              for(UInt k=0; k<pattern.pattern.columns(); k++)
                pattern.pattern(i, k) = pattern.sum(i, k)/pattern.weight(i, k);
        }
      writeFileGnssAntennaDefinition(fileNameAntennaMean, antennaList);
    }

    if(!fileNameAntennaRedundancy.empty())
    {
      logStatus<<"write redundancy <"<<fileNameAntennaRedundancy<<">"<<Log::endl;
      for(auto &antenna : antennaList)
        for(auto &pattern : antenna->patterns)
        {
          pattern.pattern *= NAN_EXPR;
          if(pattern.count.size())
            for(UInt i=0; i<pattern.pattern.rows(); i++)
              for(UInt k=0; k<pattern.pattern.columns(); k++)
                pattern.pattern(i, k) = pattern.redundancy(i, k);
        }
      writeFileGnssAntennaDefinition(fileNameAntennaRedundancy, antennaList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
