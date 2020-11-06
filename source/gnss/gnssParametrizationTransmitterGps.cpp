/***********************************************/
/**
* @file gnssParametrizationTransmitterGps.cpp
*
* @brief GPS satellites (transmitter).
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-04-27
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileStringTable.h"
#include "classes/timeSeries/timeSeries.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrizationTransmitter.h"
#include "gnss/gnssParametrizationTransmitterGps.h"

/***********************************************/

GnssParametrizationTransmitterGps::GnssParametrizationTransmitterGps(Config &config)
{
  try
  {
    FileName    fileNameTransmitterList, fileNameTransmitterInfo, fileNameAntennaDef, fileNameReceiverDef;
    UInt        interpolationDegree;
    std::string choice;

    readConfig(config, "inputfileTransmitterList",      fileNameTransmitterList,       Config::MUSTSET,  "{groopsDataDir}/gnss/transmitterGps/transmitterList.gps.txt", "ascii file with transmitter PRNs, used to loop variable {prn}");
    readConfig(config, "inputfileTransmitterInfo",      fileNameTransmitterInfo,       Config::MUSTSET,  "{groopsDataDir}/gnss/transmitterGps/transmitterInfo/igs/igs14/transmitterInfo_igs14.{prn}.xml", "variable {prn} available");
    readConfig(config, "inputfileAntennaDefintion",     fileNameAntennaDef,            Config::MUSTSET,  "{groopsDataDir}/gnss/transmitterGps/antennaDefinition/igs/igs14/antennaDefinition_igs14.xml", "phase centers and variations (ANTEX)");
    readConfig(config, "inputfileReceiverDefintion",    fileNameReceiverDef,           Config::OPTIONAL, "{groopsDataDir}/gnss/transmitterGps/receiverDefinition/receiverDefinition.xml", "transmitted signal types");
    readConfig(config, "inputfileSignalBias",           fileNameTemplateInSignal,      Config::OPTIONAL, "signalBias_{loopTime:%D}.{prn}.txt",  "variable {prn} available");
    readConfig(config, "inputfileVariational",          fileNameTemplateInVariational, Config::MUSTSET,  "variational_{loopTime:%D}.{prn}.dat", "variable {prn} available");
    readConfig(config, "inputfileClock",                fileNameTemplateInClock,       Config::MUSTSET,  "clock_{loopTime:%D}.{prn}.dat",       "variable {prn} available");
    readConfig(config, "outputfileSignalBias",          fileNameTemplateOutSignal,     Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "outputfileUsedTransmitterList", fileNameOutTransmitterList,    Config::OPTIONAL, "", "ascii file with PRNs");
    if(readConfigChoice(config, "estimateOrbit", choice, Config::OPTIONAL, "", ""))
    {
      if(readConfigChoiceElement(config, "dynamic", choice, "initial state and model parameters are estimated (in CRF)"))
      {
        estimateDynamicOrbit = TRUE;
        TimeSeriesPtr stochasticPulsePtr;

        renameDeprecatedConfig(config, "parameter", "parametrizationAcceleration", date2time(2020, 6, 3));

        readConfig(config, "outputfileVariational",       fileNameTemplateOutVariational, Config::OPTIONAL, "", "variable {prn} available");
        readConfig(config, "outputfileOrbit",             fileNameTemplateOutOrbit,       Config::OPTIONAL, "", "variable {prn} available");
        readConfig(config, "outputfileParameters",        fileNameTemplateOutParameter,   Config::OPTIONAL, "", "variable {prn} available");
        readConfig(config, "parametrizationAcceleration", parametrizationAcceleration,    Config::DEFAULT,  "", "orbit force parameters");
        if(readConfigSequence(config, "stochasticPulse", Config::OPTIONAL, "", "[mum/s] parametrization of stochastic pulses"))
        {
          readConfig(config, "timeSeries",    stochasticPulsePtr,           Config::MUSTSET, "",  "");
          readConfig(config, "onlyEclipsing", onlyEclipsingStochasticPulse, Config::DEFAULT, "0", "only apply stochastic pulses to eclipsing satellites");
          endSequence(config);
        }
        readConfig(config, "minEstimableEpochsRatio", minEstimableEpochsRatio, Config::MUSTSET, "0.75", "drop satellites with lower ratio of estimable epochs to total epochs");

        if(!isCreateSchema(config))
          if(stochasticPulsePtr)
            stochasticPulse = stochasticPulsePtr->times();
      }
      endChoice(config);
    }
    readConfig(config, "estimateCodeBias",           estimateCodeBias,   Config::DEFAULT,  "0", "");
    readConfig(config, "estimatePhaseBias",          estimatePhaseBias,  Config::DEFAULT,  "1", "");
    readConfig(config, "supportsIntegerAmbiguities", integerAmbiguities, Config::DEFAULT,  "1", "receiver tracks full cycle integer ambiguities");
    readConfig(config, "biasModel",                  biasModel,          Config::OPTIONAL, "",  "estimate time-variable signal bias");
    readConfig(config, "timeVariableBias",           timeVariableBiases, Config::OPTIONAL, "",  "a priori time-variable signal bias");
    if(readConfigChoice(config, "estimateClockError", choice, Config::OPTIONAL, "", ""))
    {
      if(readConfigChoiceElement(config, "epochWise", choice, "clock error is estimated for each epoch"))
      {
        estimateClockError = EstimateClockError::EPOCH;
        readConfig(config, "outputfileClock",         fileNameOutClock,   Config::OPTIONAL, "",       "variable {prn} available");
        readConfig(config, "sigmaZeroMeanConstraint", sigmaClockZeroMean, Config::DEFAULT,  "0.0001", "(0 = unconstrained) sigma [m] for zero-mean constraint over all satellite clocks");
      }
      if(readConfigChoiceElement(config, "noiseModel", choice, "clock errors are estimated based on a noise model"))
      {
        estimateClockError = EstimateClockError::NOISEMODEL;
        readConfig(config, "outputfileClock",          fileNameOutClock,   Config::OPTIONAL,  "",  "variable {prn} available");
        readConfig(config, "inputfileClockNoiseModel", fileNameClockNoise, Config::MUSTSET,   "",  "variable {prn} available");
        readConfig(config, "estimateVarianceFactor",   estimateClockSigma, Config::DEFAULT,   "1", "via Variance Component Estimation");
        readConfig(config, "clockModel",               clockModel,         Config::OPTIONAL,  "",  "");
        readConfig(config, "sigmaZeroMeanConstraint",  sigmaClockZeroMean, Config::DEFAULT,   "0.0001", "(0 = unconstrained) sigma [m] for zero-mean constraint over all satellite clocks");
      }
      endChoice(config);
    }
    readConfig(config, "antennaCenter",                   antennaCenterVariations,         Config::OPTIONAL, "",  "estimate antenna center variations");
    readConfig(config, "integrationDegree",               integrationDegree,               Config::MUSTSET,  "7", "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",             interpolationDegree,             Config::MUSTSET,  "7", "for orbit interpolation and velocity calculation");
    readConfig(config, "disableShadowEpochs",             disableShadowEpochs,             Config::DEFAULT,  "0", "disable epochs if satellite is in Earth's/Moon's shadow");
    readConfig(config, "disablePostShadowRecoveryEpochs", disablePostShadowRecoveryEpochs, Config::DEFAULT,  "1", "disable epochs if satellite is in post-shadow recovery maneuver (e.g. GPS block IIA)");
    readConfig(config, "ephemerides",                     ephemerides,                     Config::MUSTSET,  "", "");
    readConfig(config, "eclipse",                         eclipse,                         Config::MUSTSET,  "", "eclipse model used to determine if a satellite is in Earth's shadow");
    if(readConfigChoice(config, "noAntennaPatternFound", choice, Config::MUSTSET, "useNearestFrequency", "what should happen is no antenna pattern is found for an observation"))
    {
      if(readConfigChoiceElement(config, "ignoreObservation",   choice, "ignore observation if no matching pattern is found"))
        noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::IGNORE_OBSERVATION;
      if(readConfigChoiceElement(config, "useNearestFrequency", choice, "use pattern of nearest frequency if no matching pattern is found"))
        noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::USE_NEAREST_FREQUENCY;
      if(readConfigChoiceElement(config, "throwException",      choice, "throw exception if no matching pattern is found"))
        noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::THROW_EXCEPTION;
      endChoice(config);
    }
    if(isCreateSchema(config)) return;

    init(fileNameTransmitterList, fileNameTransmitterInfo, fileNameAntennaDef, fileNameReceiverDef, interpolationDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
