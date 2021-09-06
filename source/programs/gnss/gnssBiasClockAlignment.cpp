/***********************************************/
/**
* @file gnssBiasClockAlignment.cpp
*
* @brief Align GNSS transmitter clocks to reference clocks and adjust related receiver signal biases as well as GLONASS transmitter biases.
*
* @author Sebastian Strasser
* @date 2019-09-03
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program can be used to absolutely align GNSS transmitter clocks to reference clocks (i.e. broadcast clocks).
Each 'group' of \config{transmitter}s, usually a system like GPS or Galileo, is aligned individually by a constant shift over all transmitters.
If \config{alignClocksByFreqNo} is set, GLONASS transmitters will be divided by frequency number into groups of nominally two transmitters.
The offset between clocks and reference clocks will be shifted into receiver code biases, if \config{receiver} is provided."

By setting \config{alignFreqNoBiasesAtReceiver} and providing \config{receiver}, this program can further align GLONASS transmitter signal
biases so that the differences between frequency number-dependent receiver signal biases are minimal, which helps if PPP users don't set
up individual signal biases per frequency number at the receiver. Alignment is done by computing signal bias residuals to the mean over all
frequency numbers of a signal type at each receiver and then computing the means over all receivers for each frequency number and shifting
those from the receiver signal biases to the transmitter signal biases. Internal consistency of the biases is not affected by this.

If you only want to align GLONASS frequency numbers, provide the same clocks in
\configFile{inputfileClock}{instrument} and \configFile{inputfileReferenceClock}{instrument}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Align GNSS transmitter clocks to reference clocks and adjust related receiver signal biases as well as GLONASS transmitter biases.
* @ingroup programsGroup */
class GnssBiasClockAlignment
{
public:
  class Transmitter
  {
  public:
    FileName outNameClock, outNameBias, inNameClock, inNameReferenceClock, inNameBias, inNameTransmitterInfo;
    MiscValueArc clockArc, referenceClockArc;
    GnssSignalBias biases;
    GnssStationInfo info;
  };

  class Receiver
  {
  public:
    FileName outNameBias, inNameBias;
    GnssSignalBias biases;
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssBiasClockAlignment, SINGLEPROCESS, "Align GNSS transmitter clocks to reference clocks and adjust related receiver signal biases as well as GLONASS transmitter biases.", Gnss, Instrument)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssBiasClockAlignment::Transmitter &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "outputfileClock",          var.outNameClock,          Config::MUSTSET,  "", "aligned clock instrument file");
  readConfig(config, "outputfileSignalBias",     var.outNameBias,           Config::OPTIONAL, "", "(GLONASS only) aligned signal bias file");
  readConfig(config, "inputfileClock",           var.inNameClock,           Config::MUSTSET,  "", "clock instrument file");
  readConfig(config, "inputfileReferenceClock",  var.inNameReferenceClock,  Config::MUSTSET,  "", "reference clock instrument file");
  readConfig(config, "inputfileSignalBias",      var.inNameBias,            Config::OPTIONAL, "", "(GLONASS only) signal bias file");
  readConfig(config, "inputfileTransmitterInfo", var.inNameTransmitterInfo, Config::MUSTSET,  "{groopsDataDir}/gnss/transmitter/transmitterInfo/igs/igs14/transmitterInfo_igs14.{prn}.xml", "transmitter info file");
  endSequence(config);
  return TRUE;
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssBiasClockAlignment::Receiver &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "outputfileSignalBias",    var.outNameBias,          Config::MUSTSET, "", "aligned signal bias file");
  readConfig(config, "inputfileSignalBias",     var.inNameBias,           Config::MUSTSET, "", "signal bias file");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void GnssBiasClockAlignment::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    std::vector<Transmitter> transmitters;
    std::vector<Receiver> receivers;
    Bool alignClocksByFreqNo, alignFreqNoBiasesAtReceiver;

    readConfig(config, "transmitter",                 transmitters,                Config::MUSTSET,  "",  "one element per satellite");
    readConfig(config, "receiver",                    receivers,                   Config::OPTIONAL, "",  "one element per station");
    readConfig(config, "alignClocksByFreqNo",         alignClocksByFreqNo,         Config::DEFAULT,  "1", "align clocks for each GLONASS frequency number separately");
    readConfig(config, "alignFreqNoBiasesAtReceiver", alignFreqNoBiasesAtReceiver, Config::DEFAULT,  "1", "align frequency number-dependent code biases for each receiver");
    if(isCreateSchema(config)) return;

    logStatus<<"read transmitter input files"<<Log::endl;
    std::map<GnssType, std::vector<UInt>> group2DataIds;
    Single::forEach(transmitters.size(), [&](UInt i)
    {
      Transmitter *trans = &transmitters.at(i);
      try
      {
        trans->clockArc          = InstrumentFile::read(trans->inNameClock);
        trans->referenceClockArc = InstrumentFile::read(trans->inNameReferenceClock);
        if(!trans->inNameBias.empty())
          readFileGnssSignalBias(trans->inNameBias, trans->biases);
        readFileGnssStationInfo(trans->inNameTransmitterInfo, trans->info);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<" continue..."<<Log::endl;
        return;
      }

      UInt idRecv = trans->info.findReceiver(trans->clockArc.at(trans->clockArc.size()/2).time);
      if(idRecv == NULLINDEX)
      {
        logWarning<<"no receiver info found in <"<<trans->inNameTransmitterInfo<<">, skipping satellite..."<<Log::endl;
        return;
      }

      GnssReceiverInfo info = trans->info.receiver.at(idRecv);
      GnssType group("***"s+info.serial.at(0));
      if(alignClocksByFreqNo && group == GnssType::GLONASS && !info.version.empty())
        group.setFrequencyNumber(std::stoi(info.version));
      group2DataIds[group].push_back(i);
    });

    if(receivers.size())
    {
      logStatus<<"read receiver input files"<<Log::endl;
      Single::forEach(receivers.size(), [&](UInt i)
      {
        try
        {
          readFileGnssSignalBias(receivers.at(i).inNameBias, receivers.at(i).biases);
        }
        catch(std::exception &e)
        {
          logWarning<<e.what()<<" continue..."<<Log::endl;
        }
      });
    }

    // shift mean of transmitter clock diffs per group to receiver code biases
    logInfo<<"clock shifts per group:"<<Log::endl;
    for(const auto &group : group2DataIds)
    {
      // mean over all satellite clocks in group
      std::vector<Double> clockDiffs;
      for(UInt i : group.second)
        for(UInt idEpoch = 0; idEpoch < transmitters.at(i).clockArc.size(); idEpoch++)
          clockDiffs.push_back(transmitters.at(i).clockArc.at(idEpoch).value - transmitters.at(i).referenceClockArc.at(idEpoch).value);
      const Double mean = std::accumulate(clockDiffs.begin(), clockDiffs.end(), 0.0)/clockDiffs.size();
      logInfo<<group.first.str()<<group.second.size()%" (%i sat)"s<<" = "<<(mean*1e9)%"%10.3f ns"s<<" = "<<(mean*LIGHT_VELOCITY)%"%10.3f m"s<<Log::endl;

      // update transmitter clocks
      for(UInt i : group.second)
        for(UInt idEpoch = 0; idEpoch < transmitters.at(i).clockArc.size(); idEpoch++)
          transmitters.at(i).clockArc.at(idEpoch).value -= mean;

      // update receiver code biases
      for(auto &&rec : receivers)
        for(UInt i = 0; i < rec.biases.types.size(); i++)
          if(rec.biases.types.at(i) == group.first + GnssType::RANGE)
            rec.biases.biases.at(i) -= mean * LIGHT_VELOCITY;
    }

    if(alignFreqNoBiasesAtReceiver)
    {
      if(!receivers.size())
        throw(Exception("no receivers found, cannot align frequency number-dependent biases"));

      std::map<GnssType, std::vector<Double>> bias2Residuals;

      // compute residuals of receiver biases to mean bias over all GLONASS frequency numbers per station
      for(auto &&rec : receivers)
      {
        // GLONASS signal bias types independent of frequency number
        std::set<GnssType> types;
        for(const auto &type : rec.biases.types)
          if(type == GnssType::GLONASS)
            types.insert(type & ~GnssType::FREQ_NO & (type == GnssType::PHASE ? ~GnssType::ATTRIBUTE : GnssType::ALL));

        for(const auto &type : types)
        {
          // mean over all frequency numbers for this type
          std::vector<Double> values;
          for(UInt i = 0; i < rec.biases.types.size(); i++)
            if(rec.biases.types.at(i) == type)
              values.push_back(rec.biases.biases.at(i));
          const Double mean = ::mean(Vector(values));

          // bias residuals per frequency number for this type
          for(UInt i = 0; i < rec.biases.types.size(); i++)
            if(rec.biases.types.at(i) == type)
            {
              GnssType t = rec.biases.types.at(i) & (type == GnssType::PHASE ? ~GnssType::ATTRIBUTE : GnssType::ALL);
              bias2Residuals[t].push_back(rec.biases.biases.at(i) - mean);
            }
        }
      }

      // shift mean of receiver bias residuals over all stations to transmitter bias per type and frequency number
      logInfo<<"bias shifts per frequency number and type:"<<Log::endl;
      for(const auto &bias : bias2Residuals)
      {
        const GnssType type = bias.first;
        const Double mean = ::mean(Vector(bias.second));
        logInfo<<" "<<type.str()<<" = "<<mean%"%10.3f m"s<<Log::endl;

        // remove mean from receiver bias of all stations
        for(auto &&rec : receivers)
          for(UInt i = 0; i < rec.biases.types.size(); i++)
            if(rec.biases.types.at(i) == type)
              rec.biases.biases.at(i) -= mean;

        // add mean to transmitter bias
        for(auto &trans : transmitters)
          for(UInt i = 0; i < trans.biases.types.size(); i++)
            if(trans.biases.types.at(i) == type)
              trans.biases.biases.at(i) += mean;
      }
    }

    logStatus<<"writing aligned clock files <"<<transmitters.at(0).outNameClock<<">"<<Log::endl;
    for(const auto &trans : transmitters)
      if(trans.clockArc.size())
        InstrumentFile::write(trans.outNameClock, trans.clockArc);

    if(receivers.size())
    {
      logStatus<<"writing aligned receiver signal bias files <"<<receivers.at(0).outNameBias<<">"<<Log::endl;
      for(const auto &rec : receivers)
        if(rec.biases.types.size())
          writeFileGnssSignalBias(rec.outNameBias, rec.biases);

      if(alignFreqNoBiasesAtReceiver)
      {
        logStatus<<"writing aligned transmitter signal bias files <"<<transmitters.at(0).outNameBias<<">"<<Log::endl;
        for(const auto &trans : transmitters)
          if(trans.outNameBias)
            writeFileGnssSignalBias(trans.outNameBias, trans.biases);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
