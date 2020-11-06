/***********************************************/
/**
* @file graceSequenceOfEvents.cpp
*
* @brief Read GRACE SOE file.
*
* @author Beate Klinger
* @author Torsten Mayer-Guerr
* @date 2015-01-09
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts the GRACE SOE (sequence of events) file/format into \file{instrument file (MISCVALUES)}{instrument}.
The GRACE SOE format is described in "GRACE Level 1B Data Product User Handbook JPL D-22027" and "TN-03\_SOE\_format.txt"
given at \url{http://podaac.jpl.nasa.gov/grace/documentation.html}.
The output is one arc of satellite data which can include data gaps.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read GRACE SOE file.
* @ingroup programsConversionGroup */
class GraceSequenceOfEvents
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceSequenceOfEvents, SINGLEPROCESS, "read GRACE SOE file", Conversion, Grace)

/***********************************************/

void GraceSequenceOfEvents::run(Config &config)
{
  try
  {
    FileName    fileNameGraceA, fileNameGraceB;
    FileName    fileName;
    std::string eventType, mode;
    Double      sampling = 0;

    std::string choice;
    readConfig(config, "outputfileGraceA", fileNameGraceA, Config::OPTIONAL, "", "");
    readConfig(config, "outputfileGraceB", fileNameGraceB, Config::OPTIONAL, "", "");
    readConfig(config, "inputfile",        fileName,       Config::MUSTSET,  "{groopsDataDir}/grace/TN-01_SOE.txt", "SoE file");
    if(readConfigChoice(config, "events", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "ACCT", choice, "DSHL HeaterDisconnect"))
      {
        eventType = "ACCT";
        std::string choice;
        if(readConfigChoice(config, "mode", choice, Config::OPTIONAL, "", ""))
        {
          if(readConfigChoiceElement(config, "Heater",   choice, "DSHL HeaterDisconnect")) {mode = "DSHL-HeaterDisconnect-Event";}
          if(readConfigChoiceElement(config, "SetPoint", choice, "temperature set point")) {mode = "SetPoint";}
          endChoice(config);
        }
      }
      if(readConfigChoiceElement(config, "AOCS", choice, "coarse pointing mode or attitude hold mode"))
      {
        eventType = "AOCS";
        std::string choice;
        if(readConfigChoice(config, "mode", choice, Config::OPTIONAL, "", ""))
        {
          if(readConfigChoiceElement(config, "CPM", choice, "coarse pointing mode")) {mode = "CPM";}
          if(readConfigChoiceElement(config, "AHM", choice, "attitude hold mode"))   {mode = "AHM";}
          if(readConfigChoiceElement(config, "SM",  choice, "science mode"))         {mode = "SM";}
          endChoice(config);
        }
      }
      if(readConfigChoiceElement(config, "ACCR",  choice, "ACCR")) {eventType = "ACCR";}
      if(readConfigChoiceElement(config, "CMCAL", choice, "CoM calibration maneuver"))
      {
        eventType = "CMCAL";
        readConfig(config, "sampling", sampling, Config::OPTIONAL, "", "[seconds] create events between start and end of maneuver");
      }
      if(readConfigChoiceElement(config, "KBRCAL", choice, "KBR calibration maneuver"))
      {
        eventType = "KBRCAL";
        readConfig(config, "sampling", sampling, Config::OPTIONAL, "", "[seconds] create events between start and end of maneuver");
      }
      if(readConfigChoiceElement(config, "VCM",   choice, "CoM coordinates in SRF (m)"))              {eventType = "VCM";}
      if(readConfigChoiceElement(config, "VKB",   choice, "KBR phase center coordinates in SRF (m)")) {eventType = "VKB";}
      if(readConfigChoiceElement(config, "ICUVP", choice, "ICUVP"))                                   {eventType = "ICUVP";}
      if(readConfigChoiceElement(config, "IPU",   choice, "IPU"))                                     {eventType = "IPU";}
      if(readConfigChoiceElement(config, "IPUR",  choice, "IPUR"))                                    {eventType = "IPUR";}
      if(readConfigChoiceElement(config, "KAMI",  choice, "KAMI: time tag offset to Ka-phase meas.")) {eventType = "KAMI";}
      if(readConfigChoiceElement(config, "KMI",   choice, "K_MI: time tag offset to K-phase meas."))  {eventType = "K_MI";}
      if(readConfigChoiceElement(config, "KTOFF", choice, "KTOFF: time tag offset to KBR meas."))     {eventType = "KTOFF";}
      if(readConfigChoiceElement(config, "MANV",  choice, "MANV"))                                    {eventType = "MANV";}
      if(readConfigChoiceElement(config, "MTE1",  choice, "MTE1"))                                    {eventType = "MTE1";}
      if(readConfigChoiceElement(config, "MTE2",  choice, "MTE2"))                                    {eventType = "MTE2";}
      if(readConfigChoiceElement(config, "OCC",   choice, "OCC"))                                     {eventType = "OCC";}
      if(readConfigChoiceElement(config, "QSA",   choice, "SCA to SRF frame rotation"))               {eventType = "QSA";}
      if(readConfigChoiceElement(config, "QKS",   choice, "SCA to KBR frame rotation"))               {eventType = "QKS";}
      endChoice(config);
    }
    if(isCreateSchema(config)) return;

    logStatus<<"read input file <"<<fileName<<">"<<Log::endl;
    MiscValuesArc arcGraceA;
    MiscValuesArc arcGraceB;

    InFile file(fileName);

    std::string line;
    while(std::getline(file, line))
    {
      if(line.empty())
        continue;
      if(line.at(0) == 'x')
        continue; // the line has been deleted from active use.

      std::stringstream ss(line);
      ss.exceptions(std::ios::badbit | std::ios::failbit);

      Double          seconds;
      std::string     spaceCraft;
      std::string     event;
      UInt            count;
      std::string     comment;
      ss>>seconds>>spaceCraft>>event>>count;
      MiscValuesEpoch epoch(count);
      epoch.time = mjd2time(51544.5) + seconds2time(seconds);
      for(UInt i=0; i<count; i++)
        ss>>epoch.values(i);
      std::getline(ss, comment);

      if(event != eventType)
        continue;

      if((!mode.empty()) && (comment.find(mode) == std::string::npos))
        continue;

      if((event == "CMCAL") || (event == "KBRCAL"))
        if((sampling > 0) && (epoch.values(0) == 0)) // maneuver stop
        {
          auto sample = [](const Time &timeEnd, const Time &sampling, MiscValuesArc &arc)
          {
            if(arc.size() == 0)
              return;
            for(;;)
            {
              MiscValuesEpoch epoch = arc.back();
              epoch.time += sampling;
              if(epoch.time > timeEnd)
                return;
              arc.push_back(epoch);
            }
          };

          if(spaceCraft == "GRACEA")
            sample(epoch.time, seconds2time(sampling), arcGraceA);
          else if(spaceCraft == "GRACEB")
            sample(epoch.time, seconds2time(sampling), arcGraceB);
          continue;
        }

      if(spaceCraft == "GRACEA")
        arcGraceA.push_back(epoch);
      else if(spaceCraft == "GRACEB")
        arcGraceB.push_back(epoch);
    } // for(line)

    logInfo<<"GRACEA: events = "<<arcGraceA.size()<<Log::endl;
    logInfo<<"GRACEB: events = "<<arcGraceB.size()<<Log::endl;


    // write data
    // ----------
    if(!fileNameGraceA.empty())
    {
      logInfo<<"write data to <"<<fileNameGraceA<<">"<<Log::endl;
      InstrumentFile::write(fileNameGraceA, arcGraceA);
    }

    if(!fileNameGraceB.empty())
    {
      logInfo<<"write data to <"<<fileNameGraceB<<">"<<Log::endl;
      InstrumentFile::write(fileNameGraceB, arcGraceB);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
