/***********************************************/
/**
* @file instrumentGnssReceiver2TimeSeries.cpp
*
* @brief Convert gnssReceiver file into instrument(MISCVLAUES) file.
*
* @author Torsten Mayer-Guerr
* @date 2020-04-09
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert selected GNSS observations or residuals into a simpler time series format.
The \config{outputfileTimeSeries} is an \file{instrument file}{instrument} (MISCVALUES).
For each epoch the first data column contains the PRN, the second the satellite system,
followed by a column for each GNSS \configClass{type}{gnssType}.
As normally more than one GNSS transmitter is tracked per epoch, the output file
has several lines per observed epoch (epochs with the same time, one for each transmitter).

The second data column of the output contains a number representating the system
\begin{itemize}
\item 71: 'G', GPS
\item 82: 'R', GLONASS
\item 69: 'E', GALILEO
\item 67: 'C', BDS
\item 83: 'S', SBAS
\item 74: 'J', QZSS
\item 73: 'I', IRNSS .
\end{itemize}

A \file{GNSS residual file}{instrument} includes additional information
besides the residuals, which can also be selected with \configClass{type}{gnssType}
\begin{itemize}
\item \verb|A1*|, \verb|E1*|: azimuth and elevation at receiver
\item \verb|A2*|, \verb|E2*|: azimuth and elevation at transmitter
\item \verb|I**|: Estimated slant total electron content (STEC)
\end{itemize}

Furthermore these files may include for each residual \configClass{type}{gnssType}
information about the redundancy and the accuracy relation $\sigma/\sigma_0$
of the estimated $\sigma$ versus the apriori $\sigma_0$ from the least squares adjustment.
The three values (residuals, redundancy, $\sigma/\sigma_0$) are coded with the same type.
To get acess to all values the corresponding type must be repeated in \configClass{type}{gnssType}.

Example: Selected GPS phase residuals (\configClass{type}{gnssType}='\verb|L1*G|' and \configClass{type}{gnssType}='\verb|L2*G|').
Plotted with \program{PlotGraph} with two \configClass{layer:linesAndPoints}{plotGraphLayerType}
(\config{valueX}='\verb|data0|',  \config{valueY}='\verb|100*data3+data1|' and \config{valueY}='\verb|100*data4+data1|' respectively).
\fig{!hb}{0.8}{instrumentGnssReceiver2TimeSeries}{fig:instrumentGnssReceiver2TimeSeries}{GPS residuals in cm, shifted by PRN}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Convert gnssReceiver file into instrument(MISCVLAUES) file.
* @ingroup programsGroup */
class InstrumentGnssReceiver2TimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentGnssReceiver2TimeSeries, SINGLEPROCESS, "Convert gnssReceiver file into instrument(MISCVLAUES) file", Gnss, TimeSeries, Residuals)

/***********************************************/

void InstrumentGnssReceiver2TimeSeries::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNamesIn;
    std::vector<GnssType> types;

    readConfig(config, "outputfileTimeSeries",  fileNameOut, Config::MUSTSET, "",  "Instrument (MISCVALUES): prn, system, values for each type");
    readConfig(config, "inputfileGnssReceiver", fileNamesIn, Config::MUSTSET, "",  "GNSS receiver observations or residuals");
    readConfig(config, "type",                  types,       Config::MUSTSET, "",  "");
    if(isCreateSchema(config)) return;

    // ============================

    std::vector<MiscValuesArc> arcList;
    for(auto &fileName : fileNamesIn)
    {
      logStatus<<"read GNSS receiver data <"<<fileName<<">"<<Log::endl;
      try
      {
        InstrumentFile fileReceiver(fileName);
        for(UInt arcNo=0; arcNo<fileReceiver.arcCount(); arcNo++)
        {
          GnssReceiverArc arc = fileReceiver.readArc(arcNo);
          MiscValuesArc arcNew;
          for(auto &epoch : arc)
          {
            UInt idObs = 0;
            for(GnssType typeSat : epoch.satellite)
            {
              // find first type for the satellite system
              UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), typeSat));

              MiscValuesEpoch epochNew(2+types.size()); // prn, system, types
              epochNew.time      = epoch.time;
              epochNew.values(0) = static_cast<Double>(typeSat.prn());
                   if(typeSat == GnssType::GPS)     epochNew.values(1) = static_cast<Double>('G');
              else if(typeSat == GnssType::GLONASS) epochNew.values(1) = static_cast<Double>('R');
              else if(typeSat == GnssType::GALILEO) epochNew.values(1) = static_cast<Double>('E');
              else if(typeSat == GnssType::BDS)     epochNew.values(1) = static_cast<Double>('C');
              else if(typeSat == GnssType::SBAS)    epochNew.values(1) = static_cast<Double>('S');
              else if(typeSat == GnssType::QZSS)    epochNew.values(1) = static_cast<Double>('J');
              else if(typeSat == GnssType::IRNSS)   epochNew.values(1) = static_cast<Double>('I');

              // loop over all obs for this satellite
              Bool                  found    = FALSE;
              std::vector<GnssType> typesTmp = types;
              while((idType<epoch.obsType.size()) && (idObs<epoch.observation.size()) && (epoch.obsType.at(idType) == typeSat))
              {
                const GnssType type  = epoch.obsType.at(idType++);
                const Double   value = epoch.observation.at(idObs++);
                const UInt     idx   = GnssType::index(typesTmp, type);
                if((idx != NULLINDEX) && value)
                {
                  epochNew.values(2+idx) = value;
                  typesTmp.at(idx) = GnssType(static_cast<UInt64>(-1));
                  found = TRUE;
                }
              }

              if(found)
                arcNew.push_back(epochNew);
            } // for(idTrans)
          } // for(idEpoch)

          if(arcNew.size())
            arcList.push_back(arcNew);
        } // for(idArc)
      }
      catch(...)
      {
        logError<<"error by opening file, continue..."<<Log::endl;
      }
    } // for(fileName)

    // ============================

    if(!arcList.size())
      logWarning<<"empty arc"<<Log::endl;

    logStatus<<"write time series file <"<<fileNameOut<<">"<<Log::endl;
    InstrumentFile::write(fileNameOut, arcList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
