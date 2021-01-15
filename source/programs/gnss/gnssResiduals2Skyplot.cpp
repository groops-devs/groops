/***********************************************/
/**
* @file gnssResiduals2Skyplot.cpp
*
* @brief Convert residuals into griddedData format for plotting.
*
* @author Torsten Mayer-Guerr
* @date 2013-07-12
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write GNSS residuals together with azimuth and elevation to be plotted with \program{PlotMap}.
Azimuth and elevation are written as ellipsoidal longitude and latitude in a \file{griddedData file}{griddedData}.
The choosen ellipsoid parameters \config{R} and \config{inverseFlattening} are arbitrary but should be the same
as in \program{PlotMap}. If with \configClass{typeTransmitter}{gnssType} (e.g. '\verb|***G18|')
a single transmitter is selected the azimuth and elevation are computed from the transmitter point of view.

For each GNSS \configClass{type}{gnssType} an extra data column is created.

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
The 3 values (residuals, redundancy, $\sigma/\sigma_0$) are coded with the same type.
To get acess to all values the coresponding type must be repeated in \configClass{type}{gnssType}.

\fig{!hb}{0.5}{gnssResiduals2Skyplot}{fig:gnssResiduals2Skyplot}{GPS C2W residuals of GRAZ station at 2012-01-01}
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/system.h"
#include "files/fileInstrument.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Convert residuals into griddedData format for plotting.
* @ingroup programsGroup */
class GnssResiduals2Skyplot
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssResiduals2Skyplot, SINGLEPROCESS, "Convert residuals into griddedData format for plotting", Gnss, Grid)
GROOPS_RENAMED_PROGRAM(GnssResiduals2GriddedData, GnssResiduals2Skyplot, date2time(2019, 9, 9))

/***********************************************/

void GnssResiduals2Skyplot::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameGriddedData;
    std::vector<FileName> fileNameResiduals;
    std::vector<GnssType> types;
    GnssType              typeTransmitter;
    Double                a, f;

    readConfig(config, "outputfileGriddedData", fileNameGriddedData, Config::MUSTSET,  "", "");
    readConfig(config, "type",                  types,               Config::MUSTSET,  "", "");
    readConfig(config, "typeTransmitter",       typeTransmitter,     Config::OPTIONAL, "", "choose transmitter view, e.g. '***G18'");
    readConfig(config, "inputfileResiduals",    fileNameResiduals,   Config::MUSTSET,  "", "GNSS receiver residuals");
    readConfig(config, "R",                     a,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",     f,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    // ============================

    Ellipsoid                        ellipsoid(a, f);
    std::vector<Vector3d>            points;
    std::vector<std::vector<Double>> values(types.size());
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
            // find type for the satellite system
            UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), satType));

            // azimuth and elevation
            if((epoch.obsType.at(idType+0) != (GnssType::AZIMUT    + GnssType::L1)) ||
               (epoch.obsType.at(idType+1) != (GnssType::ELEVATION + GnssType::L1)) ||
               (epoch.obsType.at(idType+2) != (GnssType::AZIMUT    + GnssType::L2)) ||
               (epoch.obsType.at(idType+3) != (GnssType::ELEVATION + GnssType::L2)))
              throw(Exception("azimuth and elevation expected"));

            Double azimuth   = epoch.observation.at(idObs+0);
            Double elevation = epoch.observation.at(idObs+1);
            if(!typeTransmitter.hasWildcard(GnssType::PRN)) // isTransmitter
            {
              azimuth   = epoch.observation.at(idObs+2);
              elevation = epoch.observation.at(idObs+3);
            }
            const Vector3d point = ellipsoid(Angle(azimuth), Angle(elevation), 0);

            idObs  += 4;  // skip azimuth and elevation
            idType += 4;

            Bool                  found = FALSE;
            std::vector<Double>   valuesPerPoint(types.size(), NAN_EXPR);
            std::vector<GnssType> typesTmp = types;
            while((idType<epoch.obsType.size()) && (idObs<epoch.observation.size()) && (epoch.obsType.at(idType) == satType))
            {
              const GnssType type  = epoch.obsType.at(idType++);
              const Double   value = epoch.observation.at(idObs++);
              const UInt     idx   = GnssType::index(typesTmp, type);
              if((idx != NULLINDEX) && (type == typeTransmitter) && value)
              {
                valuesPerPoint.at(idx) = value;
                typesTmp.at(idx) = GnssType(static_cast<UInt64>(-1));
                found = TRUE;
              }
            } // while()

            if(found)
            {
              points.push_back(point);
              for(UInt i=0; i<values.size(); i++)
                values.at(i).push_back(valuesPerPoint.at(i));
            }
          } // for(satType)
        } // for(epoch)
      } // for(arcNo)
    } // for(idFile)

    // ============================

    logStatus<<"save values to file <"<<fileNameGriddedData<<">"<<Log::endl;
    GriddedData griddedData(ellipsoid, points, std::vector<Double>(points.size(), 1), values);
    writeFileGriddedData(fileNameGriddedData, griddedData);
    MiscGriddedData::printStatistics(griddedData);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
