/***********************************************/
/**
* @file gnssAttitudeInfoCreate.cpp
*
* @brief Creates attitude info file used by SimulateStarCameraGnss.
*
* @author Sebastian Strasser
* @date 2020-11-26
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Creates attitude info file (\file{Instrument(MISCVALUES)}{instrument})
used by \program{SimulateStarCameraGnss}. One or more \config{attitudeInfo}s can be specified.
They are valid from \config{timeStart} until the start of the subsequent \config{attitudeInfo}.
\config{maxManeuverTime} is used by \program{SimulateStarCameraGnss} to look
for ongoing orbit maneuvers before/after the given orbit that might affect the attitude at
the beginning or end of a given orbit.

\fig{!hb}{0.9}{gnssAttitudeModes}{fig:gnssAttitudeModes2}{Overview of attitude modes used by GNSS satellites}

Here is a list of GNSS satellite types for which the attitude behavior is known and their
respective attitude modes and required parameters:
\begin{itemize}
\item \textbf{GPS-II/IIA} [1]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: shadowMaxYawSteeringAndRecovery
  \item \config{noonMode}: catchUpYawSteering
  \item \config{maxYawRate}: 0.12~deg/s
  \item \config{yawBias}: 0.5~deg
  \item \config{maxManeuverTime}: 2~h
\end{itemize}
\item \textbf{GPS-IIR/IIR-M} [1]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: catchUpYawSteering
  \item \config{noonMode}: catchUpYawSteering
  \item \config{maxYawRate}: 0.2~deg/s
  \item \config{maxManeuverTime}: 30~min
\end{itemize}
\item \textbf{GPS-IIF} [2]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: shadowConstantYawSteering
  \item \config{noonMode}: catchUpYawSteering
  \item \config{maxYawRate}: 0.11~deg/s
  \item \config{yawBias}: -0.7~deg
  \item \config{maxManeuverTime}: 1.5~h
\end{itemize}
\item \textbf{GLO-M} [3]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: shadowMaxYawSteeringAndStop
  \item \config{noonMode}: centeredMaxYawSteering
  \item \config{maxYawRate}: 0.25~deg/s
  \item \config{noonBetaThreshold}: 2~deg
  \item \config{maxManeuverTime}: 1.5~h
\end{itemize}
\item \textbf{GAL-1} [4]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: smoothedYawSteering1
  \item \config{noonMode}: smoothedYawSteering1
  \item \config{maxManeuverTime}: 1.5~h
\end{itemize}
\item \textbf{GAL-2} [4]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: smoothedYawSteering2
  \item \config{noonMode}: smoothedYawSteering2
  \item \config{midnightBetaThreshold}: 4.1~deg
  \item \config{noonBetaThreshold}: 4.1~deg
  \item \config{activationThreshold}: 10~deg
  \item \config{maxManeuverTime}: 5656~s
\end{itemize}
\item \textbf{BDS-2G/3G} [5, 6]
\begin{itemize}
  \item \config{defaultMode}: orbitNormal
  \item \config{midnightMode}: orbitNormal
  \item \config{noonMode}: orbitNormal
\end{itemize}
\item \textbf{BDS-2I} [5]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: betaDependentOrbitNormal
  \item \config{noonMode}: betaDependentOrbitNormal
  \item \config{maxYawRate}: 0.085~deg/s
  \item \config{midnightBetaThreshold}: 4~deg
  \item \config{noonBetaThreshold}: 4~deg
  \item \config{activationThreshold}: 5~deg
  \item \config{maxManeuverTime}: 24~h
\end{itemize}
\item \textbf{BDS-2M} [5]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: betaDependentOrbitNormal
  \item \config{noonMode}: betaDependentOrbitNormal
  \item \config{maxYawRate}: 0.159~deg/s
  \item \config{midnightBetaThreshold}: 4~deg
  \item \config{noonBetaThreshold}: 4~deg
  \item \config{activationThreshold}: 5~deg
  \item \config{maxManeuverTime}: 13~h
\end{itemize}
\item \textbf{BDS-3I/3SI} [6]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: smoothedYawSteering2
  \item \config{noonMode}: smoothedYawSteering2
  \item \config{midnightBetaThreshold}: 3~deg
  \item \config{noonBetaThreshold}: 3~deg
  \item \config{activationThreshold}: 6~deg
  \item \config{maxManeuverTime}: 5740~s
\end{itemize}
\item \textbf{BDS-3M/3SM} [6]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: smoothedYawSteering2
  \item \config{noonMode}: smoothedYawSteering2
  \item \config{midnightBetaThreshold}: 3~deg
  \item \config{noonBetaThreshold}: 3~deg
  \item \config{activationThreshold}: 6~deg
  \item \config{maxManeuverTime}: 3090~s
\end{itemize}
\item \textbf{QZS-1} [7]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: betaDependentOrbitNormal
  \item \config{noonMode}: betaDependentOrbitNormal
  \item \config{maxYawRate}: 0.01~deg/s
  \item \config{yawBias}: 180~deg
  \item \config{midnightBetaThreshold}: 20~deg
  \item \config{noonBetaThreshold}: 20~deg
  \item \config{activationThreshold}: 18.5~deg
  \item \config{maxManeuverTime}: 24~h
\end{itemize}
\item \textbf{QZS-2G} [7]
\begin{itemize}
  \item \config{defaultMode}: orbitNormal
  \item \config{midnightMode}: orbitNormal
  \item \config{noonMode}: orbitNormal
  \item \config{yawBias}: 180~deg
\end{itemize}
\item \textbf{QZS-2I} [7]
\begin{itemize}
  \item \config{defaultMode}: nominalYawSteering
  \item \config{midnightMode}: centeredMaxYawSteering
  \item \config{noonMode}: centeredMaxYawSteering
  \item \config{maxYawRate}: 0.055~deg/s
  \item \config{midnightBetaThreshold}: 5~deg
  \item \config{noonBetaThreshold}: 5~deg
  \item \config{maxManeuverTime}: 1.5~h
\end{itemize}
\end{itemize}

Some specific satellites may deviate in their attitude behavior or parameters
(e.g. G013-G040, R713, C005, C015, C017, J001).

References for the attitude behavior information:
\begin{enumerate}
\item \href{https://doi.org/10.1007/s10291-008-0092-1}{Kouba (2009)}
\item \href{https://doi.org/10.1007/s10291-016-0562-9}{Kuang et al. (2017)}
\item \href{https://doi.org/10.1016/j.asr.2010.09.007}{Dilssner et al. (2011)}
\item \url{https://www.gsc-europa.eu/support-to-developers/galileo-satellite-metadata#3}
\item \href{https://doi.org/10.1007/s10291-018-0783-1}{Wang et al. (2018)}
\item \href{https://doi.org/10.1017/S0373463318000103}{Li et al. (2018)}
\item \url{https://qzss.go.jp/en/technical/qzssinfo/index.html}
\end{enumerate}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Creates attitude info file used by SimulateStarCameraGnss.
* @ingroup programsGroup */
class GnssAttitudeInfoCreate
{
public:
  enum AttitudeMode
  {
    NOMINAL_YAW_STEERING = 0,
    ORBIT_NORMAL = 1,                         // e.g. BDS-2G, BDS-3G, QZS-2G
    CATCH_UP_YAW_STEERING = 2,                // e.g. GPS-* (noon), GPS-IIR (midnight)
    SHADOW_MAX_YAW_STEERING_AND_RECOVERY = 3, // e.g. GPS-IIA (midnight)
    SHADOW_MAX_YAW_STEERING_AND_STOP = 4,     // e.g. GLO-M (midnight)
    SHADOW_CONSTANT_YAW_STEERING = 5,         // e.g. GPS-IIF (midnight)
    CENTERED_MAX_YAW_STEERING = 6,            // e.g. QZS-2I, GLO-M (noon)
    SMOOTHED_YAW_STEERING_1 = 7,              // e.g. GAL-1
    SMOOTHED_YAW_STEERING_2 = 8,              // e.g. GAL-2, BDS-3M, BDS-3I
    BETA_DEPENDENT_ORBIT_NORMAL = 9           // e.g. BDS-2M, BDS-2I, QZS-1
  };

  class AttitudeInfo
  {
  public:
    Time         timeStart;
    AttitudeMode defaultMode;
    AttitudeMode midnightMode;
    AttitudeMode noonMode;
    Double       maxYawRate;
    Double       yawBias;
    Double       midnightBetaThreshold;
    Double       noonBetaThreshold;
    Double       activationThreshold;
    Double       maxManeuverTime;
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssAttitudeInfoCreate, SINGLEPROCESS, "Creates attitude info file used by SimulateStarCameraGnss.", Gnss, Instrument)

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, GnssAttitudeInfoCreate::AttitudeInfo &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    auto readMode = [&](const std::string &name, const std::string &annotation)
    {
      GnssAttitudeInfoCreate::AttitudeMode mode = GnssAttitudeInfoCreate::AttitudeMode::NOMINAL_YAW_STEERING;
      std::string choice;
      readConfigChoice(config, name,  choice, Config::MUSTSET, "nominalYawSteering",  annotation);
      if(readConfigChoiceElement(config, "nominalYawSteering", choice, "yaw to keep solar panels aligned to Sun (e.g. most GNSS satellites outside eclipse)"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::NOMINAL_YAW_STEERING;
      if(readConfigChoiceElement(config, "orbitNormal", choice, "keep fixed yaw angle, for example point X-axis in flight direction (e.g. BDS-2G, BDS-3G, QZS-2G)"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::ORBIT_NORMAL;
      if(name != "defaultMode" && readConfigChoiceElement(config, "catchUpYawSteering", choice, "yaw at maximum yaw rate to catch up to nominal yaw angle (e.g. GPS-* (noon), GPS-IIR (midnight))"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::CATCH_UP_YAW_STEERING;
      if(name == "midnightMode" && readConfigChoiceElement(config, "shadowMaxYawSteeringAndRecovery", choice, "yaw at maximum yaw rate from shadow start to end, recover after shadow (e.g. GPS-IIA (midnight))"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::SHADOW_MAX_YAW_STEERING_AND_RECOVERY;
      if(name == "midnightMode" && readConfigChoiceElement(config, "shadowMaxYawSteeringAndStop", choice, "yaw at maximum yaw rate from shadow start until nominal yaw angle at shadow end is reached, then stop (e.g. GLO-M (midnight))"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::SHADOW_MAX_YAW_STEERING_AND_STOP;
      if(name == "midnightMode" && readConfigChoiceElement(config, "shadowConstantYawSteering", choice, "yaw at constant yaw rate from shadow start to end (e.g. GPS-IIF (midnight))"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::SHADOW_CONSTANT_YAW_STEERING;
      if(name != "defaultMode" && readConfigChoiceElement(config, "centeredMaxYawSteering", choice, "yaw at maximum yaw rate centered around noon/midnight (e.g. QZS-2I, GLO-M (noon))"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::CENTERED_MAX_YAW_STEERING;
      if(name != "defaultMode" && readConfigChoiceElement(config, "smoothedYawSteering1", choice, "yaw based on an auxiliary Sun vector for a smooth yaw maneuver (e.g. GAL-1)"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::SMOOTHED_YAW_STEERING_1;
      if(name != "defaultMode" && readConfigChoiceElement(config, "smoothedYawSteering2", choice, "yaw based on a modified yaw-steering law for a smooth yaw maneuver (e.g. GAL-2, BDS-3M, BDS-3I)"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::SMOOTHED_YAW_STEERING_2;
      if(name != "defaultMode" && readConfigChoiceElement(config, "betaDependentOrbitNormal", choice, "switch to orbit normal mode if below beta angle threshold (e.g. BDS-2M, BDS-2I, QZS-1)"))
        mode = GnssAttitudeInfoCreate::AttitudeMode::BETA_DEPENDENT_ORBIT_NORMAL;
      endChoice(config);

      return mode;
    };

    readConfig(config, "timeStart",             var.timeStart,             Config::MUSTSET, "0",   "");
    var.defaultMode  = readMode("defaultMode",  "default attitude mode");
    var.midnightMode = readMode("midnightMode", "attitude mode for maneuvers around orbit midnight");
    var.noonMode     = readMode("noonMode",     "attitude mode for maneuvers around orbit noon");
    readConfig(config, "maxYawRate",            var.maxYawRate,            Config::DEFAULT, "0",  "[degree/s] maximum yaw rate of the satellite");
    readConfig(config, "yawBias",               var.yawBias,               Config::DEFAULT, "0",  "[degree] yaw bias applied in satellite attitude control system");
    readConfig(config, "midnightBetaThreshold", var.midnightBetaThreshold, Config::DEFAULT, "0",  "[degree] limit midnight maneuver to this absolute angle of the Sun above/below the satellite orbital plane");
    readConfig(config, "noonBetaThreshold",     var.noonBetaThreshold,     Config::DEFAULT, "0",  "[degree] limit noon maneuver to this absolute angle of the Sun above/below the satellite orbital plane");
    readConfig(config, "activationThreshold",   var.activationThreshold,   Config::DEFAULT, "0",  "[degree] limit maneuver to this yaw/Earth-spacecraft-Sun angle (depending on mode)");
    readConfig(config, "maxManeuverTime",       var.maxManeuverTime,       Config::DEFAULT, "0",  "[s] maximum duration of maneuver or maximum maneuver lookup time before/after orbit start/end");
    endSequence(config);

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssAttitudeInfoCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut;
    std::vector<AttitudeInfo> attitudeInfos;

    readConfig(config, "outputfileAttitudeInfo", fileNameOut,   Config::MUSTSET,  "", "");
    readConfig(config, "attitudeInfo",           attitudeInfos, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    MiscValuesArc arc;
    for(auto &info : attitudeInfos)
    {
      MiscValuesEpoch epoch(9);
      epoch.time      = info.timeStart;
      epoch.values(0) = info.defaultMode;
      epoch.values(1) = info.midnightMode;
      epoch.values(2) = info.noonMode;
      epoch.values(3) = std::fabs(info.maxYawRate);
      epoch.values(4) = info.yawBias;
      epoch.values(5) = std::fabs(info.midnightBetaThreshold);
      epoch.values(6) = std::fabs(info.noonBetaThreshold);
      epoch.values(7) = std::fabs(info.activationThreshold);
      epoch.values(8) = std::fabs(info.maxManeuverTime);
      arc.push_back(epoch);
    }

    logStatus<<"write attitude info file <"<<fileNameOut<<">"<<Log::endl;
    writeFileInstrument(fileNameOut, arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
