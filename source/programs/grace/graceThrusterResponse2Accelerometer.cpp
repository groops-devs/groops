/***********************************************/
/**
* @file graceThrusterResponse2Accelerometer.cpp
*
* @brief Add modeled thruster responses to accelerometer.
*
* @author Torsten Mayer-Guerr
* @date 2022-07-31
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Add modeled thruster responses to accelerometer data.
The epochs and durations are given in the \configFile{inputfileThruster}{instrument} (THRUSTER).

The \configFile{inputfileThrusterResponse}{matrix} is a $(6\times 3)$ matrix with
the linear accelerations in the SRF ($x, y, z$) in one line per pair:
\begin{enumerate}
\item Negative Yaw,
\item Positive Pitch,
\item Positive Yaw,
\item Negative Pitch,
\item Negative Roll,
\item Positive Roll.
\end{enumerate}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Add modeled thruster responses to accelerometer.
* @ingroup programsGroup */
class GraceThrusterResponse2Accelerometer
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceThrusterResponse2Accelerometer, SINGLEPROCESS, "Add modeled thruster responses to accelerometer", Grace, Instrument)

/***********************************************/

void GraceThrusterResponse2Accelerometer::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameIn;
    FileName fileNameThruster;
    FileName fileNameResponse;

    readConfig(config, "outputfileAccelerometer",   fileNameOut,      Config::MUSTSET, "", "ACCELEROMETER");
    readConfig(config, "inputfileAccelerometer",    fileNameIn,       Config::MUSTSET, "", "ACCELEROMETER");
    readConfig(config, "inputfileThruster",         fileNameThruster, Config::MUSTSET, "", "THRUSTER");
    readConfig(config, "inputfileThrusterResponse", fileNameResponse, Config::MUSTSET, "", "thruster model (matrix with one line per pair)");
    if(isCreateSchema(config)) return;

    logStatus<<"read accelerometer data <"<<fileNameIn<<">"<<Log::endl;
    AccelerometerArc  acc      = InstrumentFile::read(fileNameIn);
    std::vector<Time> timesAcc = acc.times();
    ThrusterArc       thruster = InstrumentFile::read(fileNameThruster);
    Matrix            responses;
    readFileMatrix(fileNameResponse, responses);

    logStatus<<"add modeled thruster responses"<<Log::endl;
    auto iter = timesAcc.begin();
    for(UInt i=0; i<thruster.size(); i++)
    {
      // find acc epoch after thruster start
      iter = std::upper_bound(iter, timesAcc.end(), thruster.at(i).time);
      if(iter == timesAcc.end())   break;
      if(iter == timesAcc.begin()) iter++;

      Vector thrusterDuration = 1e-3*thruster.at(i).data(); // duration for each thruster
      for(UInt k=0; k<responses.rows(); k++)
        if(thrusterDuration(k) > 0)
        {
          auto   accEpoch = acc.begin() + std::distance(timesAcc.begin(), iter);
          Double fraction = (accEpoch->time - thruster.at(i).time).seconds() / (accEpoch->time - (accEpoch-1)->time).seconds();  // fractional part at begin
          do
          {
            if((accEpoch->time - thruster.at(i).time).seconds() - thrusterDuration(k) > 0)
              fraction -= ((accEpoch->time - thruster.at(i).time).seconds() - thrusterDuration(k)) / (accEpoch->time - (accEpoch-1)->time).seconds(); // fractional part at end
            accEpoch->acceleration += fraction * Vector3d(responses.row(k).trans());
            fraction = 1;
          }
          while(((accEpoch++)->time - thruster.at(i).time).seconds() < thrusterDuration(k));
        }
    }

    logInfo<<"write accelerometer data to <"<<fileNameOut<<">"<<Log::endl;
    InstrumentFile::write(fileNameOut, acc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
