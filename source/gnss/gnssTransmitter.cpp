/***********************************************/
/**
* @file gnssTransmitter.cpp
*
* @brief GNSS transmitter.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2013-06-28
*
*/
/***********************************************/

#include "base/import.h"
#include "inputOutput/logging.h"
#include "gnss.h"
#include "gnssTransmitter.h"

/***********************************************/
/***********************************************/

Gnss::Transmitter::Transmitter() : _gnss(nullptr), _idTrans(NULLINDEX)
{
}

/***********************************************/

Gnss::Transmitter::~Transmitter()
{
}

/***********************************************/

Gnss &Gnss::Transmitter::gnss() const
{
  if(_gnss == nullptr)
    throw(Exception("Transmitter is not registered in Gnss class"));
  return *_gnss;
}

/***********************************************/

void Gnss::Transmitter::transmitTime(UInt idEpoch, const Time &timeRecv, const Vector3d &posRecv, Time &timeTrans, Vector3d &posTrans) const
{
  try
  {
    // 0. iteration
    posTrans = position(idEpoch, timeRecv-seconds2time(20200e3/LIGHT_VELOCITY));
    // iteration
    Vector3d posOld;
    for(UInt i = 0; i < 10; i++)
    {
      timeTrans = timeRecv - seconds2time((posTrans-posRecv).r()/LIGHT_VELOCITY);
      posOld    = posTrans;
      posTrans  = position(idEpoch, timeTrans);

      if((posTrans-posOld).r() <= 0.0001) // 0.1 mm
        return;
    }

    logWarning << name() << " " << gnss().times.at(idEpoch).dateTimeStr() <<  ": no convergence for transmit time: posDiff = " << (posTrans-posOld).r()*1000%"%.1f mm"s << Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
