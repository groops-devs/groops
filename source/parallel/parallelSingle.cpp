/***********************************************/
/**
* @file parallelSingle.cpp
*
* @brief Wrapper for Message Passing Interface (MPI).
* All functions are empty statements for single processor version.
*
* @author Torsten Mayer-Guerr
* @date 2004-11-13
*
*/
/***********************************************/

#include "base/import.h"
#include "base/gnssType.h"
#include "parallel/parallel.h"

/***** FUNCTIONS *******************************/

namespace Parallel
{
/***********************************************/

CommunicatorPtr init(int /*argc*/, char */*argv*/[]) {return nullptr;}
std::function<void(UInt type, const std::string &str)> addChannel(std::function<void(UInt type, const std::string &str)> recieve, CommunicatorPtr /*comm*/) {return recieve;}
CommunicatorPtr splitCommunicator(UInt /*color*/, UInt /*key*/, CommunicatorPtr /*comm*/) {return nullptr;}
CommunicatorPtr createCommunicator(std::vector<UInt> /*ranks*/, CommunicatorPtr /*comm*/) {return nullptr;}
CommunicatorPtr selfCommunicator() {return nullptr;}
UInt myRank(CommunicatorPtr /*comm*/) {return 0;}
UInt size(CommunicatorPtr /*comm*/)   {return 1;}
void barrier(CommunicatorPtr /*comm*/) {}
void broadCastExceptions(CommunicatorPtr comm, std::function<void(CommunicatorPtr)> func) {func(comm);}
void send(const Byte */*x*/, UInt /*size*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void receive  (Byte  */*x*/, UInt /*size*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void broadCast(Byte  */*x*/, UInt /*size*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void send(const UInt     &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void send(const Double   &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void send(const Bool     &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void send(const Angle    &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void send(const Time     &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void send(const GnssType &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void send(const Vector3d &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void send(const Vector   &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void send(const Matrix   &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void receive(UInt        &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void receive(Double      &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void receive(Bool        &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void receive(Angle       &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void receive(Time        &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void receive(GnssType    &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void receive(Vector3d    &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void receive(Vector      &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void receive(Matrix      &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void broadCast(UInt      &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void broadCast(Double    &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void broadCast(Bool      &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void broadCast(Angle     &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void broadCast(Time      &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void broadCast(GnssType  &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void broadCast(Vector3d  &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void broadCast(Vector    &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
template<> void broadCast(Matrix    &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void reduceSum(UInt     &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void reduceSum(Double   &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void reduceSum(Bool     &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void reduceSum(Matrix   &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void reduceMin(UInt     &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void reduceMin(Double   &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void reduceMax(UInt     &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}
void reduceMax(Double   &/*x*/, UInt /*process*/, CommunicatorPtr /*comm*/) {}

/***********************************************/

} // namespace Parallel
